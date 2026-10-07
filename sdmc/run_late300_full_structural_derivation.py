#!/usr/bin/env python3
from pathlib import Path
import json,re,math
import numpy as np
from scipy.optimize import least_squares, linprog
from scipy.integrate import solve_ivp, cumulative_trapezoid

ROOT=Path("derivation_inputs")
OUT=Path("output/late300_full_structural_derivation")
OUT.mkdir(parents=True,exist_ok=True)
cand=json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]
bgp=next(ROOT.rglob("linear_cov_target_00_background.dat"))

def table(path):
    ls=Path(path).read_text().splitlines()
    hdr=[x for x in ls if x.startswith("#") and re.search(r"(?:^|\s)1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    if a.ndim==1:a=a[None,:]
    return names,{n:a[:,i] for i,n in enumerate(names)}

names,d=table(bgp)
z=np.asarray(d["z"],float); N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"],float)
o=np.argsort(N); z=z[o]; N=N[o]; H=H[o]
H0=float(np.interp(0.0,N,H)); E=H/H0

h=cand["H0"]/100.; ORPH=4.17998772e-5
Om=(cand["omega_b"]+cand["omega_cdm"])/(h*h)
Or=ORPH/(h*h)
Ox=1-(cand["omega_b"]+cand["omega_cdm"]+ORPH)/(h*h)
rm=Om*np.exp(-3*N); rr=Or*np.exp(-4*N)
Hflat=np.sqrt(rm+rr+Ox)
wb=(rr/3)/(rm+rr)
x_r=3*wb
hlog=np.gradient(np.log(H),N,edge_order=2)
q=-1-hlog

def findcol(*needles):
    for n in names:
        low=n.lower()
        if all(x.lower() in low for x in needles):
            return n
    return None

def exactcol(name):
    for n in names:
        if n.strip().lower()==name.lower():
            return n
    return None
Fcol=exactcol("M2_smg")
deltaFcol=exactcol("delta_M2_smg")
if Fcol is None and deltaFcol is not None:
    # handled below as 1+delta_M2
    Fcol=deltaFcol
Dcol=exactcol("kin (D)") or findcol("kin","D")
cscol=exactcol("c_s^2") or findcol("c_s")
acol=exactcol("M2_running_smg")

out={"status":"full late300 structural-coordinate derivation and non-derivability audit",
     "columns":{"F":Fcol,"D":Dcol,"cs2":cscol,"alphaM":acol}}

# 1. Exact present closure.
out["Omega_x0"]={"candidate":cand["Omega_x0"],"derived":Ox,
                 "abs_residual":abs(cand["Omega_x0"]-Ox),
                 "status":"derived exactly from flat present closure"}

# 2. Exact action-trajectory decomposition of F and D. This is not first principles;
# it tests whether the coordinates are independently encoded after the action is specified.
if Fcol or acol:
    if Fcol:
        F=np.asarray(d[Fcol],float)[o]
        if Fcol.strip().lower()=="delta_m2_smg":
            F=1.0+F
        F_source="direct M2 column"
    else:
        alpha=np.asarray(d[acol],float)[o]
        lnF=cumulative_trapezoid(alpha,N,initial=0.0)
        F=np.exp(lnF)
        F_source="integrated alpha_M with early F=1 normalization"
    mask=np.isfinite(F)&(F>0)&(z<200)
    def Fmodel(p):
        AF,zc,w=p
        Nc=-np.log1p(zc)
        S=.5*(1+np.tanh((N[mask]-Nc)/(2*w)))
        return np.exp(AF*S)
    fit=least_squares(lambda p:(np.log(Fmodel(p))-np.log(F[mask])),
                      [cand["A_F"],cand["z_c"],cand["width"]],
                      bounds=([1e-6,0.05,0.01],[1.0,100.,5.0]),
                      xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=5000)
    AFf,zcf,wf=fit.x
    Nc=-np.log1p(zcf); Sall=.5*(1+np.tanh((N-Nc)/(2*wf)))
    out["F_trajectory_decomposition"]={
      "derived_from_frozen_trajectory":{"A_F":float(AFf),"z_c":float(zcf),"width":float(wf)},
      "candidate":{"A_F":cand["A_F"],"z_c":cand["z_c"],"width":cand["width"]},
      "max_logF_residual":float(np.max(np.abs(np.log(Fmodel(fit.x))-np.log(F[mask])))),
      "F_source":F_source,
      "status":"exact/near-exact decomposition of the already selected frozen action; not a first-principles prediction"
    }
else:
    AFf,zcf,wf=cand["A_F"],cand["z_c"],cand["width"]
    Nc=-np.log1p(zcf); Sall=.5*(1+np.tanh((N-Nc)/(2*wf)))

if Dcol:
    D=np.asarray(d[Dcol],float)[o]
    md=np.isfinite(D)&(Sall>1e-8)&(Sall<1-1e-8)&(z<1e5)
    def Dmodel(p):
        floor,d0,powr=p
        return floor+d0*np.power(Sall[md],powr)
    fd=least_squares(lambda p:Dmodel(p)-D[md],
                     [cand["D_floor"],cand["D0"],cand["kinetic_power"]],
                     bounds=([1e-8,1e-8,0.05],[2,3,5]),
                     xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=5000)
    out["D_trajectory_decomposition"]={
      "derived_from_frozen_trajectory":{"D_floor":float(fd.x[0]),"D0":float(fd.x[1]),"power":float(fd.x[2])},
      "candidate":{"D_floor":cand["D_floor"],"D0":cand["D0"],"power":cand["kinetic_power"]},
      "max_abs_residual":float(np.max(np.abs(Dmodel(fd.x)-D[md]))),
      "status":"decomposition of the selected kinetic trajectory, not an origin calculation"
    }

# 3. Expansion-sector decomposition. Fixed family widths used in the late300 lineage.
dN=0.5; tauA=.25; tauB=1.5
def expansion_model(p,NN=N):
    lam,zt,A,B=p
    zz=np.exp(-NN)-1
    rm_=Om*np.exp(-3*NN); rr_=Or*np.exp(-4*NN)
    wb_=(rr_/3)/(rm_+rr_)
    f0=3*(1+wb_)/(lam*lam)
    Nt=-np.log1p(zt)
    W=.5*(1-np.tanh((NN-Nt)/dN))
    f=W*f0
    delta=A*zz*np.exp(-zz/tauA)-B*zz*zz*np.exp(-zz/tauB)
    return np.sqrt((rm_+rr_+Ox)*(1+delta)**2/(1-f))
mx=np.isfinite(E)&(z>=0)&(z<200)
fe=least_squares(lambda p:(np.log(expansion_model(p)[mx])-np.log(E[mx])),
                 [cand["lambda_e"],cand["z_t"],cand["A_late"],cand["B_late"]],
                 bounds=([5,.5,-.2,0],[100,100,.2,.2]),
                 xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=10000)
out["background_trajectory_decomposition"]={
  "fixed_family_shape":{"DeltaN":dN,"tau_A":tauA,"tau_B":tauB},
  "derived_from_frozen_H":{"lambda_e":float(fe.x[0]),"z_t":float(fe.x[1]),"A_late":float(fe.x[2]),"B_late":float(fe.x[3])},
  "candidate":{"lambda_e":cand["lambda_e"],"z_t":cand["z_t"],"A_late":cand["A_late"],"B_late":cand["B_late"]},
  "max_logH_residual":float(np.max(np.abs(np.log(expansion_model(fe.x)[mx])-np.log(E[mx])))),
  "jacobian_singular_values":[float(x) for x in np.linalg.svd(fe.jac,compute_uv=False)],
  "status":"decomposition of accepted H(N); demonstrates identifiability within the ansatz, not first-principles derivation"
}

# 4. Conditional A_late kinematic relation.
q0=float(np.interp(0,N,q))
qc0=0.5*Om+Or-Ox
Aq=q0-qc0
out["A_late_kinematic_closure"]={"q0_target":q0,"q0_flat_control":qc0,
  "derived":Aq,"candidate":cand["A_late"],"abs_residual":abs(Aq-cand["A_late"]),
  "status":"conditional: derived if present q0 is independently predicted"}

# 5. Early tracker fraction and lambda_e extracted directly from accepted H.
# At high z delta_H is negligible. Invert f = 1-Hflat^2/H^2 and lambda^2=3(1+wb)/f
early=(z>100)&(z<1e5)&np.isfinite(E)
fobs=1-(Hflat[early]**2)/(E[early]**2)
lams=np.sqrt(np.maximum(3*(1+wb[early])/np.maximum(fobs,1e-30),0))
good=np.isfinite(lams)&(fobs>1e-6)&(fobs<.2)
out["lambda_tracker_inversion"]={
  "candidate":cand["lambda_e"],
  "median_from_accepted_early_background":float(np.median(lams[good])) if np.any(good) else None,
  "p16_p84":[float(np.percentile(lams[good],16)),float(np.percentile(lams[good],84))] if np.any(good) else None,
  "candidate_f_r":4/cand["lambda_e"]**2,"candidate_f_m":3/cand["lambda_e"]**2,
  "old_manuscript_f_r_anchor_0p01_implies_lambda":20.0,
  "status":"lambda follows algebraically once the early tracker fraction is predicted; current late300 fraction is not the old 1% anchor"
}

# 6. Recover tracker window W and infer z_t, DeltaN from accepted background after removing late two-lobe term.
zz=z
delta=cand["A_late"]*zz*np.exp(-zz/tauA)-cand["B_late"]*zz*zz*np.exp(-zz/tauB)
flate=(rm+rr+Ox)*(1+delta)**2
f_infer=1-flate/(E*E)
f0=3*(1+wb)/(cand["lambda_e"]**2)
Winfer=f_infer/f0
mw=(Winfer>1e-5)&(Winfer<1-1e-5)&np.isfinite(Winfer)&(z>1)&(z<100)
Y=np.arctanh(1-2*Winfer[mw])
X=N[mw]
if len(X)>5:
    co=np.polyfit(X,Y,1); slope,inter=co
    dNfit=1/slope; Ntfit=-inter/slope; ztfit=np.exp(-Ntfit)-1
    out["handoff_trajectory_inversion"]={
      "derived_from_accepted_H":{"z_t":float(ztfit),"DeltaN":float(dNfit)},
      "candidate_z_t":cand["z_t"],"fixed_DeltaN":dN,
      "status":"exact handoff extraction from accepted background, not a dynamical origin law"
    }

# 7. Minimum positive-energy F envelope and best smooth-window projection.
Rmr=(rm+rr)/(E*E)
idx=np.where(Rmr>1)[0]
if len(idx)>5:
    i0,i1=idx[0],idx[-1]
    Ni=N[i0:i1+1]; Ri=Rmr[i0:i1+1]
    sol=solve_ivp(lambda x,y:[float(np.interp(x,Ni,Ri))-y[0]],(Ni[0],Ni[-1]),[1.0],
                  dense_output=True,rtol=1e-11,atol=1e-13,max_step=.005)
    Fenv=sol.sol(Ni)[0]
    # fit monotonic logistic window to lnF envelope
    def lm(p):
        af,nc,w=p
        S=.5*(1+np.tanh((Ni-nc)/(2*w)))
        return af*S
    lf=least_squares(lambda p:lm(p)-np.log(Fenv),[.006,np.median(Ni),.5],
                     bounds=([1e-8,Ni[0]-5,.01],[.5,Ni[-1]+5,5]),max_nfev=10000)
    afm,ncm,wm=lf.x
    out["minimum_health_F_projection"]={
      "Fmax":float(np.max(Fenv)),"lnFmax":float(np.max(np.log(Fenv))),
      "projected":{"A_F":float(afm),"z_c":float(np.exp(-ncm)-1),"width":float(wm)},
      "accepted":{"A_F":cand["A_F"],"z_c":cand["z_c"],"width":cand["width"]},
      "rms_logF":float(np.sqrt(np.mean((lm(lf.x)-np.log(Fenv))**2))),
      "status":"health condition gives only a lower envelope; accepted F has additional physical/observational amplitude"
    }

# 8. Minimum kinetic support within D=floor+D0*S^p from native numerator Ns=D cs2.
if Dcol and cscol:
    cs=np.asarray(d[cscol],float)[o]
    Ns=D*cs
    mk=np.isfinite(Ns)&np.isfinite(Sall)&(Sall>=0)&(Sall<=1)&(Ns>0)
    X=Sall[mk]; Y=Ns[mk]
    best=None
    for powr in np.linspace(.1,4,157):
        xp=X**powr
        # floor + d0*x >= Y, nonnegative. Minimize mean support.
        A=-np.column_stack([np.ones_like(xp),xp]); b=-Y
        c=np.array([1.0,np.mean(xp)])
        lp=linprog(c,A_ub=A,b_ub=b,bounds=[(0,None),(0,None)],method="highs")
        if lp.success:
            floor,d0=lp.x
            val=lp.fun
            if best is None or val<best[0]: best=(val,powr,floor,d0)
    if best:
        _,powr,floor,d0=best
        support=floor+d0*X**powr
        out["minimum_kinetic_support"]={
          "minimum_mean_support_profile":{"D_floor":float(floor),"D0":float(d0),"power":float(powr),
                                          "max_constraint_slack":float(np.max(support-Y)),
                                          "min_constraint_slack":float(np.min(support-Y))},
          "accepted":{"D_floor":cand["D_floor"],"D0":cand["D0"],"power":cand["kinetic_power"]},
          "accepted_mean_D_over_sample":float(np.mean(cand["D_floor"]+cand["D0"]*X**cand["kinetic_power"])),
          "minimum_mean_D_over_sample":float(np.mean(support)),
          "status":"stability+subluminality alone do not select the accepted kinetic normalization"
        }

# 9. Minimal analytic power principle.
out["kinetic_power_minimality"]={
  "candidate_power":cand["kinetic_power"],
  "statement":"p=1 is the leading analytic linear response D(S)=D_floor+D0*S+O(S^2), but this is a model-definition/minimality principle rather than a numerical derivation of late300 because other powers were historically tested."
}

# 10. Coordinate verdict.
out["coordinate_verdict"]={
 "Omega_x0":"CLOSED: exact structural/background algebraic closure; no independent parameter.",
 "A_late":"CONDITIONAL: essentially fixed by q0-q_control, but q0 is not yet independently predicted by current SDMC.",
 "B_late":"CONDITIONAL: fixed once acoustic distance and late-profile widths are imposed; observational if the acoustic target is data supplied.",
 "lambda_e":"CONDITIONAL: lambda=sqrt(4/f_r)=sqrt(3/f_m); requires a first-principles prediction of the tracker fraction. The older f_r=0.01 anchor predicts 20, not late300's 18.40625.",
 "z_t":"OPEN AT THEORY LEVEL: extractable from H(N), but no current SDMC identity uniquely predicts the matter-to-structural handoff time.",
 "A_F_z_c_width":"OPEN AT THEORY LEVEL: positive energy and health constrain and localize the window but do not uniquely reproduce its accepted amplitude/timing/width.",
 "D_floor_D0":"OPEN AT THEORY LEVEL: native stability supplies lower bounds, not the accepted normalization.",
 "kinetic_power":"PARTLY REDUCIBLE BY MINIMALITY: p=1 is the leading analytic response, but was historically selected among alternatives and is not yet a microscopic theorem.",
 "DeltaN_tauA_tauB":"FIXED ANSATZ SHAPE choices in this family, not late300 fit coordinates; current manuscripts do not derive them uniquely."
}
out["model_selection_consequence"]={
 "fully_derived_sdmc_specific_coordinates_now":["Omega_x0"],
 "conditionally_reducible":["A_late","B_late","lambda_e","kinetic_power"],
 "still_require_new_dynamics_or_UV":["z_t","A_F","z_c","width","D_floor","D0"],
 "bottom_line":"The present manuscript/action identities are insufficient to derive all exact late300 coordinates. A zero full-family delta_k cannot yet be claimed. The correct next theoretical closure is a dynamical handoff plus microscopic/action law for F(S) and D(S)."
}
(OUT/"late300_full_structural_derivation.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_FULL_STRUCTURAL_DERIVATION",json.dumps(out,sort_keys=True))
