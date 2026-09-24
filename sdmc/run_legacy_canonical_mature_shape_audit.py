#!/usr/bin/env python3
"""
Audit the legacy canonical mature shape implied by the earlier SDMC
lambda_l=sqrt(2) endpoint against the released structural-field action.

For the mature quadratic normal form
    G2 = sigma^-2 [kappa1 Z + kappa2 Z^2 - U0]
with
    kappa1 = 2 F_inf (1-r),
    kappa2 = r F_inf/Zstar,
    U0 = F_inf Zstar (4-r),

the canonical limit is r=0:
    kappa2=0,
    kappa1=2F_inf,
    U0=4F_inf Zstar.

With phi=sqrt(2F_inf) ln sigma this becomes a canonical scalar with an
exponential potential whose slope in effective-Planck units is sqrt(2).
Therefore r=0 is the exact structural-field representation of the earlier
Part-VI late canonical coasting endpoint.

This audit keeps the covariant Balanced-Identity normalization
chi_inf=F_inf, hence N_inf=sqrt(8*pi/3), and scans only the future G2 stitch
rate muQ.  The accepted z>=0 action is unchanged.
"""
from pathlib import Path
import json,math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/legacy_canonical_mature_shape_audit.json")

R_CANONICAL=0.0
MU_GRID=[1.0,2.0,3.753897593021392,5.0,10.0,20.0,40.0,80.0]
NMAX=6.0
NPTS=1201
NU0=[0.0,0.1,1.0,10.0]

def exact_cs2(diag,NN):
    H=np.asarray([d["H"] for d in diag])
    F=np.asarray([d["F"] for d in diag])
    aB=np.asarray([d["alphaB"] for d in diag])
    aM=np.asarray([d["alphaM"] for d in diag])
    D=np.asarray([d["D_full"] for d in diag])
    aBN=np.asarray([d["alphaB_N"] for d in diag])
    c2=np.empty_like(NN)
    for j,d in enumerate(diag):
        b=aB[j]; mm=aM[j]; Hn=H[j]; Fv=F[j]
        xp=d["rho_smg_class"]+d["p_smg_class"]
        matter=(d["rho_m_action"]+4.*d["rho_r_action"]/3.)/3.
        num=((2.-b)*(b+2.*mm)/2.
             +1.5*(2.-b)*xp/(Hn*Hn)
             -1.5*(2.-2.*Fv+b*Fv)*matter/(Hn*Hn*Fv)
             +aBN[j])
        c2[j]=num/D[j]
    return c2

def background(model):
    sol=solve_ivp(lambda ne,y:aud.rhs_diag(ne,y,model)[0],
                  (0.,NMAX),[1.,aud.v0],
                  rtol=3e-9,atol=[1e-11,1e-14],
                  max_step=.01,dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    NN=np.linspace(0.,NMAX,NPTS)
    yy=sol.sol(NN)
    diag=[aud.rhs_diag(float(ne),[float(s),float(v)],model)[1]
          for ne,(s,v) in zip(NN,yy.T)]
    H=np.asarray([d["H"] for d in diag])
    F=np.asarray([d["F"] for d in diag])
    aB=np.asarray([d["alphaB"] for d in diag])
    D=np.asarray([d["D_full"] for d in diag])
    c2=exact_cs2(diag,NN)
    Qs=2.*F*D/(2.-aB)**2
    h=CubicSpline(NN,np.log(H))(NN,1)
    dlnQs=CubicSpline(NN,np.log(Qs))(NN,1)
    friction=3.+h+dlnQs
    chi=np.asarray([d["chi_action"] for d in diag]) if "chi_action" in diag[0] else None
    return NN,H,c2,Qs,friction,D,diag,yy

def mode(NN,H,c2,friction,nu0):
    if nu0==0.:
        return {"nu0":0.0,"success":True,"max_abs_zeta":1.0,"final_abs_zeta":1.0}
    sH=CubicSpline(NN,np.log(H)); sc=CubicSpline(NN,c2); sf=CubicSpline(NN,friction)
    H0=H[0]
    def rhs(ne,y):
        hn=math.exp(float(sH(ne)))
        nu=nu0*H0/(math.exp(ne)*hn)
        return [y[1],-float(sf(ne))*y[1]-float(sc(ne))*nu*nu*y[0]]
    sol=solve_ivp(rhs,(NN[0],NN[-1]),[1.,0.],rtol=2e-8,atol=2e-10,max_step=.01)
    z=np.asarray(sol.y[0])
    return {
      "nu0":nu0,"success":bool(sol.success),
      "max_abs_zeta":float(np.max(np.abs(z))),
      "final_abs_zeta":float(abs(z[-1]))
    }

rows=[]
for mu in MU_GRID:
    spec={"name":f"legacy_canonical_mu_{mu:g}",
          "chi_inf":aud.FINF,"r":R_CANONICAL,"muQ":mu}
    try:
        model=aud.build_refined_noslip_candidate(spec,iterations=2,blend_rate=500.)
        result=aud.evolve(model,Nmax=10.,npts=2001)
        NN,H,c2,Qs,friction,D,diag,yy=background(model)
        # rhs_diag alone does not attach chi_action; reconstruct it exactly as
        # evolve() does from the positive structural source.
        d2=[aud.rhs_diag(float(ne),[float(s),float(v)],model)[1]
            for ne,(s,v) in zip(NN,yy.T)]
        rhoXa=np.asarray([x["rhoX_structural"] for x in d2])
        chi=aud.CHI0*(rhoXa/rhoXa[0])*yy[0]*yy[0]
        gamma=chi/np.asarray([x["F"] for x in d2])
        modes=[mode(NN,H,c2,friction,x) for x in NU0]
        row={
          "muQ":mu,"status":"ok",
          "min_D":float(np.min(D)),
          "min_cs2":float(np.min(c2)),
          "max_cs2":float(np.max(c2)),
          "min_Qs":float(np.min(Qs)),
          "min_friction":float(np.min(friction)),
          "Gamma_max":float(np.max(gamma)),
          "Gamma_final":float(gamma[-1]),
          "Gamma_transient_overshoot_over_final":
            float(np.max(gamma)/max(gamma[-1],1e-300)-1.),
          "N_final":float(result["final"]["N_struct"]),
          "p_final":float(result["final"]["p"]),
          "q_final":float(result["final"]["q"]),
          "q_max":float(result["future_kinematics"]["q_max"]),
          "max_abs_alphaB_plus_2alphaM":
            float(result["health_diagnostics"]["max_abs_alphaB_plus_2alphaM"]),
          "modes":modes,
          "all_envelope_modes_non_growing":
            all(m["success"] and m["max_abs_zeta"]<=1.000001 for m in modes),
          "healthy":
            bool(np.min(D)>0. and np.min(c2)>0. and np.all(np.isfinite(c2))),
        }
    except Exception as e:
        row={"muQ":mu,"status":"failure","error":str(e),"healthy":False}
    rows.append(row)

healthy=[r for r in rows if r.get("healthy")]
best=min(healthy,key=lambda r:(r.get("Gamma_transient_overshoot_over_final",1e99),
                              abs(r.get("q_max",1e99)))) if healthy else None

out={
  "status":(
    "legacy canonical late-shape audit. r=0 is the exact canonical "
    "lambda_l=sqrt(2) mature normal form; only the unobserved future stitch "
    "rate is scanned and the accepted z>=0 action is unchanged."
  ),
  "analytic_equivalence":{
    "r":0.0,
    "kappa1":"2 F_inf",
    "kappa2":"0",
    "U0":"4 F_inf Zstar",
    "canonical_field":"phi=sqrt(2 F_inf) ln(sigma)",
    "effective_planck_field":"phi_hat=phi/sqrt(F_inf)=sqrt(2) ln(sigma)",
    "potential":"V proportional to exp(-sqrt(2)*phi_hat)",
    "lambda_l":"sqrt(2)",
    "cs2_inf":1.0,
    "D_inf":2.0,
    "N_inf":aud.XI
  },
  "rows":rows,
  "healthy_muQ":[r["muQ"] for r in healthy],
  "best_scanned_healthy":best
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LEGACY_CANONICAL_MATURE_SHAPE_AUDIT")
for r in rows:
    print("CANONICAL_SCAN",json.dumps(r,sort_keys=True))
print("BEST_HEALTHY",json.dumps(best,sort_keys=True))
