#!/usr/bin/env python3
"""
Test whether the successful refined canonical future trajectory admits a
single positive quadratic Lyapunov distance in the two closure coordinates

    u = ln Gamma,
    v = ln(1+2 r_cone) = -ln(c_s^2).

For
    L_lambda = 1/2 [u^2 + lambda v^2],
we seek one constant lambda>0 such that dL_lambda/dN <= 0 along the full
canonical action-level trajectory.  This is a path-level existence test, not a
proof that the microscopic SDMC beta field is globally a gradient flow.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_closure_lyapunov_audit.json")

SPEC={
  "name":"legacy_canonical_lyapunov",
  "chi_inf":aud.FINF,
  "r":0.0,
  "muQ":50.0,
  "mu_k1":50.0,
  "mu_k2":5.0,
  "mu_V":120.0,
}

model=aud.build_refined_noslip_candidate(SPEC,iterations=2,blend_rate=500.)
NN,yy,diag=aud.release_trajectory(model,Nmax=10.,npts=4001)

N=np.asarray(NN)
rhoXa=np.asarray([d["rhoX_structural"] for d in diag])
F=np.asarray([d["F"] for d in diag])
chi=aud.CHI0*(rhoXa/rhoXa[0])*yy[0]*yy[0]
Gamma=chi/F

# Exact Bellini-Sawicki sound speed, same compact expression used by evolve().
cs2=np.empty_like(N)
for j,d in enumerate(diag):
    H=float(d["H"]); Fv=float(d["F"])
    b=float(d["alphaB"]); m=float(d["alphaM"])
    D=float(d["D_full"])
    dB=float(d["alphaB_N"])
    hln=-1.-float(d["q"])
    matter=float(d["rho_m_action"]+4.*d["rho_r_action"]/3.)
    num=((2.-b)*(-hln+.5*b+m)
         -matter/(H*H*Fv)
         +dB)
    cs2[j]=num/D

if np.any(Gamma<=0.) or np.any(cs2<=0.):
    raise RuntimeError("closure coordinates require Gamma>0 and cs2>0")

u=np.log(Gamma)
v=-np.log(cs2)  # exactly ln(1+2 r_cone)
rcone=.5*(1./cs2-1.)

du=np.gradient(u,N,edge_order=2)
dv=np.gradient(v,N,edge_order=2)
A=u*du
B=v*dv

# Need A + lambda B <= tol for a single lambda>0.
# Exclude only the last few points where both coordinates and derivatives are
# at roundoff scale; retain the entire physical transition otherwise.
mask=(N<=9.8)
scale=max(1.,float(np.max(np.abs(A[mask]))),float(np.max(np.abs(B[mask]))))
tol=2e-8*scale

lo=0.0
hi=float("inf")
incompatible=[]
constraints=[]
for i in np.where(mask)[0]:
    a=float(A[i]); b=float(B[i])
    if abs(b)<1e-14:
        if a>tol:
            incompatible.append({"ln_a":float(N[i]),"A":a,"B":b,
                                 "reason":"B~0 while normalization term increases"})
        continue
    bound=(-a+tol)/b
    if b>0.:
        hi=min(hi,bound)
        constraints.append(("upper",float(N[i]),bound))
    else:
        lo=max(lo,bound)
        constraints.append(("lower",float(N[i]),bound))

lo=max(lo,0.)
feasible=(not incompatible) and hi>0. and lo<=hi
if feasible:
    # Prefer lambda=1 if allowed; otherwise use a comfortably interior point.
    if lo<=1.<=hi:
        lam=1.0
    elif math.isinf(hi):
        lam=max(1.,1.2*lo)
    elif lo<=0.:
        lam=.5*hi
    else:
        lam=math.sqrt(lo*hi)
else:
    lam=1.0

L=.5*(u*u+lam*v*v)
dL=A+lam*B

def sign_changes(arr,eps=1e-10):
    s=np.sign(np.where(np.abs(arr)<eps,0.,arr))
    nz=s[s!=0.]
    return int(np.sum(nz[1:]*nz[:-1]<0.)) if len(nz)>1 else 0

def extrema(y):
    dy=np.gradient(y,N,edge_order=2)
    idx=np.where(dy[1:]*dy[:-1]<0.)[0]+1
    return [{"ln_a":float(N[i]),"value":float(y[i])} for i in idx[:20]]

samples=[]
for nt in [0.,.25,.5,.75,1.,1.5,2.,3.,5.,8.,10.]:
    i=int(np.argmin(np.abs(N-nt)))
    samples.append({
      "ln_a":float(N[i]),
      "Gamma":float(Gamma[i]),
      "r_cone":float(rcone[i]),
      "cs2":float(cs2[i]),
      "u_lnGamma":float(u[i]),
      "v_minus_lncs2":float(v[i]),
      "Lyapunov":float(L[i]),
      "dLyapunov_dN":float(dL[i]),
    })

# Test whether Gamma alone can serve as a global order parameter: if Gamma is
# non-monotone, the answer is no.  Also find near-equal Gamma pairs with very
# different r_cone to demonstrate path hysteresis in the closure plane.
pairs=[]
for i in range(0,len(N),8):
    dg=np.abs(Gamma-Gamma[i])
    cand=np.where((dg<2e-3)&(np.abs(N-N[i])>0.5))[0]
    if len(cand):
        j=cand[np.argmax(np.abs(rcone[cand]-rcone[i]))]
        dr=abs(float(rcone[j]-rcone[i]))
        if dr>0.05:
            pairs.append({
              "Gamma_1":float(Gamma[i]),"r_cone_1":float(rcone[i]),"ln_a_1":float(N[i]),
              "Gamma_2":float(Gamma[j]),"r_cone_2":float(rcone[j]),"ln_a_2":float(N[j]),
              "abs_delta_Gamma":float(abs(Gamma[j]-Gamma[i])),
              "abs_delta_r_cone":dr,
            })
pairs=sorted(pairs,key=lambda x:x["abs_delta_r_cone"],reverse=True)[:5]

out={
  "status":(
    "path-level Lyapunov audit of the successful refined canonical closure "
    "trajectory in (Gamma,r_cone); not a proof of a global microscopic "
    "gradient flow"
  ),
  "coordinates":{
    "Gamma":"chi/F",
    "r_cone":"(1/c_s^2-1)/2",
    "u":"ln Gamma",
    "v":"ln(1+2 r_cone)=-ln c_s^2",
  },
  "trajectory":{
    "Gamma_initial":float(Gamma[0]),
    "Gamma_final":float(Gamma[-1]),
    "Gamma_max":float(np.max(Gamma)),
    "Gamma_max_at_ln_a":float(N[int(np.argmax(Gamma))]),
    "r_cone_initial":float(rcone[0]),
    "r_cone_final":float(rcone[-1]),
    "r_cone_min":float(np.min(rcone)),
    "r_cone_min_at_ln_a":float(N[int(np.argmin(rcone))]),
    "Gamma_derivative_sign_changes":sign_changes(np.gradient(Gamma,N)),
    "r_cone_derivative_sign_changes":sign_changes(np.gradient(rcone,N)),
    "Gamma_extrema":extrema(Gamma),
    "r_cone_extrema":extrema(rcone),
    "same_Gamma_different_shape_examples":pairs,
  },
  "quadratic_lyapunov":{
    "definition":"L=1/2[(ln Gamma)^2 + lambda (-ln c_s^2)^2]",
    "tolerance":tol,
    "lambda_lower_bound":float(lo),
    "lambda_upper_bound":None if math.isinf(hi) else float(hi),
    "positive_lambda_interval_exists":bool(feasible),
    "selected_lambda":float(lam),
    "max_dL_dN_through_ln_a_9p8":float(np.max(dL[mask])),
    "fraction_points_dL_positive_beyond_tolerance":
      float(np.mean(dL[mask]>tol)),
    "L_initial":float(L[0]),
    "L_final":float(L[-1]),
    "L_drop_factor":float(L[-1]/L[0]) if L[0]>0 else None,
    "incompatible_constraints":incompatible[:20],
  },
  "samples":samples,
  "interpretation":{
    "Gamma_alone_is_global_order_parameter":
      bool(sign_changes(np.gradient(Gamma,N))==0),
    "r_cone_alone_is_global_order_parameter":
      bool(sign_changes(np.gradient(rcone,N))==0),
    "meaning":(
      "If the positive-lambda interval exists, the demonstrated canonical "
      "trajectory admits a single monotone quadratic distance even when its "
      "individual coordinates overshoot. This supports, but does not prove, "
      "a common relaxation potential."
    )
  }
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("CANONICAL_CLOSURE_LYAPUNOV_AUDIT")
print("TRAJECTORY",json.dumps(out["trajectory"],sort_keys=True))
print("LYAPUNOV",json.dumps(out["quadratic_lyapunov"],sort_keys=True))
for row in samples:
    print("SAMPLE",json.dumps(row,sort_keys=True))
