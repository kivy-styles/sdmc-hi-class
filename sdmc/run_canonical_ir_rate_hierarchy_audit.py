#!/usr/bin/env python3
"""
Measure the tangent/transverse decay hierarchy of the canonical infrared
closure manifold.

Tangent coordinate:
    eps_parallel = 1 - Gamma/C_cal
which is the residual matter-loaded normalization displacement.

Transverse coordinate:
    eps_perp = r_cone = (1/c_s^2 - 1)/2
which measures departure from the canonical luminal kinetic shape.

A normally attractive one-dimensional IR manifold should show the transverse
coordinate decaying parametrically faster than the tangent matter mode.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_ir_rate_hierarchy_audit.json")

SPEC={
  "name":"legacy_canonical_ir_rates",
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

cs2=np.empty_like(N)
for j,d in enumerate(diag):
    H=float(d["H"]); Fv=float(d["F"])
    b=float(d["alphaB"]); m=float(d["alphaM"])
    D=float(d["D_full"]); dB=float(d["alphaB_N"])
    hln=-1.-float(d["q"])
    matter=float(d["rho_m_action"]+4.*d["rho_r_action"]/3.)
    num=((2.-b)*(-hln+.5*b+m)-matter/(H*H*Fv)+dB)
    cs2[j]=num/D

rcone=.5*(1./cs2-1.)
Nrel=np.asarray([d["N_struct"] for d in diag])/aud.XI
Gamma_ideal=(2.+Nrel*Nrel)/3.

late=(N>=5.)&(N<=10.)
Ccal=float(np.median(Gamma[late]/Gamma_ideal[late]))
Gamma_hat=Gamma/Ccal

eps_parallel=np.abs(1.-Gamma_hat)
eps_perp=np.abs(rcone)

def fit_rate(arr,n1,n2,floor=1e-14):
    m=(N>=n1)&(N<=n2)&np.isfinite(arr)&(arr>floor)
    if np.sum(m)<20:
        return {"window":[n1,n2],"status":"insufficient_points",
                "npoints":int(np.sum(m))}
    x=N[m]; y=np.log(arr[m])
    slope,intercept=np.polyfit(x,y,1)
    pred=intercept+slope*x
    ssr=float(np.sum((y-pred)**2))
    sst=float(np.sum((y-np.mean(y))**2))
    r2=1.-ssr/sst if sst>0 else 1.
    return {
      "window":[n1,n2],
      "npoints":int(np.sum(m)),
      "decay_exponent":float(-slope),
      "intercept":float(intercept),
      "R2":float(r2),
      "value_start":float(arr[m][0]),
      "value_end":float(arr[m][-1]),
    }

windows=[(1.,2.),(1.5,2.5),(2.,3.),(2.5,3.5),(3.,4.),(3.5,4.5),(4.,5.)]
parallel_fits=[fit_rate(eps_parallel,*w,floor=1e-12) for w in windows]
perp_fits=[fit_rate(eps_perp,*w,floor=1e-13) for w in windows]

# Pointwise logarithmic decay rates where numerically reliable.
def local_rate(arr):
    out=np.full_like(arr,np.nan,dtype=float)
    good=arr>1e-13
    l=np.full_like(arr,np.nan,dtype=float)
    l[good]=np.log(arr[good])
    # compute on contiguous reliable interior using gradient of a clipped log
    idx=np.where(good)[0]
    if len(idx)>4:
        # reliability is contiguous over the physical interval here
        lo,hi=idx[0],idx[-1]
        out[lo:hi+1]=-np.gradient(l[lo:hi+1],N[lo:hi+1],edge_order=2)
    return out

nu_par=local_rate(eps_parallel)
nu_perp=local_rate(eps_perp)

samples=[]
for nt in [1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,6.]:
    i=int(np.argmin(np.abs(N-nt)))
    samples.append({
      "ln_a":float(N[i]),
      "Gamma_hat":float(Gamma_hat[i]),
      "eps_parallel":float(eps_parallel[i]),
      "r_cone":float(rcone[i]),
      "eps_perp":float(eps_perp[i]),
      "nu_parallel_local":None if not np.isfinite(nu_par[i]) else float(nu_par[i]),
      "nu_perp_local":None if not np.isfinite(nu_perp[i]) else float(nu_perp[i]),
      "rate_ratio_perp_over_parallel":
        None if not (np.isfinite(nu_par[i]) and np.isfinite(nu_perp[i]) and abs(nu_par[i])>1e-12)
        else float(nu_perp[i]/nu_par[i]),
    })

# Lyapunov onset result from the preceding audit: re-evaluate the simple
# lambda=1 distance and find the earliest grid point after which it stays
# decreasing through ln a=9.8.
u=np.log(Gamma)
v=-np.log(cs2)
L=.5*(u*u+v*v)
dL=np.gradient(L,N,edge_order=2)
end=np.where(N<=9.8)[0][-1]
tol=2e-8*max(1.,float(np.max(np.abs(dL[:end+1]))))
onset=None
for i in range(end+1):
    if np.all(dL[i:end+1] <= tol):
        onset=i
        break

out={
  "status":(
    "canonical IR rate hierarchy audit: compares residual matter-normalization "
    "decay tangent to the manifold with causal-shape decay transverse to it"
  ),
  "calibration":{"C_cal":Ccal},
  "coordinates":{
    "parallel":"abs(1-Gamma/C_cal)",
    "perpendicular":"abs(r_cone)",
  },
  "parallel_fits":parallel_fits,
  "perpendicular_fits":perp_fits,
  "samples":samples,
  "lyapunov_lambda1":{
    "definition":"L=1/2[(ln Gamma)^2+(-ln cs2)^2]",
    "earliest_monotone_grid_ln_a":
      None if onset is None else float(N[onset]),
    "earliest_monotone_scale_factor":
      None if onset is None else float(math.exp(N[onset])),
    "grid_step_ln_a":float(N[1]-N[0]),
    "tolerance":float(tol),
  },
  "analytic_expectation":{
    "parallel_IR":"eps_parallel ~ m1 a^-1 -> exponent 1",
    "perpendicular_transition":"higher operators decay much faster; k2 relative tail ~ a^-5 up to log-polynomial/mixing corrections",
    "interpretation":(
      "A transverse exponent substantially larger than one demonstrates "
      "rapid canonicalization onto a slowly evolving matter-loaded IR manifold."
    ),
  },
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("CANONICAL_IR_RATE_HIERARCHY_AUDIT")
print("LYAPUNOV_ONSET",json.dumps(out["lyapunov_lambda1"],sort_keys=True))
for r in parallel_fits:
    print("PARALLEL_FIT",json.dumps(r,sort_keys=True))
for r in perp_fits:
    print("PERP_FIT",json.dumps(r,sort_keys=True))
for r in samples:
    print("RATE_SAMPLE",json.dumps(r,sort_keys=True))
