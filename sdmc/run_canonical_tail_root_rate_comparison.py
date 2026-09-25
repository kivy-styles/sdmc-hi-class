#!/usr/bin/env python3
"""
Compare the infrared tangent/transverse decay rates of the two positive
canonical Planck/braiding C3 tail roots.

The matter-loaded tangent mode is physical and should retain exponent ~1 on
both branches.  The transverse causal-shape relaxation can retain memory of
the chosen finite Planck/braiding interpolation.  This audit measures whether
the slow mu_F root changes the approach exponent while leaving the mature
canonical fixed point unchanged.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_tail_root_rate_comparison.json")
ROOTS=[
  ("slow",2.8610674139117047),
  ("fast",5.573218997214239),
]
BASE={
  "chi_inf":aud.FINF,
  "r":0.0,
  "muQ":50.0,
  "mu_k1":50.0,
  "mu_k2":5.0,
  "mu_V":120.0,
}

def fit_rate(N,arr,n1,n2,floor=1e-13):
    m=(N>=n1)&(N<=n2)&np.isfinite(arr)&(np.abs(arr)>floor)
    if np.count_nonzero(m)<20:
        return {"window":[n1,n2],"status":"insufficient_points",
                "npoints":int(np.count_nonzero(m))}
    x=N[m]; y=np.log(np.abs(arr[m]))
    slope,intercept=np.polyfit(x,y,1)
    pred=intercept+slope*x
    ssr=float(np.sum((y-pred)**2))
    sst=float(np.sum((y-np.mean(y))**2))
    return {
      "window":[n1,n2],
      "npoints":int(np.count_nonzero(m)),
      "decay_exponent":float(-slope),
      "R2":float(1.-ssr/sst if sst>0 else 1.),
      "value_start":float(abs(arr[m][0])),
      "value_end":float(abs(arr[m][-1])),
    }

def branch(label,muF):
    spec={**BASE,"name":f"canonical_rate_{label}","muF_override":muF}
    model=aud.build_refined_noslip_candidate(spec,iterations=2,blend_rate=500.)
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
    cal=(N>=5.)&(N<=10.)
    Ccal=float(np.median(Gamma[cal]/Gamma_ideal[cal]))
    epar=1.-Gamma/Ccal
    eperp=rcone

    wins=[(1.,2.),(1.5,2.5),(2.,3.),(2.5,3.5),(3.,4.),
          (3.5,4.5),(4.,5.),(4.5,5.5),(5.,6.)]
    return {
      "label":label,
      "mu_F":muF,
      "mu_g":muF+1.,
      "C_cal":Ccal,
      "tangent_fits":[fit_rate(N,epar,*w,floor=1e-12) for w in wins],
      "transverse_fits":[fit_rate(N,eperp,*w,floor=1e-13) for w in wins],
      "final":{
        "Gamma":float(Gamma[-1]),
        "r_cone":float(rcone[-1]),
        "cs2":float(cs2[-1]),
        "N_struct":float(diag[-1]["N_struct"]),
        "p":float(diag[-1]["p"]),
        "q":float(diag[-1]["q"]),
      }
    }

rows=[branch(*x) for x in ROOTS]
out={
  "status":"two-positive-root canonical tangent/transverse decay comparison",
  "branches":rows,
  "analytic_reference":{
    "physical_tangent_exponent":1.0,
    "k2_relative_tail_exponent":5.0,
    "slow_muF":ROOTS[0][1],
    "fast_muF":ROOTS[1][1],
    "interpretation":(
      "The tangent exponent is fixed by residual matter. The transverse "
      "approach can depend on whichever noncanonical coefficient correction "
      "is slowest after cancellations imposed by No-Slip."
    )
  }
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("CANONICAL_TAIL_ROOT_RATE_COMPARISON")
for row in rows:
    print("ROOT_RATE",row["label"],json.dumps({
      "mu_F":row["mu_F"],"mu_g":row["mu_g"],"final":row["final"]
    },sort_keys=True))
    for x in row["tangent_fits"]:
        print("TANGENT",row["label"],json.dumps(x,sort_keys=True))
    for x in row["transverse_fits"]:
        print("TRANSVERSE",row["label"],json.dumps(x,sort_keys=True))
