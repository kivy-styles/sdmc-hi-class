#!/usr/bin/env python3
"""
C3 No-Slip fixed-point refinement for the released future SDMC tails.

The baseline future action preserves the accepted z>=0 branch and approaches
the inverse-square mature fixed point.  This audit uses the shared C3
refinement implemented in run_future_homogeneous_covariant_stitch_audit.py:

    g_required(sigma) = -F_,sigma / (2 Z_released(sigma)),

with the future-only switch

    W(u)=1-exp(-u)(1+u+u^2/2+u^3/6),
    u=lambda ln(sigma).

Because W and its first three derivatives vanish at sigma=1, the accepted
present G3 jet is preserved through third field derivative.  This is the
smoothness level required by the hi_class perturbation gravity functions.

The default audit performs two release/reconstruction iterations with
lambda=500 and validates the resulting fixed action with a fresh homogeneous
evolution, including exact Bellini-Sawicki D, c_s^2, No-Slip, G_eff and slip
diagnostics.
"""
from pathlib import Path
import json
import numpy as np

import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/future_noslip_fixedpoint_refinement.json")
LAMBDA_SWITCH=500.0
NITER=2
NMAX=10.0

def one_release(model):
    NN,yy,diag=aud.release_trajectory(model,Nmax=NMAX,npts=2001)
    ns=np.asarray([x["noslip_alpha_combo"] for x in diag])
    imax=int(np.argmax(np.abs(ns)))
    return NN,yy,diag,{
      "max_abs_alphaB_plus_2alphaM":float(abs(ns[imax])),
      "at_ln_a":float(NN[imax]),
      "N_at_final":float(diag[-1]["N_struct"]),
      "q_max":float(max(x["q"] for x in diag)),
      "D_min":float(min(x["D_full"] for x in diag)),
    }

def refine(spec):
    model=aud.build_candidate(spec)
    model["_base_action"]=model["action"]
    history=[]

    for iteration in range(NITER+1):
        NN,yy,diag,row=one_release(model)
        row["iteration"]=iteration
        history.append(row)
        if iteration==NITER:
            break
        model=aud.refine_noslip_model(
          model,yy,blend_rate=LAMBDA_SWITCH
        )

    validation=aud.evolve(model,Nmax=NMAX,npts=2001)
    return {
      "lambda_switch":LAMBDA_SWITCH,
      "iterations":NITER,
      "history":history,
      "validation":validation,
    }

results={}
for spec in aud.CANDIDATES:
    results[spec["name"]]=refine(spec)

out={
  "status":(
    "future-only C3 inverse reconstruction of the linear-G3 coefficient; "
    "the accepted z>=0 action jet is unchanged through third derivative at "
    "sigma=1. The refined profile is validated by a fresh homogeneous "
    "evolution with exact background Horndeski health diagnostics. Native "
    "future hi_class perturbation propagation remains a separate gate."
  ),
  "results":results,
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_NOSLIP_FIXEDPOINT_REFINEMENT")
for name,r in results.items():
    v=r["validation"]
    print("REFINED",name,json.dumps({
      "history":r["history"],
      "health":v["health_diagnostics"],
      "future":v["future_kinematics"],
      "final":v["final"],
    },sort_keys=True))
