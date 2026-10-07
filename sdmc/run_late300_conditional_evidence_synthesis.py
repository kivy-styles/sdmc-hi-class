#!/usr/bin/env python3
from pathlib import Path
import json,math

def load_one(root,name):
    hits=list(Path(root).rglob(name))
    if len(hits)!=1:
        raise RuntimeError((root,name,[str(x) for x in hits]))
    return json.loads(hits[0].read_text())

p=load_one("inputs/planck","planck_conditional_bayes_factor.json")
hl=load_one("inputs/hsc_late","late300.json")
hc=load_one("inputs/hsc_local","local021.json")

dh=float(hl["logZ_laplace"])-float(hc["logZ_laplace"])
dp=float(p["delta_logZ_late300_minus_local021"])
joint=dp+dh
out={
  "status":"provisional fixed-cosmology conditional Planck+HSC evidence synthesis",
  "planck":{
    "delta_logZ_late300_minus_local021":dp,
    "BF_late300_over_local021":float(p["BF_late300_over_local021"]),
    "method":"high-precision nested nuisance evidence"
  },
  "hsc":{
    "delta_logZ_late300_minus_local021":dh,
    "BF_late300_over_local021":math.exp(dh),
    "late300_logZ_laplace":float(hl["logZ_laplace"]),
    "local021_logZ_laplace":float(hc["logZ_laplace"]),
    "method":"14D Laplace nuisance evidence around corrected-ReACT us_native MAP",
    "late300_hessian_positive_definite":bool(hl["hessian_positive_definite"]),
    "local021_hessian_positive_definite":bool(hc["hessian_positive_definite"])
  },
  "joint":{
    "delta_logZ_late300_minus_local021":joint,
    "BF_late300_over_local021":math.exp(joint)
  },
  "qualification":[
    "Planck and HSC nuisance sectors are independent conditional on the two fixed cosmologies, so their conditional log-evidence ratios may be added.",
    "The HSC factor is presently a Laplace approximation, not the still-running nested HSC evidence.",
    "No SDMC structural/action/cosmological prior-volume penalty is included; this is not full model-family evidence."
  ],
  "verdict":"At the present conditional level the high-precision Planck preference for local021 and the HSC Laplace preference for late300 nearly cancel. Replace the HSC term with nested evidence when that run completes."
}
Path("output").mkdir(exist_ok=True)
Path("output/late300_conditional_evidence_synthesis.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_CONDITIONAL_EVIDENCE_SYNTHESIS",json.dumps(out,sort_keys=True))
