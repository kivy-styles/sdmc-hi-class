#!/usr/bin/env python3
"""
Audit the causal-cone shape variable associated with the SDMC scalar sector.

For the mature constant-F, G3=0 k-essence form,
    G2 = sigma^-2 f(Z),
the exact sound speed is
    c_s^2 = f_Z/(f_Z + 2 Z f_ZZ).

Define
    r_cone = (c_s^{-2}-1)/2 = Z f_ZZ/f_Z.

For the quadratic mature family, r_cone equals the shape modulus r exactly at
the fixed point.  Away from that limit (e.g. at the accepted present point,
where F and braiding still run), r_cone is only an equivalent causal-shape
diagnostic and must not be confused with the action coefficient r.
"""
from pathlib import Path
import json,math

SRC=Path("output/future_homogeneous_covariant_stitch_audit.json")
OUT=Path("output/causal_shape_modulus_audit.json")

d=json.loads(SRC.read_text())
acc=d["accepted_present"]

def rcone(cs2):
    return 0.5*(1./float(cs2)-1.)

present_cs2=float(acc["hi_class_saved_health"]["c_s^2"])
present_rcone=rcone(present_cs2)

candidates={}
for name,row in d["noslip_refined_candidates"].items():
    r=float(row["tail_parameters"]["r"])
    csinf=1./(1.+2.*r)
    csfinal=float(row["final"]["cs2_horndeski"])
    candidates[name]={
      "r_action":r,
      "cs2_exact_fixed_point":csinf,
      "r_from_exact_fixed_point_cs2":rcone(csinf),
      "cs2_at_ln_a_10":csfinal,
      "r_cone_at_ln_a_10":rcone(csfinal),
      "r_cone_minus_r_at_ln_a_10":rcone(csfinal)-r,
      "interpretation":(
        "In the mature constant-F, G3=0 limit r_cone is the actual quadratic "
        "shape modulus. The small finite-ln(a) mismatch is the residual "
        "matter/transition correction."
      )
    }

out={
  "status":(
    "causal-cone shape audit: r_cone=(c_s^{-2}-1)/2; exact shape modulus "
    "only in the mature constant-F, unbraided k-essence limit"
  ),
  "definition":{
    "r_cone":"(1/c_s^2 - 1)/2",
    "mature_identity":"r_cone=Z f_ZZ/f_Z",
    "quadratic_fixed_point":"r_cone=r",
    "flow_identity":"d ln(1+2 r_cone)/dN = - d ln(c_s^2)/dN",
  },
  "present":{
    "cs2":present_cs2,
    "r_cone_equivalent":present_rcone,
    "warning":(
      "The accepted present action still has alpha_M and braiding. Therefore "
      "3.114... is an equivalent causal-cone curvature, not the mature action "
      "shape modulus r."
    )
  },
  "candidates":candidates,
  "canonical":{
    "r":0.0,
    "cs2_inf":1.0,
    "r_cone_inf":0.0,
    "meaning":"robust luminality is equivalent to vanishing mature kinetic curvature"
  }
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("CAUSAL_SHAPE_MODULUS_AUDIT")
print("PRESENT",json.dumps(out["present"],sort_keys=True))
for name,row in candidates.items():
    print("CANDIDATE",name,json.dumps(row,sort_keys=True))
print("CANONICAL",json.dumps(out["canonical"],sort_keys=True))
