#!/usr/bin/env python3
from pathlib import Path
import json,re,math

PATCH=Path("sdmc/apply_sdmc_null_g2_stabilizer.py").read_text()

# Verify that the exploratory completion is implemented only through G2 and its
# derivatives.  This is the repository-level statement needed for the closure
# audit; it is not a claim that the coefficient epsilon is derived.
assigns=re.findall(r'pgf->(G[2345][A-Za-z0-9_]*)\s*\+=',PATCH)
if not assigns:
    raise RuntimeError("no additive Horndeski-function modifications found")
nong2=[x for x in assigns if not x.startswith("G2")]
if nong2:
    raise RuntimeError("null completion unexpectedly modifies non-G2 sector: "+repr(nong2))

# The implemented deformation is epsilon Y^4 / Xs^3, Y=X-Xs(phi).
# Any mixed derivative d_X^a d_phi^b with total order <=3 still contains at
# least one power of Y and therefore vanishes at Y=0.  Record this directly.
derivative_orders=[]
for a in range(4):
    for b in range(4-a):
        remaining_power=4-a-b
        derivative_orders.append({
            "dX_order":a,"dphi_order":b,"total_order":a+b,
            "minimum_remaining_Y_power":remaining_power,
            "vanishes_on_trajectory":remaining_power>=1
        })
if not all(r["vanishes_on_trajectory"] for r in derivative_orders):
    raise RuntimeError("quartic null condition failed")

out={
  "status":"trajectory-null G2 compatibility audit against accepted SDMC closures",
  "deformation":"Delta G2 = epsilon*(X-Xstar(phi))^4/Xstar(phi)^3",
  "repository_checks":{
    "modified_horndeski_objects":sorted(set(assigns)),
    "non_G2_modifications_detected":nong2,
    "trajectory_derivatives_through_total_order_3_vanish":True,
    "derivative_table":derivative_orders
  },
  "closure_algebra":{
    "exact_no_slip":{
      "accepted_relation":"X*G3_X = -G4_phi  (equivalently alpha_B=-2 alpha_M on the accepted trajectory)",
      "delta_G3_X":0.0,
      "delta_G4_phi":0.0,
      "delta_no_slip_relation":0.0,
      "result":"unchanged by a pure Delta G2 deformation"
    },
    "tensor_luminality":{
      "accepted_sector":"G4=F(phi)/2 with no X dependence and G5=0",
      "delta_G4_X":0.0,
      "delta_G4_XX":0.0,
      "delta_G5":0.0,
      "result":"unchanged by a pure Delta G2 deformation"
    },
    "accepted_trajectory_background_and_linear_tangent":{
      "result":"unchanged to the derivative order used by the reconstructed background/linear system, because the deformation and mixed derivatives through total order 3 vanish at X=Xstar(phi)"
    }
  },
  "scientific_status":{
    "compatible_with_tested_closures":True,
    "selected_by_those_closures":False,
    "epsilon_derived":False,
    "reason":"No-Slip and tensor-luminality conditions constrain the G3/G4/G5 sectors or the accepted trajectory; they do not fix an off-trajectory G2 coefficient that is exactly null on that trajectory.",
    "required_for_adoption":"an additional physical principle, microscopic derivation, or an explicitly off-trajectory empirical/initial-condition criterion that fixes or constrains epsilon"
  },
  "model_selection_note":{
    "accepted_trajectory_identifiability":"flat in epsilon under the already-audited accepted-trajectory likelihood",
    "regular_parameter_penalty":"not justified for that exactly flat conditional likelihood alone",
    "full_family_warning":"epsilon can become identifiable if initial conditions/off-trajectory or nonlinear observables are part of the statistical model"
  }
}
Path("output").mkdir(exist_ok=True)
Path("output/null_g2_closure_compatibility.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NULL_G2_CLOSURE_COMPATIBILITY",json.dumps(out,sort_keys=True))
