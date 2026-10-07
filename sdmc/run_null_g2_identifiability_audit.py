#!/usr/bin/env python3
"""
Conditional statistical-identifiability audit for the exploratory trajectory-null G2 completion.

This does NOT adopt the completion into SDMC and does NOT determine epsilon.
It records the model-selection consequence only if:
  (i) the accepted late300 trajectory is exactly unchanged by epsilon,
 (ii) the likelihood being evaluated depends only on observables on that trajectory,
(iii) the epsilon prior is normalized on an allowed interval.
"""
import json, math
from pathlib import Path

out={
  "status":"conditional identifiability audit for exploratory trajectory-null G2 completion",
  "completion":"Delta G2 = epsilon*(X-Xstar(phi))^4/Xstar(phi)^3",
  "assumptions":[
    "Delta G2 and all mixed derivatives required through total order 3 vanish on X=Xstar(phi).",
    "The accepted late300 background and linear observables are therefore independent of epsilon within numerical reconstruction tolerance.",
    "The likelihood considered is conditional on the accepted trajectory and does not integrate over off-trajectory initial conditions.",
    "The prior pi(epsilon) is normalized over whatever physically allowed interval is adopted."
  ],
  "evidence_identity":{
    "equation":"Z_epsilon = integral L0*pi(epsilon) d epsilon = L0",
    "normalized_prior_integral":1.0,
    "occam_factor_from_exactly_flat_epsilon_direction":1.0,
    "delta_logZ_from_epsilon_direction":0.0
  },
  "local_information_geometry":{
    "score_dlogL_depsilon":0.0,
    "fisher_eigenvalue_epsilon":0.0,
    "hessian_rank_added_by_epsilon":0,
    "regular_BIC_klogn_count_for_this_direction":"not applicable under exact non-identifiability"
  },
  "interpretation":[
    "epsilon changes off-trajectory phase-space stability and is therefore a real completion coordinate, not a gauge label.",
    "But it is not identifiable from an accepted-trajectory likelihood if that likelihood is exactly epsilon-independent.",
    "Counting it as one ordinary fitted cosmological parameter in regular AIC/BIC would over-penalize that conditional likelihood.",
    "epsilon can become identifiable if the statistical model integrates over initial conditions, uses observables sensitive to off-trajectory evolution, or the trajectory-null property is broken at higher/nonlinear order."
  ],
  "qualification":[
    "This result is conditional and does not make epsilon an SDMC prediction.",
    "The null-G2 completion remains exploratory and is not adopted into the accepted SDMC action.",
    "Full model-family evidence still requires priors for any identifiable SDMC structural/action coordinates."
  ]
}
Path("output").mkdir(exist_ok=True)
Path("output/null_g2_identifiability_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NULL_G2_IDENTIFIABILITY_AUDIT",json.dumps(out,sort_keys=True))
