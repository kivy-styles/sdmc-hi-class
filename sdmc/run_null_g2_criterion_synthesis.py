#!/usr/bin/env python3
from pathlib import Path
import json

ROOT=Path("criterion_inputs")
OUT=Path("output/null_g2_criterion_synthesis")
OUT.mkdir(parents=True,exist_ok=True)

def load(name):
    hits=list(ROOT.rglob(name))
    if len(hits)!=1:
        raise RuntimeError(f"expected one {name}, got {hits}")
    return json.loads(hits[0].read_text())

window=load("null_g2_stabilizer_scan.json")
basin=load("null_g2_basin_classification.json")
ident=load("null_g2_identifiability_audit.json")

rows=basin["rows"]
best_grad=max(rows,key=lambda r:(r["gradient_stable_fraction_total"],r["subluminal_fraction_total"]))
best_sub=max(rows,key=lambda r:(r["subluminal_fraction_total"],r["gradient_stable_fraction_total"],-r["n_integration_failure"]))

out={
  "status":"trajectory-null G2 criterion-dependence synthesis",
  "challenge_grid_subluminal_window":{
    "smallest_strict_passing_epsilon":window["smallest_strict_passing_epsilon"],
    "largest_strict_passing_epsilon":window["largest_strict_passing_epsilon"],
    "definition":"selected challenge points require D>0, c_s^2>0, c_s^2<=1, stable base, and trajectory invariance tolerance",
    "scope":"challenge-grid window, not universal physical parameter interval"
  },
  "finite_basin_classification":{
    "grid":basin["grid"],
    "first_sampled_all_ghost_gradient_stable":basin["first_sampled_epsilon_all_20_ghost_gradient_stable"],
    "best_sampled_gradient_coverage":best_grad,
    "best_sampled_subluminal_coverage":best_sub,
    "important_distinction":"gradient stability and subluminality are separate criteria; maximizing one does not uniquely maximize the other"
  },
  "identifiability":{
    "score_dlogL_depsilon":ident["local_information_geometry"]["score_dlogL_depsilon"],
    "fisher_eigenvalue_epsilon":ident["local_information_geometry"]["fisher_eigenvalue_epsilon"],
    "hessian_rank_added_by_epsilon":ident["local_information_geometry"]["hessian_rank_added_by_epsilon"],
    "delta_logZ_from_epsilon_direction":ident["evidence_identity"]["delta_logZ_from_epsilon_direction"]
  },
  "verdict":{
    "unique_epsilon_selected":False,
    "reason":"Different finite-amplitude criteria prefer different sampled epsilon values, while the accepted-trajectory likelihood is exactly non-identifying in epsilon.",
    "allowed_statement":"The null-G2 family contains finite epsilon ranges that enlarge the stable phase-space basin while leaving the accepted trajectory unchanged.",
    "disallowed_statement":"Current tests do not derive a unique fundamental epsilon or justify adopting one as part of accepted SDMC.",
    "next_theory_requirement":"A new physical principle or microscopic completion must fix the off-trajectory G2 structure; alternatively an explicit initial-condition/off-trajectory statistical model can constrain epsilon empirically."
  }
}
(OUT/"null_g2_criterion_synthesis.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NULL_G2_CRITERION_SYNTHESIS",json.dumps(out,sort_keys=True))
