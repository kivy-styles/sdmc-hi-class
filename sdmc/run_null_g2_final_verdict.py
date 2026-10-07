#!/usr/bin/env python3
from pathlib import Path
import json

ROOT=Path("inputs")
OUT=Path("output/null_g2_final_verdict")
OUT.mkdir(parents=True,exist_ok=True)

def load_one(d):
    hits=list((ROOT/d).rglob("*.json"))
    if len(hits)!=1:
        raise RuntimeError(f"expected one json in {d}, got {hits}")
    return json.loads(hits[0].read_text())

basin=load_one("basin")
criterion=load_one("criterion")
field=load_one("field")
clock=load_one("clock")
ir=load_one("ir")
mature=load_one("mature")
model=load_one("model")

base=next(r for r in basin["rows"] if abs(r["epsilon"])<1e-15)
best_grad=max(basin["rows"],key=lambda r:r["gradient_stable_fraction_total"])
best_sub=max(basin["rows"],key=lambda r:r["subluminal_fraction_total"])

out={
 "status":"final scientific synthesis of exploratory trajectory-null G2 stabilizer",
 "accepted_action_baseline":{
   "epsilon":0.0,
   "n_total":base["n_total"],
   "n_gradient_stable":base["n_gradient_stable"],
   "n_subluminal_admissible":base["n_subluminal_admissible"],
   "n_superluminal":base["n_superluminal"],
   "n_integration_failure":base["n_integration_failure"],
   "interpretation":"The accepted reconstructed action already possesses a broad finite initial-condition basin; the exploratory completion is not required to reproduce late300 or its accepted-trajectory observables."
 },
 "finite_basin_tradeoff":{
   "best_sampled_gradient_coverage":{
     "epsilon":best_grad["epsilon"],
     "gradient_stable_fraction":best_grad["gradient_stable_fraction_total"],
     "subluminal_fraction":best_grad["subluminal_fraction_total"],
     "n_superluminal":best_grad["n_superluminal"],
     "n_integration_failure":best_grad["n_integration_failure"]
   },
   "best_sampled_subluminal_coverage":{
     "epsilon":best_sub["epsilon"],
     "gradient_stable_fraction":best_sub["gradient_stable_fraction_total"],
     "subluminal_fraction":best_sub["subluminal_fraction_total"],
     "n_superluminal":best_sub["n_superluminal"],
     "n_integration_failure":best_sub["n_integration_failure"]
   },
   "strict_dominance_found":False,
   "reason":"No sampled epsilon simultaneously improves every tested physical criterion over the unmodified accepted action; gradient-stability gains trade against superluminality and/or integration failure."
 },
 "accepted_trajectory_statistics":{
   "epsilon_identifiable":False,
   "score_dlogL_depsilon":criterion["identifiability"]["score_dlogL_depsilon"],
   "fisher_eigenvalue":criterion["identifiability"]["fisher_eigenvalue_epsilon"],
   "hessian_rank_added":criterion["identifiability"]["hessian_rank_added_by_epsilon"],
   "delta_logZ_from_epsilon":criterion["identifiability"]["delta_logZ_from_epsilon_direction"],
   "regular_AIC_BIC_charge":"not applicable for the exactly flat accepted-trajectory direction"
 },
 "field_redefinition":{
   "constant_epsilon_coordinate_invariant":field["interpretation"]["constant_epsilon_field_coordinate_invariant"],
   "structural_clock_resolution":clock["verdict"],
   "remaining_normalization_ambiguity":clock["verdict"]["normalization_ambiguity"]
 },
 "infrared_and_mature_limits":{
   "trajectory_ir_decay":ir["quartic_member"],
   "mature_robust_luminality":mature["verdict"],
   "combined_interpretation":"The term can be dynamically irrelevant along the mature attractor yet a constant nonzero coefficient is incompatible with the exact mature robust-luminality/minimal-derivative branch in an open off-trajectory neighborhood."
 },
 "model_selection":{
   "conditional_null_g2_occam_factor":model["null_g2_completion"]["conditional_occam_factor"],
   "conditional_delta_logZ":model["null_g2_completion"]["conditional_delta_logZ_from_epsilon"],
   "remaining_theory_family_issue":model["theory_family_complexity"]["remaining_issue"]
 },
 "final_verdict":{
   "adopt_into_accepted_late300_action":False,
   "accepted_action_status":"unchanged",
   "classification":"exploratory off-trajectory completion / diagnostic only",
   "reasons":[
     "late300 free evolution and accepted-trajectory observables do not require the term",
     "no tested epsilon strictly dominates the unmodified action across gradient stability, subluminality and integration success",
     "epsilon is exactly non-identifiable in the accepted-trajectory conditional likelihood",
     "a constant numerical epsilon is not generic-field-coordinate invariant unless normalized to the manuscript structural clock",
     "even in structural-clock normalization the overall eps_* is not derived",
     "constant nonzero eps_* is incompatible with exact mature robust luminality/minimal-derivative canonicalization"
   ],
   "scientific_use":"retain the null-G2 construction as proof that off-trajectory Horndeski completion freedom can reshape finite initial-condition stability without altering the accepted late300 trajectory",
   "next_theory_problem":"derive the microscopic transverse flow of the accepted action itself; do not use the null-G2 stabilizer as a substitute for that missing derivation"
 }
}

(OUT/"null_g2_final_verdict.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NULL_G2_FINAL_VERDICT",json.dumps(out,sort_keys=True))
