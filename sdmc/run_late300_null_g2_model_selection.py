#!/usr/bin/env python3
from pathlib import Path
import json, math

ROOT=Path("inputs")
OUT=Path("output"); OUT.mkdir(exist_ok=True)

def find_json(root, needle):
    hits=[p for p in Path(root).rglob("*.json") if needle in p.name]
    if len(hits)!=1:
        raise RuntimeError(f"expected one {needle} json under {root}, got {hits}")
    return json.loads(hits[0].read_text())

ms=find_json(ROOT/"model_selection","late300_model_selection_highprecision")
ident=find_json(ROOT/"identifiability","null_g2_identifiability_audit")
cond=find_json(ROOT/"conditional","late300_conditional_evidence_synthesis")

if ident["local_information_geometry"]["hessian_rank_added_by_epsilon"] != 0:
    raise RuntimeError("null-G2 epsilon unexpectedly identifiable")
if abs(ident["evidence_identity"]["delta_logZ_from_epsilon_direction"]) > 1e-12:
    raise RuntimeError("null-G2 epsilon has nonzero accepted-trajectory evidence contribution")

out={
  "status":"late300 model-selection update with trajectory-null G2 identifiability",
  "accepted_action_replay":{
    "post_reconstruction_delta_k":0,
    "note":"No post-reconstruction structural retuning is performed in the frozen-action replay."
  },
  "null_g2_completion":{
    "completion":ident["completion"],
    "real_completion_coordinate":True,
    "accepted_trajectory_identifiable":False,
    "score_dlogL_depsilon":ident["local_information_geometry"]["score_dlogL_depsilon"],
    "fisher_eigenvalue_epsilon":ident["local_information_geometry"]["fisher_eigenvalue_epsilon"],
    "hessian_rank_added":ident["local_information_geometry"]["hessian_rank_added_by_epsilon"],
    "conditional_delta_logZ_from_epsilon":ident["evidence_identity"]["delta_logZ_from_epsilon_direction"],
    "conditional_occam_factor":ident["evidence_identity"]["occam_factor_from_exactly_flat_epsilon_direction"],
    "regular_AIC_BIC_parameter_charge":"not applicable for this exactly flat direction in the accepted-trajectory conditional likelihood",
    "qualification":"This does not make epsilon a prediction. It becomes statistically relevant if one integrates over off-trajectory initial conditions, uses observables that probe the off-trajectory completion, or breaks the exact trajectory-null property."
  },
  "fixed_cosmology_conditional_evidence":{
    "planck_delta_logZ":cond["planck"]["delta_logZ_late300_minus_local021"],
    "hsc_delta_logZ_laplace":cond["hsc"]["delta_logZ_late300_minus_local021"],
    "joint_delta_logZ":cond["joint"]["delta_logZ_late300_minus_local021"],
    "joint_BF_late300_over_local021":cond["joint"]["BF_late300_over_local021"],
    "joint_step_sensitivity_delta_logZ":cond["joint"]["step_sensitivity_delta_logZ"],
    "scope":"Conditional fixed-cosmology nuisance evidence; no SDMC structural/action/cosmological prior-volume integration."
  },
  "completed_hsc_raw_fit":ms["HSC"],
  "theory_family_complexity":{
    "historically_explored_union_count":ms["parameter_provenance"]["historically_explored_union_count"],
    "model_family_delta_k":ms["parameter_provenance"]["model_family_delta_k"],
    "remaining_issue":"The null-G2 flat direction should not be added as one regular BIC parameter in the accepted-trajectory conditional likelihood, but the effective dimensionality/prior volume of the identifiable SDMC structural/action family remains unresolved."
  },
  "verdict":{
    "conditional_data_fit":"Planck and HSC nuisance evidence nearly cancel for the two fixed cosmologies; the current conditional joint Bayes factor is essentially unity.",
    "null_g2":"The exploratory trajectory-null completion repairs off-trajectory stability without changing accepted-trajectory observables, and therefore adds zero identifiable dimension to that conditional likelihood.",
    "final_model_family":"Still unresolved. A full SDMC-vs-control model-family evidence requires a declared prior or theoretical derivation for the identifiable structural/action coordinates. The null-G2 epsilon is not the source of the outstanding BIC/Bayesian complexity penalty."
  }
}
Path(OUT/"late300_null_g2_model_selection_update.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_NULL_G2_MODEL_SELECTION_UPDATE",json.dumps(out,sort_keys=True))
