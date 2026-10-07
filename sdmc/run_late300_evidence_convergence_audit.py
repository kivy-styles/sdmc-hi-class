#!/usr/bin/env python3
from pathlib import Path
import json, math

ROOT=Path("evidence_inputs")
OUT=Path("output/evidence_convergence_audit")
OUT.mkdir(parents=True,exist_ok=True)

def load_json(root, model):
    hits=list(root.rglob(f"{model}.json"))
    if len(hits)!=1:
        raise RuntimeError(f"expected one {model}.json under {root}, got {hits}")
    return json.loads(hits[0].read_text())

pilot_l=load_json(ROOT/"pilot_late300","late300")
pilot_c=load_json(ROOT/"pilot_local021","local021")
conv_l=load_json(ROOT/"conv_late300","late300")
conv_c=load_json(ROOT/"conv_local021","local021")

def pair(a,b):
    d=a["logZ_shared_prior_normalization_omitted"]-b["logZ_shared_prior_normalization_omitted"]
    e=math.hypot(a["logZerr"],b["logZerr"])
    return {
      "logZ_late300":a["logZ_shared_prior_normalization_omitted"],
      "logZerr_late300":a["logZerr"],
      "logZ_local021":b["logZ_shared_prior_normalization_omitted"],
      "logZerr_local021":b["logZerr"],
      "delta_logZ_late300_minus_local021":d,
      "delta_logZ_error_quadrature":e,
      "BF_late300_over_local021":math.exp(d),
      "late300_status":a.get("evidence_status","unrecorded"),
      "local021_status":b.get("evidence_status","unrecorded"),
      "late300_ncall":a.get("ncall"),
      "local021_ncall":b.get("ncall"),
      "late300_maxcall":a.get("maxcall"),
      "local021_maxcall":b.get("maxcall"),
      "late300_terminated_by_maxcall":a.get("terminated_by_maxcall"),
      "local021_terminated_by_maxcall":b.get("terminated_by_maxcall"),
    }

pilot=pair(pilot_l,pilot_c)
conv=pair(conv_l,conv_c)

# Mechanical validation of the scientific distinction.
if abs(conv["delta_logZ_late300_minus_local021"] - (-1.1869226315016022)) > 1e-9:
    raise RuntimeError("converged delta logZ changed")
if conv_l.get("evidence_status")!="requested_dlogz_reached" or conv_c.get("evidence_status")!="requested_dlogz_reached":
    raise RuntimeError("converged jobs did not reach requested dlogz")
if conv_l.get("terminated_by_maxcall") or conv_c.get("terminated_by_maxcall"):
    raise RuntimeError("converged jobs unexpectedly hit maxcall")

out={
 "status":"conditional Planck nested-evidence convergence audit",
 "pilot_run":37621316765,
 "converged_run":37629170856,
 "pilot":pilot,
 "converged":conv,
 "change":{
   "delta_logZ_shift":conv["delta_logZ_late300_minus_local021"]-pilot["delta_logZ_late300_minus_local021"],
   "BF_ratio_converged_over_pilot":conv["BF_late300_over_local021"]/pilot["BF_late300_over_local021"],
   "combined_sigma_from_zero":abs(conv["delta_logZ_late300_minus_local021"])/conv["delta_logZ_error_quadrature"],
 },
 "verdict":{
   "pilot_status":"retired_for_model_interpretation_due_to_nonconvergence",
   "converged_status":"current_conditional_fixed_spectrum_planck_nuisance_evidence",
   "interpretation":"The converged conditional Planck nuisance evidence mildly favors local021, but the shift is comparable to the nested-evidence uncertainty and is not a decisive fixed-spectrum Bayes-factor result.",
   "scope_warning":"This is not final SDMC-vs-control model-family evidence: cosmological and SDMC structural/action coordinates are fixed and their prior volumes are not integrated."
 }
}
(OUT/"late300_planck_evidence_convergence_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_EVIDENCE_CONVERGENCE_AUDIT",json.dumps(out,sort_keys=True))
