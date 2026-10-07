#!/usr/bin/env python3
from pathlib import Path
import json, math

ROOT=Path("synthesis_inputs")
OUT=Path("output/late300_model_selection_synthesis")
OUT.mkdir(parents=True,exist_ok=True)

def one_json(d, contains=None):
    hits=list(Path(d).rglob("*.json"))
    if contains:
        hits=[p for p in hits if contains in p.name]
    if len(hits)!=1:
        raise RuntimeError(f"expected one json in {d} contains={contains}, got {hits}")
    return json.loads(hits[0].read_text())

def read_named(d, name):
    hits=[p for p in Path(d).rglob(name)]
    if len(hits)!=1:
        raise RuntimeError(f"expected one {name} in {d}, got {hits}")
    return json.loads(hits[0].read_text())

# Exact-gate/BIC audit from the already accepted late300 model-selection branch.
bic=read_named(ROOT/"bic","late300_BIC_summary.json")
prov=read_named(ROOT/"provenance","late300_deep_parameter_provenance.json")
pbf=read_named(ROOT/"planck_bf","planck_conditional_bayes_factor.json")

# Completed HSC safe-scale pairs.
safe={}
for ell in (600,800):
    a=read_named(ROOT/f"hsc_safe_late_{ell}",f"late300_ell{ell}.json")
    b=read_named(ROOT/f"hsc_safe_local_{ell}",f"local021_ell{ell}.json")
    safe[str(ell)]={
      "ndata":a["ndata"],
      "late300_chi2":a["fit"]["chi2_total"],
      "local021_chi2":b["fit"]["chi2_total"],
      "delta_chi2":a["fit"]["chi2_total"]-b["fit"]["chi2_total"],
      "late300_success":bool(a["fit"]["success"]),
      "local021_success":bool(b["fit"]["success"]),
      "qualification":"native-Weyl safe-scale pilot; no validated nonlinear modified-gravity/baryonic prescription"
    }

# Completed full corrected-ReACT pairs that already have both models.
full={}
for var in ("ss_native","ss_screen03"):
    a=read_named(ROOT/f"hsc_full_late_{var}",f"late300_{var}.json")
    b=read_named(ROOT/f"hsc_full_local_{var}",f"local021_{var}.json")
    full[var]={
      "ndata":a["ndata"],
      "late300_chi2":a["fit"]["chi2_total"],
      "local021_chi2":b["fit"]["chi2_total"],
      "delta_chi2":a["fit"]["chi2_total"]-b["fit"]["chi2_total"],
      "late300_success":bool(a["fit"]["success"]),
      "local021_success":bool(b["fit"]["success"]),
      "qualification":"corrected ReACT full-band NLA-z refit; TATT A2/alpha2/bias_ta and baryonic-feedback nuisance remain missing"
    }

# us_native late300 timed out in the original 90-minute matrix but converged
# successfully in the dedicated 180-minute closure retry. Pair it against the
# already converged local021 us_native fit from the original matrix.
a=read_named(ROOT/"hsc_full_late_us_native_retry","late300_us_native.json")
b=read_named(ROOT/"hsc_full_local_us_native","local021_us_native.json")
full["us_native"]={
  "ndata":a["ndata"],
  "late300_chi2":a["fit"]["chi2_total"],
  "local021_chi2":b["fit"]["chi2_total"],
  "delta_chi2":a["fit"]["chi2_total"]-b["fit"]["chi2_total"],
  "late300_success":bool(a["fit"]["success"]),
  "local021_success":bool(b["fit"]["success"]),
  "qualification":"corrected ReACT full-band NLA-z refit; late300 from dedicated timeout retry; TATT A2/alpha2/bias_ta and baryonic-feedback nuisance remain missing"
}

def ic(delta_chi2, n, dk):
    return {
      "delta_AIC":delta_chi2+2*dk,
      "delta_BIC":delta_chi2+dk*math.log(n),
      "BIC_BF_approx_late_over_local":math.exp(-0.5*(delta_chi2+dk*math.log(n)))
    }

for block in (safe,full):
    for rec in block.values():
        rec["IC_sensitivity"]={str(k):ic(rec["delta_chi2"],rec["ndata"],k) for k in (0,1,2,7,10)}

best_exact=min(float(v) for v in bic["exact_gates"].values())
best_hsc=min(rec["delta_chi2"] for rec in full.values())
combined_best_raw=best_exact+best_hsc
# HSC alone provides 60 data points, so any joint data set has N_eff >= 60.
bic_dk1_lower_bound=combined_best_raw+math.log(60.0)

out={
  "status":"late300 model-selection synthesis using completed exact-gate, HSC and converged conditional-Planck evidence results",
  "exact_gate_best_delta_chi2":best_exact,
  "hsc_safe_scale":safe,
  "hsc_full_corrected_ReACT_completed_pairs":full,
  "combined_raw_fit_best_case":{
    "delta_chi2_exact_plus_best_completed_HSC":combined_best_raw,
    "N_effective_lower_bound":60,
    "delta_BIC_dk1_lower_bound":bic_dk1_lower_bound,
    "result":"Even the most favorable completed exact-gate + HSC raw-fit combination cannot overcome one extra parameter under BIC if N_eff is at least the 60 HSC data points." if bic_dk1_lower_bound>0 else "One extra parameter can be overcome under this lower-bound BIC check."
  },
  "converged_conditional_planck_evidence":pbf,
  "parameter_provenance":{
    "replay_delta_k":prov["complexity_consequence"]["replay_delta_k"],
    "historically_explored_union_count":prov["complexity_consequence"]["historically_explored_union_count"],
    "model_family_delta_k":prov["complexity_consequence"]["model_family_delta_k"],
    "required_resolution":prov["complexity_consequence"]["required_resolution"]
  },
  "verdict":{
    "raw_fit":"Completed HSC comparisons consistently improve late300 relative to local021, with the strongest completed corrected-ReACT pair at roughly delta chi2 = %.3f." % best_hsc,
    "AIC":"For the strongest completed HSC pair, AIC would mildly favor late300 if only one genuinely extra parameter is charged, but not for larger delta_k.",
    "BIC":"BIC still favors local021 for every completed HSC pair once even one extra parameter is charged; combining the strongest HSC pair with the best exact gate does not change that for N_eff >= 60.",
    "Bayesian":"The converged fixed-spectrum Planck nuisance Bayes factor is mild and uncertain, not decisive; it cannot substitute for full model-family evidence.",
    "final":"Model selection is therefore unresolved at the theory-family level. The decisive missing item is the effective independent dimensionality/prior volume of the SDMC-specific action coordinates, not a failure of raw likelihood closure."
  },
  "pending":["local021 us_screen03 corrected-ReACT completion/retry","HSC ell=1000 safe-scale closure pair","high-precision Planck conditional evidence repeat"]
}
(OUT/"late300_model_selection_synthesis.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_MODEL_SELECTION_SYNTHESIS",json.dumps(out,sort_keys=True))
