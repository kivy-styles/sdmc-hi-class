#!/usr/bin/env python3
from pathlib import Path
import json

ROOT=Path("deep_provenance_inputs")
OUT=Path("output/late300_deep_parameter_provenance")
OUT.mkdir(parents=True,exist_ok=True)

d0=(ROOT/"d0_audit.py").read_text()
dp=(ROOT/"d0_power_screen.yml").read_text()
lz=(ROOT/"lambda_zt_grid.yml").read_text()
late=(ROOT/"late300_candidate.json").read_text()
candidate=json.loads(late)
zt_summary=json.loads((ROOT/"lambda_selection"/"partial_zeq_edge_best_zt_closure_summary.json").read_text())

def need(cond,msg):
    if not cond:
        raise RuntimeError(msg)

# D0 provenance: repository explicitly describes the canonical value as inherited
# old-r032 and compares it to independently matched/sensitivity values.
need("inherited old-r032 D0=0.34231919445927034" in d0,
     "canonical D0 is not identified as inherited old-r032 in fetched audit")
need("D0_VALUES=[0.24,0.265,D0_MATCH,0.305,0.325,D0_OLD,0.36]" in d0,
     "D0 sensitivity scan not found")

# Kinetic response exponent was subsequently treated as a tunable/sensitivity
# coordinate rather than a uniquely action-derived number.
for p in ['power: "0.30"','power: "0.45"','power: "0.65"','power: "0.90"','power: "1.20"']:
    need(p in dp, f"kinetic exponent scan missing {p}")
for d in ['d0: "0.24"','d0: "0.28"','d0: "0.32"','d0: "0.36"','d0: "0.40"','d0: "0.44"']:
    need(d in dp, f"D0-power matrix missing {d}")

# lambda_e was explicitly scanned jointly with z_t in the earlier r032 family.
need("LAMBDAS=[18.0,18.5,19.0,19.5,20.0]" in lz,
     "lambda_e scan grid not found")
need("ZTS=[18.5,19.0,19.5,20.0,20.5]" in lz,
     "z_t scan grid not found")
need("Omega_X0 is index 0" in lz and "not" in lz,
     "Omega_X0 non-independent scan note missing")

p=candidate["parameters"]
need(abs(p["D0"]-0.34231919445927034)<1e-15,"late300 D0 changed")
need(p["kinetic_power"]==1,"late300 kinetic power changed")
need(abs(p["lambda_e"]-18.40625)<1e-15,"late300 lambda_e changed")
best_zt=zt_summary["best"]
need(best_zt["id"]=="sobol008","historical lambda selection best id changed")
need(abs(best_zt["lambda_e"]-18.40625)<1e-15,"historical winning lambda is not 18.40625")
need(abs(best_zt["chi2_planck"]-1027.4233457902783)<1e-8,"historical lambda-selection score changed")

out={
  "status":"late300 deep inherited-coordinate provenance audit",
  "canonical":{
    "D0":p["D0"],
    "kinetic_power":p["kinetic_power"],
    "lambda_e":p["lambda_e"]
  },
  "repository_evidence":{
    "D0":{
      "classification":"historically inherited and explicitly sensitivity-scanned",
      "evidence":"The physical-D0 audit labels 0.34231919445927034 as inherited old-r032 D0 and scans it against matched and neighboring values.",
      "first_principles_derivation_found":False
    },
    "kinetic_power":{
      "classification":"historically treated as a tunable/sensitivity coordinate",
      "tested_values":[0.30,0.45,0.65,0.90,1.20],
      "canonical_value":1.0,
      "first_principles_derivation_found":False
    },
    "lambda_e":{
      "classification":"explicitly selected in an acoustic-locked Planck closure screen, then frozen in the late300 lineage",
      "example_earlier_scan_values":[18.0,18.5,19.0,19.5,20.0],
      "canonical_value":18.40625,
      "selection_run":35623640746,
      "selection_candidate":"sobol008",
      "selection_chi2_planck":best_zt["chi2_planck"],
      "selection_delta_planck_vs_geom":best_zt["delta_planck_vs_geom"],
      "selection_z_t":best_zt["z_t"],
      "exact_selection_origin_of_18p40625_closed_by_this_audit":True,
      "first_principles_derivation_found":False
    }
  },
  "complexity_consequence":{
    "direct_late300_selection_path_explicit_sdmc_coordinates":7,
    "additional_historically_explored_sdmc_coordinates":["D0","kinetic_power","lambda_e"],
    "historically_explored_union_count":10,
    "warning":"The historically explored union count is NOT an effective AIC/BIC delta_k. Coordinates were explored in staged and partially different model generations, and theory relations may reduce the actual dimensionality.",
    "replay_delta_k":"0 after the accepted action is reconstructed and frozen.",
    "model_family_delta_k":"not uniquely identified from repository history alone.",
    "required_resolution":"Derive these coordinates from SDMC identities/action closure, or include them in a predeclared model-family prior and integrate them in a full evidence calculation."
  }
}
(OUT/"late300_deep_parameter_provenance.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_DEEP_PARAMETER_PROVENANCE",json.dumps(out,sort_keys=True))
