#!/usr/bin/env python3
from pathlib import Path
import json, re

ROOT=Path("provenance_inputs")
OUT=Path("output/late300_parameter_provenance")
OUT.mkdir(parents=True,exist_ok=True)

candidate=json.loads(Path("sdmc/late300_candidate.json").read_text())
growth=(ROOT/"growth_pareto.py").read_text()
prim=(ROOT/"primordial_matrix.yml").read_text()
latebg=(ROOT/"late_background_screen.py").read_text()
struct=(ROOT/"structural_refine.py").read_text()
action=(ROOT/"accepted_action.py").read_text()

def need(cond,msg):
    if not cond:
        raise RuntimeError(msg)

# Mechanical checks of the historical search lineage.
need("Only four perturbative / independent-kinetic structural coordinates vary" in growth,
     "growth-Pareto four-coordinate statement missing")
for token in ["A_F, z_c, width, D_floor","sob=qmc.Sobol(d=4","run(label,AF,ZC,W,DF)"]:
    need(token in growth, f"growth-Pareto token missing: {token}")

need("Stage 1: 256 Sobol backgrounds in (A_late, B_late, H0)" in latebg,
     "late-background search dimension statement missing")
need("sob=qmc.Sobol(d=3" in latebg, "late-background Sobol dimension mismatch")
need("run_background(label,A,B,H0)" in latebg, "late-background runner signature missing")

need("sob=qmc.Sobol(d=5" in struct, "structural-refine Sobol dimension mismatch")
need("run(label,AF,ZC,W,DF,ZT)" in struct, "structural-refine runner signature missing")
for token in ["AF=0.0260","ZC=3.60","W =0.295","DF=0.0400","ZT=15.90"]:
    need(token in struct, f"structural-refine range token missing: {token}")

for label in ["a010","a034","a022","a046"]:
    need(f'label: "{label}"' in prim, f"primordial matrix missing {label}")
for token in ["A_s = ${{ matrix.As }}","n_s = ${{ matrix.ns }}","tau_reio = ${{ matrix.tau }}"]:
    need(token in prim, f"primordial matrix token missing: {token}")

# All g022 structural/background coordinates are literal constants in the a034 promotion.
p=candidate["parameters"]
literal_checks={
    "A_F":p["A_F"],"z_c":p["z_c"],"width":p["width"],"D0":p["D0"],
    "D_floor":p["D_floor"],"lambda_e":p["lambda_e"],"z_t":p["z_t"],
    "A_late":p["A_late"],"B_late":p["B_late"]
}
for name,val in literal_checks.items():
    sval=f"{val}"
    need(sval in prim, f"{name}={sval} not found literally in primordial promotion workflow")

# Accepted-action reconstruction provenance.
need('TARGET_INI = Path("output/linear_cov_target.ini")' in action,
     "accepted action no longer reads the late300 target ini")
need('output/linear_cov_target_00_background.dat' in action,
     "accepted action no longer reads the target background")
for token in ["g=-F1/(H*H)","k2=","k1=","V="]:
    need(token in action, f"accepted-action reconstruction token missing: {token}")

# Omega_x is not an independent coordinate in these search scripts: it is algebraically
# closed from H0 and the physical matter/radiation densities.
need("ox=1.-(OB+OC+OR)/(h*h)" in growth or "ox=1-(OB+OC+OR)/(h*h)" in growth,
     "growth script does not algebraically close Omega_x")
need("ox=1.-(OB+OC+OR)/(h*h)" in latebg or "ox=1-(OB+OC+OR)/(h*h)" in latebg,
     "late-background script does not algebraically close Omega_x")

# a034 changed only the shared primordial trio inside the promoted g022 background.
matrix=[]
pat=re.compile(r'\{ label: "(a\d+)", As: "([^"]+)", ns: "([^"]+)", tau: "([^"]+)" \}')
for m in pat.finditer(prim):
    matrix.append(dict(label=m.group(1),As=float(m.group(2)),ns=float(m.group(3)),tau=float(m.group(4))))
need(len(matrix)==4, f"expected 4 primordial finalists, found {len(matrix)}")

shared_standard=[
    "H0","omega_b","omega_cdm","A_s","n_s","tau_reio","YHe"
]
sdmc_specific_final=[
    "A_F","z_c","width","D0","kinetic_power","D_floor",
    "Omega_x0","lambda_e","z_t","A_late","B_late"
]
historically_scanned_sdmc=[
    "A_F","z_c","width","D_floor","z_t","A_late","B_late"
]
fixed_origin_not_closed=[
    "D0","kinetic_power","lambda_e"
]
derived_not_independent=[
    "Omega_x0"
]

out={
  "status":"late300 parameter-provenance and model-selection audit",
  "canonical_name":candidate["canonical_name"],
  "historical_name":candidate["historical_name"],
  "checks":{
    "g022_growth_search_dimension":4,
    "g022_growth_search_coordinates":["A_F","z_c","width","D_floor"],
    "pre_g022_structural_refine_dimension":5,
    "pre_g022_structural_refine_coordinates":["A_F","z_c","width","D_floor","z_t"],
    "late_background_search_dimension":3,
    "late_background_search_coordinates":["A_late","B_late","H0"],
    "a034_primordial_promotion_coordinates":["A_s","n_s","tau_reio"],
    "a034_finalists":matrix,
    "accepted_action_is_reconstructed_from_target":True,
    "omega_x0_is_algebraically_closed":True
  },
  "parameter_classes":{
    "shared_standard_cosmology":shared_standard,
    "sdmc_specific_final_coordinates":sdmc_specific_final,
    "explicitly_data_scanned_sdmc_coordinates_in_lineage":historically_scanned_sdmc,
    "fixed_in_the_audited_late300_lineage_but_deeper_origin_not_closed_here":fixed_origin_not_closed,
    "algebraically_derived_not_independent":derived_not_independent
  },
  "complexity_interpretation":{
    "free_action_replay_delta_k":"0 at the post-reconstruction replay stage: no structural retuning occurs after the action is frozen.",
    "a034_primordial_trio":"A_s, n_s and tau are ordinary cosmological parameters also present in the control family, so they are not automatically extra SDMC parameters.",
    "explicit_sdmc_search_lower_bound_if_all_remain_independent":
      len(historically_scanned_sdmc),
    "qualification_on_lower_bound":
      "This count is a provenance lower bound on explicitly searched SDMC-specific coordinates, not a final effective AIC/BIC delta_k. The search was staged and adaptive, and some coordinates may be theoretically related or later derived.",
    "why_naive_fixed_action_delta_k_zero_is_not_final":
      "Treating the frozen late300 action as delta_k=0 is valid for replay consistency, but it does not erase the fact that the frozen action was reconstructed from a target selected using cosmological data. A post-selection fixed-model AIC/BIC can therefore under-penalize discovery freedom.",
    "why_naive_delta_k_equal_7_is_not_final":
      "Counting every historically scanned coordinate as one independent final parameter can over-penalize if SDMC identities reduce dimensionality or fix them from deeper theory.",
    "statistically_clean_resolution":
      "Either derive the SDMC-specific coordinates from action/fundamental identities, or declare a reproducible pre-data prior family and compute full model-family evidence. Until then, report AIC/BIC as sensitivity bands rather than a unique final model-selection verdict."
  },
  "canonical_parameters":candidate["parameters"]
}

Path(OUT/"late300_parameter_provenance.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_PARAMETER_PROVENANCE",json.dumps(out,sort_keys=True))
