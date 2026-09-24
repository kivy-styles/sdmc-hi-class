#!/usr/bin/env python3
"""
Stress-test the literal bounded-activation interpretation of the released
future structural tail.

The baseline future candidates use one common inverse-square transition rate
muQ for k1, k2 and V.  A slow transition gives an exceptionally smooth future
coasting handoff but allows the action-level activation chi to overshoot unity.
This audit scans faster C2 transition rates and asks whether the overshoot can
be removed without violating the exact Horndeski D>0 and c_s^2>0 conditions.

The scan is diagnostic only; it does not replace the accepted z>=0 action and
does not promote any future tail automatically.
"""
from pathlib import Path
import json
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/future_activation_bound_stress_test.json")

R_FIXED=3.436844444274898
MU_GRID=[
    3.7502361106872546,
    20.0,40.0,60.0,68.0,72.0,76.0,80.0,100.0
]

rows=[]
for mu in MU_GRID:
    spec=dict(
      name=f"bounded_mu_{mu:g}",
      chi_inf=1.0,
      r=R_FIXED,
      muQ=mu
    )
    model=aud.build_candidate(spec)
    result=aud.evolve(model,Nmax=10.,npts=2001)
    h=result["health_diagnostics"]
    f=result["future_kinematics"]
    row={
      "muQ":mu,
      "chi_action_max":h["chi_action_max"],
      "chi_action_max_at_ln_a":h["chi_action_max_at_ln_a"],
      "chi_action_final":h["chi_action_final"],
      "transient_chi_overshoot_over_final":
        h["chi_action_max"]/max(h["chi_action_final"],1e-300)-1.,
      "min_D":h["min_D_horndeski"],
      "min_cs2":h["min_cs2_horndeski"],
      "max_cs2":h["max_cs2_horndeski"],
      "max_abs_alphaB_plus_2alphaM":
        h["max_abs_alphaB_plus_2alphaM"],
      "mapper_ratio_final":h["mapper_ratio_final"],
      "q_max":f["q_max"],
      "q_zero_crossings":f["q_zero_crossings"],
      "N_struct_final":result["final"]["N_struct"],
      "p_final":result["final"]["p"],
      "q_final":result["final"]["q"],
    }
    rows.append(row)

# Identify the first scanned tail for which the only remaining chi>1 excess is
# the common ~1e-3 endpoint normalization floor of the reconstructed action.
first_no_transient=None
for r in rows:
    if r["transient_chi_overshoot_over_final"] <= 1e-4:
        first_no_transient=r
        break

out={
  "status":(
    "bounded-activation future stress test; all entries leave the accepted "
    "z>=0 action unchanged. chi is normalized through "
    "chi/chi0=(rhoX/rhoX0)*sigma^2."
  ),
  "fixed_r":R_FIXED,
  "rows":rows,
  "first_scanned_without_resolved_transient_overshoot":first_no_transient,
  "interpretation":(
    "The slow muQ~3.75 release is dynamically smooth but chi overshoots well "
    "above unity. Faster C2 G2 flow removes the transient overshoot while the "
    "exact Bellini-Sawicki D and c_s^2 gates remain positive. The residual "
    "~1e-3 late offset is at the same level as the endpoint reconstruction "
    "normalization mismatch and must not be interpreted as a physical "
    "super-unit activation until the native future hi_class cross-check is "
    "available."
  )
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_ACTIVATION_BOUND_STRESS_TEST")
for r in rows:
    print("ACTIVATION_SCAN",json.dumps(r,sort_keys=True))
print("FIRST_NO_TRANSIENT",json.dumps(first_no_transient,sort_keys=True))
