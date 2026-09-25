#!/usr/bin/env python3
"""
Test the exact one-dimensional canonical infrared manifold.

For the mature r=0 action
    G2 = sigma^-2 (2 F_inf Z - 4 F_inf Z_*),
the positive structural source is
    rho_X = 2 F_inf (Z + 2 Z_*) / sigma^2.
With rho_v = 6 Z_*/sigma^2, the ideal Balanced-Identity relation is
    Gamma = (2 + Z/Z_*)/3
          = [2 + (N_struct/N_inf)^2]/3.

The reconstructed present normalization can carry a tiny constant calibration
factor C_cal.  This audit measures whether the full refined canonical
trajectory approaches
    Gamma = C_cal [2 + (N/N_inf)^2]/3
with one constant C_cal in the mature regime.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_ir_manifold_audit.json")

SPEC={
  "name":"legacy_canonical_ir_manifold",
  "chi_inf":aud.FINF,
  "r":0.0,
  "muQ":50.0,
  "mu_k1":50.0,
  "mu_k2":5.0,
  "mu_V":120.0,
}

model=aud.build_refined_noslip_candidate(SPEC,iterations=2,blend_rate=500.)
res=aud.evolve(model,Nmax=10.,npts=4001)

rows=[]
for s in res["samples"]:
    nrel=float(s["N_struct"]/aud.XI)
    base=(2.+nrel*nrel)/3.
    gamma=float(s["Gamma_action"])
    ccal=gamma/base
    cs2=float(s["cs2_horndeski"])
    rcone=.5*(1./cs2-1.)
    rows.append({
      "ln_a":float(s["ln_a"]),
      "N_over_Ninf":nrel,
      "Gamma_action":gamma,
      "Gamma_ideal_canonical":base,
      "C_cal_local":ccal,
      "r_cone":rcone,
      "cs2":cs2,
    })

late=[r for r in rows if r["ln_a"]>=5.]
Ccal=float(np.median([r["C_cal_local"] for r in late]))
for r in rows:
    r["Gamma_calibrated_canonical"]=Ccal*r["Gamma_ideal_canonical"]
    r["Gamma_minus_calibrated"]=(
      r["Gamma_action"]-r["Gamma_calibrated_canonical"])

# First-order matter prediction from the independently closed m1.
m1=0.2199284631434875
for r in rows:
    u=math.exp(-r["ln_a"])
    r["Gamma_first_order"]=Ccal*(1.-m1*u)

out={
  "status":(
    "canonical infrared manifold audit: tests Gamma=C_cal[2+(N/Ninf)^2]/3 "
    "after the finite Horndeski matching layer"
  ),
  "analytic":{
    "rhoX":"2 F_inf (Z+2 Zstar)/sigma^2",
    "rho_v":"6 Zstar/sigma^2",
    "ideal_Gamma":"[2+Z/Zstar]/3=[2+(N/Ninf)^2]/3",
    "ideal_C_cal":1.0,
    "first_order":"Gamma/C_cal = 1 - m1 a^-1 + O(a^-2 ln a)",
  },
  "fit":{
    "window_ln_a_min":5.0,
    "C_cal":Ccal,
    "C_cal_minus_1":Ccal-1.,
    "late_C_cal_min":float(min(r["C_cal_local"] for r in late)),
    "late_C_cal_max":float(max(r["C_cal_local"] for r in late)),
    "late_C_cal_span":float(max(r["C_cal_local"] for r in late)
                            -min(r["C_cal_local"] for r in late)),
    "max_abs_Gamma_manifold_residual_late":
      float(max(abs(r["Gamma_minus_calibrated"]) for r in late)),
    "m1_independent_mapper_closure":m1,
  },
  "rows":rows,
  "interpretation":{
    "IR_dimension":1,
    "meaning":(
      "Once r_cone has collapsed to zero and the higher Horndeski operators "
      "have decoupled, the remaining normalization relaxation is carried by "
      "the single kinetic/matter coordinate N/N_inf (equivalently Z/Zstar). "
      "Any constant C_cal !=1 is a reconstruction-normalization residual, not "
      "a new dynamical mode."
    )
  }
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("CANONICAL_IR_MANIFOLD_AUDIT")
print("FIT",json.dumps(out["fit"],sort_keys=True))
for r in rows:
    print("ROW",json.dumps(r,sort_keys=True))
