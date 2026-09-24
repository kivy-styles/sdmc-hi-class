#!/usr/bin/env python3
"""
Close the late canonical matter amplitude against the asymptotic mapper.

The resonant canonical audit measures
    Omega_m = m1 a^-1 + 3 m1^2 a^-2 + ...
and
    M = sigma/a -> M_inf.

For the mature r=0 branch,
    H = dot(sigma)/sigma = v_inf/(M_inf a),
so conserved matter gives
    m1 = rho_m0 M_inf^2/(3 F_inf v_inf^2)
       = Omega_m0 M_inf^2/(F_inf p0^2) (N0/N_inf)^2.

This provides an independent algebraic prediction of the fitted resonant
coefficient m1 from the present matter normalization and the mapper endpoint.
"""
from pathlib import Path
import json,math
import run_future_homogeneous_covariant_stitch_audit as aud

SRC=Path("output/canonical_resonant_asymptotic_audit.json")
OUT=Path("output/canonical_mapper_resonance_closure.json")

d=json.loads(SRC.read_text())
fit=d["fit"]
m1_fit=float(fit["m1"])
Minf=float(fit["M_inf_est"])

Omega_m0=float(aud.rho_m0/(3.*aud.H0*aud.H0))
m1_pred=(Omega_m0*Minf*Minf/(aud.FINF*aud.p0_bg*aud.p0_bg)
         *(aud.N0_STRUCT/aud.XI)**2)

lead=1.5*m1_pred

# Absolute mapper normalizations from the accepted Part-IX epoch-mapper audit.
Km0=2.1992619194e20
Kp0=1.1743136282e20
Km_inf=Km0*Minf
Kp_inf=Kp0*Minf

samples=[]
for row in d["samples"]:
    N=float(row["ln_a"])
    u=math.exp(-N)
    first=1.-lead*u
    samples.append({
      "ln_a":N,
      "u":u,
      "N_over_Ninf_actual":float(row["N_over_Ninf"]),
      "N_over_Ninf_first_order_from_present_norm":first,
      "first_order_residual":
        float(row["N_over_Ninf"])-first,
      "p_actual":float(row["p"]),
      "p_first_order_from_present_norm":first,
      "q_actual":float(row["q"]),
      "q_first_order_from_present_norm":-lead*u,
    })

out={
  "status":(
    "independent closure of the canonical resonant matter amplitude from the "
    "present matter normalization and asymptotic structural mapper"
  ),
  "inputs":{
    "Omega_m0_standard":Omega_m0,
    "F_inf":aud.FINF,
    "p0":aud.p0_bg,
    "N0":aud.N0_STRUCT,
    "N_inf":aud.XI,
    "M_inf":Minf,
    "Km0":Km0,
    "Kp0":Kp0,
  },
  "derived":{
    "m1_from_present_normalization":m1_pred,
    "m1_from_resonant_fit":m1_fit,
    "m1_absolute_difference":m1_fit-m1_pred,
    "m1_fractional_difference":m1_fit/m1_pred-1.,
    "leading_lapse_bridge_coefficient_3m1_over_2":lead,
    "Km_inf":Km_inf,
    "Kp_inf":Kp_inf,
    "Km_inf_over_Km0":Minf,
    "Kp_inf_over_Kp0":Minf,
  },
  "samples":samples,
  "interpretation":{
    "m1_not_free":(
      "once the canonical endpoint, conserved matter normalization and "
      "asymptotic mapper M_inf are fixed, the leading a^-1 approach "
      "coefficient is algebraically determined"
    ),
    "first_order_lapse":
      "N/N_inf = 1 - (3/2)m1 a^-1 + O(a^-2 ln a)",
    "first_order_mapper":
      "M/M_inf = 1 + (3/2)m1 a^-1 + O(a^-2 ln a)",
  }
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("CANONICAL_MAPPER_RESONANCE_CLOSURE")
print("DERIVED",json.dumps(out["derived"],sort_keys=True))
for row in samples:
    print("SAMPLE",json.dumps(row,sort_keys=True))
