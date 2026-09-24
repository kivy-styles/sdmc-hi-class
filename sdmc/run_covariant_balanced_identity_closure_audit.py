#!/usr/bin/env python3
"""
Covariant Balanced-Identity closure audit for the mature SDMC structural branch.

Part I defines the bare structural state by
    rho_v = rho_P / S^2
    R     = c / sqrt(G rho_v).

The covariant action defines a positive structural source
    rho_X = chi rho_v
and, in the mature No-Slip limit with constant F,
    G_eff/G = 1/F
both in the Friedmann weighting and in the quasi-static No-Slip sector.

Therefore the gravitational strength of the active source relative to the
bare Unified Balanced Identity source is
    Gamma = (G_eff rho_X)/(G rho_v) = chi/F.

The active gravitational radius is
    R_X/R = Gamma^(-1/2) = sqrt(F/chi),

while the mature structural coefficient is
    Xi_active = sqrt(8*pi/3) sqrt(Gamma).

If one imposes the additional covariant matching condition R_X = R, then
Gamma_inf=1, chi_inf=F_inf and N_inf=sqrt(8*pi/3).

This script does not declare that matching principle fundamental.  It compares
it with the alternative bare-source restoration condition chi_inf=1 and
quantifies both against the released future action.
"""
from pathlib import Path
import json,math

SRC=Path("output/future_homogeneous_covariant_stitch_audit.json")
OUT=Path("output/covariant_balanced_identity_closure_audit.json")

d=json.loads(SRC.read_text())
Xi=math.sqrt(8.*math.pi/3.)
F_inf=float(d["accepted_present"]["F_inf"])

out={
  "status":(
    "Comparison of two distinct mature normalization closures. "
    "Covariant Balanced-Identity matching requires chi/F=1; "
    "bare-source restoration requires chi=1. Neither is silently assumed."
  ),
  "identities":{
    "bare_density":"rho_v=rho_P/S^2",
    "bare_structural_radius":"R=c/sqrt(G*rho_v)",
    "active_density":"rho_X=chi*rho_v",
    "mature_no_slip_coupling":"G_eff/G=1/F",
    "gravitational_weight_ratio":"Gamma=(G_eff*rho_X)/(G*rho_v)=chi/F",
    "active_radius_ratio":"R_X/R=sqrt(F/chi)=Gamma^(-1/2)",
    "mature_structural_coefficient":"Xi_active=sqrt(8*pi/3)*sqrt(Gamma)",
    "Gamma_beta":"beta_Gamma=dlnGamma/dlnsigma=2-m_X-alpha_M/p",
    "Gamma_sum_rule":"ln(Gamma_inf/Gamma_0)=integral beta_Gamma dlnsigma",
  },
  "closures":{
    "covariant_balanced_identity":{
      "condition":"G_eff*rho_X=G*rho_v, equivalently R_X=R",
      "Gamma_inf":1.0,
      "chi_inf":F_inf,
      "N_inf":Xi,
      "interpretation":(
        "Preserves the original Unified Balanced Identity radius-density "
        "gravitational strength after the non-minimal Planck mass freezes."
      )
    },
    "bare_source_restoration":{
      "condition":"rho_X=rho_v",
      "Gamma_inf":1./F_inf,
      "chi_inf":1.0,
      "N_inf":Xi/math.sqrt(F_inf),
      "active_radius_ratio_inf":math.sqrt(F_inf),
      "interpretation":(
        "Restores the bare inverse-square source amplitude, while the "
        "frozen F_inf leaves the active gravitational strength lower by 1/F_inf."
      )
    }
  },
  "candidates":{}
}

for name,row in d["noslip_refined_candidates"].items():
    final=row["final"]
    first=min(row["samples"],key=lambda x:abs(float(x["ln_a"])))
    chi=float(final["chi_action"])
    F=float(final["F"])
    gamma=chi/F
    gamma0=float(first["chi_action"])/float(first["F"])
    rr=math.sqrt(F/chi)
    Xi_active=Xi*math.sqrt(gamma)
    N=float(final["N_struct"])
    p=float(final["p"])
    # Exact kinematic coefficient HR/c=N/p for R=l_P S.
    HRc=N/p
    out["candidates"][name]={
      "ln_a":float(final["ln_a"]),
      "chi_action":chi,
      "F":F,
      "Gamma_action":gamma,
      "Gamma_present":gamma0,
      "integrated_beta_Gamma_to_final":math.log(gamma/gamma0),
      "required_integral_for_covariant_match":-math.log(gamma0),
      "required_integral_for_bare_source_match":
        math.log((1./F_inf)/gamma0),
      "G_eff_rhoX_over_G_rhov":gamma,
      "R_active_over_R_structural":rr,
      "Xi_active_from_action":Xi_active,
      "HR_over_c_from_kinematics":HRc,
      "HRc_minus_Xi_active":HRc-Xi_active,
      "N_struct":N,
      "p":p,
      "asymptotic_target_N":float(row["asymptotic_target"]["N_inf"]),
      "distance_to_covariant_match":gamma-1.,
      "distance_to_bare_source_match":chi-1.,
    }

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("COVARIANT_BALANCED_IDENTITY_CLOSURE_AUDIT")
out["flow_sum_rule"]={
  "present_Gamma":next(iter(out["candidates"].values()))["Gamma_present"],
  "covariant_match_required_area":
    next(iter(out["candidates"].values()))["required_integral_for_covariant_match"],
  "bare_source_required_area":
    next(iter(out["candidates"].values()))["required_integral_for_bare_source_match"],
  "area_difference":
    next(iter(out["candidates"].values()))["required_integral_for_covariant_match"]
    -next(iter(out["candidates"].values()))["required_integral_for_bare_source_match"],
  "ln_F_inf":math.log(F_inf),
}
# Rewrite once more with the derived flow block included.
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

for k,v in out["closures"].items():
    print("CLOSURE",k,json.dumps(v,sort_keys=True))
print("FLOW_SUM_RULE",json.dumps(out["flow_sum_rule"],sort_keys=True))
for name,row in out["candidates"].items():
    print("CANDIDATE",name,json.dumps(row,sort_keys=True))
