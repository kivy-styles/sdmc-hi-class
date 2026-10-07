#!/usr/bin/env python3
from pathlib import Path
import json
import sympy as sp

# Manuscript structural variables:
#   psi = ln(S/S0)
#   sigma = exp(psi) = S/S0
# and kinetic invariants
#   Y = -1/2 (grad psi)^2
#   Z = -1/2 (grad sigma)^2.
#
# Because dpsi/dsigma = 1/sigma:
#   Y = Z/sigma^2,  Ystar = Zstar/sigma^2.
#
# Start from the normalized trajectory-null family written in the
# logarithmic structural clock psi:
#   Delta G2_psi = eps_* (Y-Ystar)^n / Ystar^(n-1).
# The audit asks what the same covariant deformation looks like in sigma.

eps,sigma,Z,Zs,n=sp.symbols("eps sigma Z Zs n", positive=True, nonzero=True)
Y=Z/sigma**2
Ys=Zs/sigma**2

expr=sp.simplify(eps*(Y-Ys)**n/Ys**(n-1))
expected=sp.simplify(eps/sigma**2*(Z-Zs)**n/Zs**(n-1))
general_identity=sp.simplify(expr-expected)==0

expr4=sp.simplify(expr.subs(n,4))
expected4=sp.simplify(eps/sigma**2*(Z-Zs)**4/Zs**3)
quartic_identity=sp.simplify(expr4-expected4)==0

# For a generic structural coordinate q with psi=f(q),
# E_q(q)=eps_* (dpsi/dq)^2 is the coefficient law that preserves the
# preferred-clock completion.
A=sp.symbols("A", positive=True, nonzero=True)
covariant_coefficient=sp.simplify(eps*A**2)

out={
  "status":"structural-clock normalization audit for trajectory-null G2 completion",
  "manuscript_structural_variables":{
    "log_clock":"psi = ln(S/S0)",
    "direct_field":"sigma = exp(psi) = S/S0",
    "kinetic_relation":"Y = Z/sigma^2",
    "trajectory_kinetic_relation":"Ystar = Zstar/sigma^2"
  },
  "preferred_clock_completion":{
    "general":"Delta G2_psi = eps_* (Y-Ystar)^n / Ystar^(n-1)",
    "quartic":"Delta G2_psi = eps_* (Y-Ystar)^4 / Ystar^3",
    "eps_star_role":"constant normalization defined in the logarithmic structural clock"
  },
  "direct_structural_field_form":{
    "general":"Delta G2_sigma = eps_* sigma^(-2) (Z-Zstar)^n / Zstar^(n-1)",
    "quartic":"Delta G2_sigma = eps_* sigma^(-2) (Z-Zstar)^4 / Zstar^3",
    "general_identity_verified":bool(general_identity),
    "quartic_identity_verified":bool(quartic_identity),
    "coefficient_scaling":"epsilon_sigma(sigma) = eps_* / sigma^2"
  },
  "general_structural_reparametrization":{
    "rule":"For psi=f(q), epsilon_q(q) = eps_* [dpsi/dq]^2",
    "symbolic_coefficient":str(covariant_coefficient),
    "interpretation":"Once psi is declared the physical structural clock, the coefficient function in any reparametrized structural field is fixed by covariance; it is no longer an arbitrary function of field coordinate."
  },
  "connection_to_existing_sdmc_scaling":{
    "match":"The direct-field coefficient scales as sigma^(-2), the same inverse-square structural scaling identified in the manuscript's mature structural dynamics.",
    "what_this_does_establish":"The functional field dependence required by the structural-clock completion is compatible with the manuscript's inverse-square structural scaling.",
    "what_this_does_not_establish":"It does not derive the numerical constant eps_* and does not prove that the null-G2 completion is part of accepted SDMC."
  },
  "model_selection":{
    "accepted_trajectory_identifiability":"unchanged: exactly flat in eps_* for the audited accepted-trajectory background/linear likelihood",
    "regular_parameter_charge":"not justified for that exactly flat conditional likelihood",
    "remaining_theory_question":"derive eps_* from a microscopic/structural principle or introduce an explicit off-trajectory/initial-condition statistical model"
  },
  "verdict":{
    "field_coordinate_ambiguity":"removed conditionally if psi=ln(S/S0) is taken as the preferred physical structural clock",
    "functional_form_ambiguity":"reduced: the sigma-coordinate coefficient must be proportional to sigma^(-2)",
    "normalization_ambiguity":"not removed: eps_* remains undetermined",
    "adoption_status":"exploratory completion remains unadopted pending a principle that fixes eps_*"
  }
}

Path("output/null_g2_structural_clock_normalization").mkdir(parents=True,exist_ok=True)
Path("output/null_g2_structural_clock_normalization/null_g2_structural_clock_normalization.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n"
)
print("NULL_G2_STRUCTURAL_CLOCK_NORMALIZATION",json.dumps(out,sort_keys=True))
