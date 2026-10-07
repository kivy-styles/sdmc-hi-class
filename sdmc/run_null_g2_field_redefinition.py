#!/usr/bin/env python3
from pathlib import Path
import json
import sympy as sp

# Structural-clock field redefinition used in the manuscript:
# phi = f(psi), A(psi)=dphi/dpsi, X_phi=A^2 Y_psi.
eps,A,Y,Ys,n=sp.symbols("eps A Y Ys n", nonzero=True)
X=A**2*Y
Xs=A**2*Ys

# General normalized trajectory-null family.
expr=eps*(X-Xs)**n/Xs**(n-1)
expr_simplified=sp.simplify(expr)

# Quartic member that preserves all mixed derivatives through total order 3.
expr4=sp.simplify(expr_simplified.subs(n,4))
expected4=eps*A**2*(Y-Ys)**4/Ys**3
quartic_identity=sp.simplify(expr4-expected4)==0

# The same A^2 prefactor occurs for every n in the normalized family.
general_factor=sp.simplify(expr_simplified/(eps*(Y-Ys)**n/Ys**(n-1)))
general_factor_identity=sp.simplify(general_factor-A**2)==0

# If one insists on writing the transformed completion with the same normalized
# functional form eps_tilde*(Y-Ys)^n/Ys^(n-1), then eps_tilde=eps*A^2.
eps_tilde=sp.simplify(eps*A**2)

# A constant epsilon in phi variables is therefore not a constant scalar under
# an arbitrary nonlinear field redefinition if A=A(psi) varies.
out={
  "status":"trajectory-null G2 field-redefinition covariance audit",
  "input_completion":"Delta G2_phi = epsilon_phi * (X_phi-Xstar_phi)^n / Xstar_phi^(n-1)",
  "field_redefinition":{
    "phi":"f(psi)",
    "A":"dphi/dpsi",
    "X_phi":"A(psi)^2 * Y_psi",
    "Xstar_phi":"A(psi)^2 * Ystar_psi"
  },
  "symbolic_result":{
    "general_transformed_completion":"Delta G2_psi = epsilon_phi*A(psi)^2*(Y-Ystar)^n/Ystar^(n-1)",
    "same_form_coefficient_rule":"epsilon_psi(psi) = epsilon_phi*A(psi)^2",
    "general_factor_verified":bool(general_factor_identity),
    "quartic_identity_verified":bool(quartic_identity),
    "quartic_transformed_completion":"epsilon_phi*A(psi)^2*(Y-Ystar)^4/Ystar^3"
  },
  "minimality":{
    "quartic_power_n4_preserved":True,
    "reason":"A nonsingular field redefinition rescales X and Xstar by the same A^2, so the order of vanishing in (X-Xstar) is unchanged."
  },
  "interpretation":{
    "constant_epsilon_field_coordinate_invariant":False,
    "physical_content":"The off-trajectory deformation is a legitimate covariant scalar contribution once its full coefficient function is transformed, but the numerical statement 'epsilon = constant' is tied to the chosen scalar-field coordinate unless A is constant.",
    "structural_clock_implication":"If SDMC treats psi=ln(S/S0) as the preferred physical structural clock, a completion coefficient should be defined in that preferred variable or through a redefinition-invariant normalization before a numerical epsilon can be called fundamental.",
    "model_selection_implication":"The existing accepted-trajectory non-identifiability result is unchanged. This audit further blocks interpreting the finite-basin optimum epsilon as a field-independent fundamental parameter."
  },
  "next_requirement":"Either specify the structural-clock field as the preferred physical normalization and derive its coefficient, or construct a redefinition-invariant off-trajectory completion from microscopic SDMC quantities."
}
Path("output/null_g2_field_redefinition").mkdir(parents=True,exist_ok=True)
Path("output/null_g2_field_redefinition/null_g2_field_redefinition.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NULL_G2_FIELD_REDEFINITION",json.dumps(out,sort_keys=True))
