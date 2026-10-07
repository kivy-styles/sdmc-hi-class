#!/usr/bin/env python3
from pathlib import Path
import json
import sympy as sp

# Mature structural form from the manuscript:
#   G2 = sigma^-2 f(Z)
# robust luminality on an open Z interval:
#   cs2(Z)=f_Z/(f_Z+2 Z f_ZZ)=1
# for Z>0 and f_Z != 0 => f_ZZ=0.
#
# Add the structural-clock quartic trajectory-null completion to the mature
# canonical shape and test whether a constant eps_* is compatible.

Z,Zs,A,B,eps=sp.symbols("Z Zs A B eps", positive=True, nonzero=True)
f0=A*Z-B
dq=eps*(Z-Zs)**4/Zs**3
f=sp.expand(f0+dq)
fZ=sp.diff(f,Z)
fZZ=sp.simplify(sp.diff(fZ,Z))
cs2=sp.simplify(fZ/(fZ+2*Z*fZZ))
cs2_minus_1=sp.factor(sp.together(cs2-1))

# Polynomial identity on an open interval: fZZ must vanish identically.
poly=sp.Poly(sp.expand(fZZ*Zs**3),Z)
coeffs=poly.all_coeffs()
identity_requires_eps_zero=all(sp.simplify(c/eps)!=0 for c in coeffs if c!=0) and len(coeffs)>0

# Explicitly solve the identity coefficients = 0.
sol=sp.solve([sp.Eq(c,0) for c in coeffs],eps,dict=True)

out={
  "status":"mature robust-luminality test of structural-clock null-G2 completion",
  "manuscript_mature_form":{
    "G2":"sigma^-2 f(Z)",
    "canonical_f":"A Z - B",
    "robust_luminality":"c_s^2(Z)=1 throughout an open interval",
    "theorem_consequence":"f_ZZ(Z)=0 throughout that interval"
  },
  "tested_completion":{
    "Delta_f":"eps_* (Z-Zstar)^4/Zstar^3",
    "full_f":"A Z - B + eps_* (Z-Zstar)^4/Zstar^3"
  },
  "symbolic":{
    "f_ZZ":str(fZZ),
    "cs2_minus_1":str(cs2_minus_1),
    "open_interval_identity_solution":[{str(k):str(v) for k,v in d.items()} for d in sol],
    "constant_nonzero_eps_compatible":False,
    "identity_requires_eps_zero":bool(identity_requires_eps_zero or sol==[{eps:0}])
  },
  "interpretation":{
    "trajectory_only":"At Z=Zstar the quartic completion and its first three relevant trajectory derivatives vanish, so accepted-trajectory background/linear closure is unchanged.",
    "open_interval":"Away from Zstar, f_ZZ = 12 eps_* (Z-Zstar)^2/Zstar^3. Therefore any constant nonzero eps_* violates exact robust luminality on an open interval.",
    "minimal_derivative":"The same conclusion follows from mature two-derivative minimality: the fully relaxed f(Z) must be linear, so a quartic kinetic-shape term cannot remain in the exact mature action.",
    "relation_to_ir_decay":"The previous sigma^-10 result is still correct along the attracting trajectory. It establishes dynamical irrelevance, not exact equality to the canonical mature action in a finite off-trajectory neighborhood."
  },
  "allowed_role":{
    "transient_handoff":"A nonzero null-G2 stabilizer may be used only as an exploratory finite-handoff/off-trajectory completion if its coefficient flows to zero before the exact mature canonical regime.",
    "mature_constant_completion":"ruled out if robust luminality or minimal-derivative selection is adopted as an exact mature principle."
  },
  "remaining_problem":"Derive a microscopic transverse-flow law for eps_*(psi) that can be nonzero during the finite handoff yet tends to zero in the mature canonical limit, or reject the completion.",
  "verdict":{
    "constant_eps_star_mature_status":"not compatible with exact mature robust luminality/minimal-derivative canonicalization",
    "accepted_late300_status":"unchanged; the null-G2 term remains exploratory and is not part of the accepted reconstructed action",
    "next_theory_target":"transient coupling flow eps_*(psi)->0"
  }
}
Path("output/null_g2_mature_luminality").mkdir(parents=True,exist_ok=True)
Path("output/null_g2_mature_luminality/null_g2_mature_luminality.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n"
)
print("NULL_G2_MATURE_LUMINOSITY",json.dumps(out,sort_keys=True))
