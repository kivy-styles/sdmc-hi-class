#!/usr/bin/env python3
"""
Derive the C3 No-Slip tail-root equation for the mature canonical endpoint.

For x=ln(sigma), let the present jets be
  F_n = d^n F/dx^n|_0,  g_n = d^n g/dx^n|_0.
The exponential-polynomial C3 tails have cubic leading coefficients
  A_F3(mu) = F3+3 mu F2+3 mu^2 F1+mu^3(F0-Finf),
  A_g3(mu+1) = g3+3(mu+1)g2+3(mu+1)^2 g1+(mu+1)^3 g0.

No-Slip requires g=-F_,sigma/(2 Z).  Because the sigma derivative adds one
extra inverse power of sigma, mu_g=mu_F+1, and the asymptotic leading terms
give
  mu A_F3(mu) = 2 Zstar A_g3(mu+1).

This is a quartic equation for mu_F.  The audit lists all real roots and the
positive branches.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_planck_braiding_tail_root_audit.json")

F0,F1,F2,F3=[float(x) for x in aud.bars["F"]]
g0,g1,g2,g3=[float(x) for x in aud.g_xder]

Zstar=float(aud.Z0*(aud.XI/aud.N0_STRUCT)**2)

# np.poly1d uses descending powers.
AF3=np.poly1d([F0-aud.FINF,3.*F1,3.*F2,F3])
mu=np.poly1d([1.,0.])
t=np.poly1d([1.,1.])
Ag3=g0*t**3+3.*g1*t**2+3.*g2*t+g3
P=mu*AF3-2.*Zstar*Ag3

roots=np.roots(P)
rows=[]
for z in roots:
    zr=float(np.real(z)); zi=float(np.imag(z))
    row={
      "real":zr,
      "imag":zi,
      "is_real":abs(zi)<1e-9,
      "is_positive_real":abs(zi)<1e-9 and zr>0.,
    }
    if row["is_positive_real"]:
        muv=zr
        af=float(AF3(muv)); ag=float(Ag3(muv))
        row.update({
          "mu_F":muv,
          "mu_g":muv+1.,
          "A_F3":af,
          "A_g3":ag,
          "root_residual":float(P(muv)),
          "Z_from_tail_ratio":float(muv*af/(2.*ag)),
          "relative_F_decay_exponent":muv,
          "g_decay_exponent":muv+1.,
        })
    rows.append(row)

positive=sorted([r for r in rows if r["is_positive_real"]],key=lambda r:r["real"])
selected=min(positive,key=lambda r:abs(r["real"]-5.57)) if positive else None

out={
  "status":"quartic C3 No-Slip tail-root audit for the canonical endpoint",
  "present_jets":{
    "F":[F0,F1,F2,F3],
    "g":[g0,g1,g2,g3],
    "F_inf":float(aud.FINF),
    "Zstar":Zstar,
  },
  "equation":{
    "A_F3":"F3+3 mu F2+3 mu^2 F1+mu^3(F0-Finf)",
    "A_g3":"g3+3(mu+1)g2+3(mu+1)^2 g1+(mu+1)^3 g0",
    "quartic":"mu A_F3 - 2 Zstar A_g3 = 0",
    "polynomial_coefficients_descending":[float(x) for x in P.c],
    "mu_g_relation":"mu_g=mu_F+1 from F_,sigma=sigma^-1 F_,x",
  },
  "roots":rows,
  "positive_real_roots":[float(r["real"]) for r in positive],
  "selected_fast_branch":selected,
  "interpretation":{
    "mu_F_not_free":(
      "within the chosen C3 exponential-polynomial tail ansatz, mu_F is fixed "
      "by the accepted present F/g jets, the canonical endpoint Zstar, and "
      "asymptotic No-Slip"
    ),
    "branch_nonuniqueness":(
      "multiple positive roots, if present, represent distinct interpolation "
      "branches reaching the same mature endpoint; health/perturbation tests "
      "are needed to select among them"
    )
  }
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("CANONICAL_PLANCK_BRAIDING_TAIL_ROOT_AUDIT")
print("POLYNOMIAL",json.dumps(out["equation"],sort_keys=True))
print("POSITIVE_ROOTS",json.dumps(out["positive_real_roots"]))
print("SELECTED",json.dumps(out["selected_fast_branch"],sort_keys=True))
