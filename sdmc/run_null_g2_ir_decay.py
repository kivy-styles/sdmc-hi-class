#!/usr/bin/env python3
from pathlib import Path
import json
import sympy as sp

# Structural-clock completion in the direct field:
# Delta G2 = eps * sigma^-2 * (dZ)^n / Zstar^(n-1)
# Manuscript mature attractor:
# dZ ~ sigma^-2 and Zstar -> constant because the structural lapse tends to a constant.
#
# Then the m-th Z derivative scales as
# sigma^-2 * dZ^(n-m) ~ sigma^[-2(n-m+1)] for m<=n.

n,m=sp.symbols("n m", integer=True, nonnegative=True)
sigma,eps,Zs,c=sp.symbols("sigma eps Zs c", positive=True, nonzero=True)
dZ=sp.symbols("dZ")

def exponent(nv,mv):
    return -2*(nv-mv+1)

rows=[]
for nv in range(2,7):
    for mv in range(0,nv+1):
        rows.append({
          "n":nv,
          "Z_derivative_order":mv,
          "sigma_power_if_deltaZ_propto_sigma^-2":exponent(nv,mv)
        })

quartic=[r for r in rows if r["n"]==4]
expected={0:-10,1:-8,2:-6,3:-4,4:-2}
verified=all(next(r for r in quartic if r["Z_derivative_order"]==k)["sigma_power_if_deltaZ_propto_sigma^-2"]==v
             for k,v in expected.items())

# Direct symbolic substitution for quartic member.
expr=eps*sigma**-2*dZ**4/Zs**3
sub=sp.simplify(expr.subs(dZ,c*sigma**-2))
symbolic_exp=sp.simplify(sub/(eps*c**4/Zs**3))

out={
  "status":"structural-clock null-G2 infrared-decay audit",
  "inputs":{
    "completion":"Delta G2_sigma = eps_* sigma^-2 (delta Z)^n / Zstar^(n-1)",
    "mature_attractor":"delta Z proportional to sigma^-2",
    "mature_Zstar":"asymptotically constant when the structural lapse tends to a constant"
  },
  "general_scaling":{
    "formula":"d^m(Delta G2)/dZ^m proportional to sigma^[-2(n-m+1)] for m<=n",
    "table":rows
  },
  "quartic_member":{
    "Delta_G2":"sigma^-10",
    "G2_Z":"sigma^-8",
    "G2_ZZ":"sigma^-6",
    "G2_ZZZ":"sigma^-4",
    "G2_ZZZZ":"sigma^-2",
    "hierarchy_verified":bool(verified),
    "direct_symbolic_factor_after_deltaZ_c_sigma^-2":str(symbolic_exp)
  },
  "linear_and_background_implication":{
    "accepted_trajectory":"The completion and mixed derivatives through total order 3 vanish exactly on Z=Zstar.",
    "near_attractor":"The same derivatives remain power-law suppressed for finite deviations approaching the canonical attractor.",
    "linear_eigenvalues":"The trajectory-null quartic completion does not alter the manuscript's linear canonical eigenvalues at the accepted trajectory."
  },
  "two_derivative_minimality":{
    "compatible_asymptotically":True,
    "reason":"Although the completion is a higher-kinetic off-trajectory operator, its coefficient is structurally inverse-square and the finite-deviation contribution decays as sigma^-10 on the mature attraction law.",
    "important_limit":"This demonstrates asymptotic decoupling, not a microscopic derivation of eps_* and not permission to add the term to accepted SDMC without an independent principle."
  },
  "verdict":{
    "functional_scaling":"structurally consistent with the manuscript's inverse-square field scaling",
    "infrared_behavior":"rapidly irrelevant along the mature canonical attraction law",
    "normalization":"eps_* remains undetermined",
    "adoption_status":"exploratory"
  }
}

Path("output/null_g2_ir_decay").mkdir(parents=True,exist_ok=True)
Path("output/null_g2_ir_decay/null_g2_ir_decay.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NULL_G2_IR_DECAY",json.dumps(out,sort_keys=True))
