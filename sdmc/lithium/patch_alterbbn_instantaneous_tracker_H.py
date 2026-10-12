#!/usr/bin/env python3
"""Pin-checked adaptation of historical AlterBBN *background only* to reproduce
the accepted SDMC radiation-tracker relation H²=H_std²/(1-f_r) during BBN.

The earlier rho_dark=const*T^4 approach was not a constant *fraction* of the
total relativistic plasma as e+e− annihilation changes effective g*. Instead
modify every dark-density H computation to rescale the **actual thermal energy
density at each integration step** by (1+dd0), with dd0=fr/(1-fr).

All nuclear rates and abundance equations in the legacy nuclear engine are
left completely unchanged. This is a BBN background proxy, not a new
independent microscopic SDMC action or PRIMAT run.
"""
import sys
from pathlib import Path
path=Path(sys.argv[1])
s=path.read_text()
formulations=[
  ("rho_gamma+rho_epem+rho_neutrinos+rho_neuteq+rho_baryons+rhod",2),
  ("rho_gamma+rho_epem+rho_wimp+rho_neutrinos+rho_neuteq+rho_baryons+rhod",4)
]
for old,count in formulations:
    if s.count(old)!=count:
       raise RuntimeError(f"Historical external source changed. Expected {count} instances of {old}, found {s.count(old)}")
    new="("+old.removesuffix("+rhod")+")*(1.+paramrelic.dd0)"
    s=s.replace(old,new)
assert "rho_gamma+rho_epem+rho_neutrinos+rho_neuteq+rho_baryons+rhod" not in s
assert "rho_gamma+rho_epem+rho_wimp+rho_neutrinos+rho_neuteq+rho_baryons+rhod" not in s
# All six H evaluations now multiply the full instantaneous physical
# radiation/plasma density. No changes to H_SBBN zero-extra-energy branch.
path.write_text(s)
print("CONSTANT_FRACTION_TRACKER_BBN_H2_GATE_SIX_FORMULAE",[(k,c) for k,c in formulations])
