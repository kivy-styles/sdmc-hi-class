#!/usr/bin/env python3
import json
from pathlib import Path
import sympy as sp

X,phi,eps=sp.symbols("X phi eps", real=True)
Xs=sp.Function("Xs")(phi)
Y=X-Xs

rows={}
for n in range(1,7):
    F=eps*Y**n/Xs**(n-1)
    checks=[]
    for a in range(0,4):
        for b in range(0,4-a):
            d=sp.diff(F,X,a,phi,b)
            on=sp.simplify(d.subs(X,Xs))
            checks.append({"dX":a,"dphi":b,"total_order":a+b,"value_on_trajectory":str(on),
                           "vanishes":bool(on==0)})
    all3=all(c["vanishes"] for c in checks)
    # first pure-X derivative that survives on trajectory
    first=None
    for m in range(0,n+2):
        v=sp.simplify(sp.diff(F,X,m).subs(X,Xs))
        if v!=0:
            first={"order":m,"value":str(v)}
            break
    rows[str(n)]={"all_mixed_total_order_le_3_vanish":all3,
                  "first_nonzero_pure_X_derivative":first,
                  "checks":checks}

passing=[int(n) for n,r in rows.items() if r["all_mixed_total_order_le_3_vanish"]]
out={
  "status":"trajectory-null G2 polynomial minimality audit",
  "family":"deltaG2 = epsilon*(X-Xs(phi))^n / Xs(phi)^(n-1)",
  "required_null_derivative_order":3,
  "passing_orders":passing,
  "minimal_passing_order":min(passing),
  "quartic_statement":"n=4 is the lowest polynomial power for which deltaG2 and every mixed derivative with total order <=3 vanish identically on X=Xs(phi).",
  "amplitude_statement":"The null conditions fix the minimum power but do not fix epsilon. epsilon remains an off-trajectory completion coefficient unless SDMC supplies an independent microscopic or symmetry condition.",
  "rows":rows
}
Path("output").mkdir(exist_ok=True)
Path("output/null_g2_minimality.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NULL_G2_MINIMALITY",json.dumps(out,sort_keys=True))
