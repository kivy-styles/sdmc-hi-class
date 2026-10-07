#!/usr/bin/env python3
"""
Late300 model-selection sensitivity ledger.

This does NOT claim a Bayesian evidence result.  It records the information-
criterion consequence of the exact accepted-action maximum-likelihood gains as
a function of the number of additional effective fitted parameters Delta k.

For AIC:
    Delta AIC = Delta chi2 + 2 Delta k,
where negative favors late300.

The calibrated SH0ES row is shown separately as a fixed-background diagnostic.
It must not be added to the Pantheon+ shape-only gate because the SN content
overlaps, and it is not a substitute for a common joint re-optimization.
"""
from pathlib import Path
import json,pandas as pd
OUT=Path("output/late300_model_selection"); OUT.mkdir(parents=True,exist_ok=True)

gates={
 "P+D":-0.12465336634750201,
 "P+D+PantheonPlus":-0.21081492998609974,
 "P+D+Union3":-0.11478786184716228,
 "P+D+DESY5":-0.1907355533444388,
}
shoes_fixed_late300=313.57632657154716
shoes_fixed_local021=323.9418467805046
shoes_delta=shoes_fixed_late300-shoes_fixed_local021

rows=[]
for gate,dchi in gates.items():
    for dk in range(0,13):
        rows.append(dict(dataset=gate,delta_chi2=dchi,delta_k=dk,
                         delta_AIC=dchi+2*dk,
                         AIC_favors_late300=(dchi+2*dk)<0,
                         status="exact accepted-action likelihood gate"))
for dk in range(0,13):
    rows.append(dict(dataset="calibrated_SH0ES_fixed_background_only",
                     delta_chi2=shoes_delta,delta_k=dk,
                     delta_AIC=shoes_delta+2*dk,
                     AIC_favors_late300=(shoes_delta+2*dk)<0,
                     status="diagnostic only; requires common joint re-optimization before combined model-selection claim"))
df=pd.DataFrame(rows)
df.to_csv(OUT/"late300_AIC_sensitivity.csv",index=False)
summary={
 "exact_gates":gates,
 "calibrated_SH0ES_fixed_background_delta_chi2":shoes_delta,
 "conservative_statement":
   "For the four exact Planck+DESI(+alternative SN) gates, the late300 gain is smaller than 2. Therefore AIC does not favor late300 if the current empirical formulation has even one additional effective fitted parameter relative to local021.",
 "important_qualification":
   "If future theory derives structural coordinates independently rather than fitting them, the effective parameter penalty changes. A proper conclusion requires a common joint fit/evidence calculation with explicit priors.",
 "SH0ES_qualification":
   "The fixed-background calibrated ladder favors late300 by about 10.37 chi2, but the Planck+DESI and SH0ES sectors have not yet been jointly re-optimized under a common model-selection run."
}
(OUT/"late300_model_selection_summary.json").write_text(json.dumps(summary,indent=2))
print("LATE300_MODEL_SELECTION",json.dumps(summary,sort_keys=True),flush=True)
