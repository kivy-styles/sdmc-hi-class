#!/usr/bin/env python3
"""
Late300 BIC and asymptotic Bayes-factor sensitivity audit.

Definitions:
  Delta chi2 = chi2_late300 - chi2_local021
  Delta BIC  = Delta chi2 + Delta k * ln(N)
Negative Delta BIC favors late300.

For Delta k > 0, the critical effective data count is
  Ncrit = exp(-Delta chi2 / Delta k)
such that late300 is BIC-favored only for N < Ncrit.

The BIC-derived log Bayes factor is the standard large-sample approximation
  ln BF_late/local ~= -0.5 * Delta BIC.
This is NOT a nested-sampling evidence calculation and is labelled accordingly.
"""
from pathlib import Path
import math, json, csv

OUT=Path("output/late300_model_selection_bic")
OUT.mkdir(parents=True,exist_ok=True)

gates={
 "P+D":-0.12465336634750201,
 "P+D+PantheonPlus":-0.21081492998609974,
 "P+D+Union3":-0.11478786184716228,
 "P+D+DESY5":-0.1907355533444388,
}
shoes_delta=313.57632657154716-323.9418467805046

# We do not pretend a unique iid N exists for mixed correlated cosmology likelihoods.
# Report several sensitivity values plus the exact Ncrit inequality.
Ns=[2,10,30,100,300,1000,2500,5000,10000,30000]
rows=[]
for dataset,dchi in list(gates.items())+[("calibrated_SH0ES_fixed_background_only",shoes_delta)]:
    for dk in range(0,13):
        ncrit=None if dk==0 else math.exp(-dchi/dk)
        for N in Ns:
            dbic=dchi+dk*math.log(N)
            rows.append(dict(
                dataset=dataset,delta_chi2=dchi,delta_k=dk,N_effective=N,
                delta_BIC=dbic,
                BIC_favors_late300=dbic<0,
                lnBF_late_over_local_BIC_approx=-0.5*dbic,
                BF_late_over_local_BIC_approx=math.exp(max(-700,min(700,-0.5*dbic))),
                Ncrit_for_late300_if_dk_gt_0=ncrit,
                status=("diagnostic fixed-background SH0ES only"
                        if dataset.startswith("calibrated_SH0ES")
                        else "exact accepted-action maximum-likelihood gate")
            ))

with (OUT/"late300_BIC_sensitivity.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

thresholds={}
for dataset,dchi in list(gates.items())+[("calibrated_SH0ES_fixed_background_only",shoes_delta)]:
    thresholds[dataset]={}
    for dk in range(1,13):
        thresholds[dataset][str(dk)]=math.exp(-dchi/dk)

primary_max_gain=max(-x for x in gates.values())
summary={
 "definition":{
   "delta_BIC":"delta_chi2 + delta_k * ln(N_effective)",
   "negative_favors":"late300",
   "BIC_Bayes_factor_approximation":"ln BF_late/local ~= -0.5 delta_BIC",
 },
 "exact_gates":gates,
 "critical_N":thresholds,
 "primary_result":(
   "For every exact Planck+DESI(+alternative SN) gate, even delta_k=1 is BIC-disfavored "
   "for any N_effective >= 2. The largest exact likelihood gain is only %.6f, whereas ln(2)=%.6f."
   %(primary_max_gain,math.log(2))
 ),
 "delta_k_zero_result":(
   "If a future derivation makes all late300-specific coordinates fixed predictions so that delta_k=0, "
   "BIC reduces to the raw delta_chi2 and each current exact gate favors late300 slightly."
 ),
 "SH0ES_qualification":(
   "The fixed-background calibrated SH0ES gain is shown only as a sensitivity diagnostic. "
   "It is not a common jointly re-optimized evidence calculation and must not be naively added to Pantheon+."
 ),
 "Bayesian_warning":(
   "The reported BF values are BIC/Laplace asymptotic approximations, not nested-sampling Bayesian evidences."
 )
}
(OUT/"late300_BIC_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print("LATE300_BIC_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
