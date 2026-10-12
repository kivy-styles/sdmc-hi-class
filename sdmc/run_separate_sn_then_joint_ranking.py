#!/usr/bin/env python3
"""Re-rank completed paired CLASS candidate histories in the user's requested
scientific order: three independent SN-only fits FIRST, then their three
independent Planck+DESI-BAO+SN joint scores. No double-counting overlap
between Pantheon+, Union3 and DES-Y5. Same spectra/likelihood archive.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
DATASETS={"pantheonplus":"Pantheon+","union3":"Union3","desy5":"DES-Y5"}
FAMILIES=("lcdm","sdmc")
BASE_KEYS=("family","label","H0","omega_b","omega_cdm","A_s","n_s","tau_reio")

def choose(rows,key,stage):
    df=rows[rows.status=="OK"].copy()
    if not len(df):raise ValueError("No successful candidates")
    out={}
    for family in FAMILIES:
        m=df[df.family==family].sort_values(key,kind="stable")
        if not len(m):raise ValueError("No successful "+family+" candidates")
        baseline=m[m.label=="baseline"]
        best=m.iloc[0]; out[family]={
           "id":str(best["label"]),"H0":float(best["H0"]),
           "score":float(best[key]),
           "baseline_score":(float(baseline.iloc[0][key]) if len(baseline) else None),
           "count":int(len(m))
        }
    out["delta_chi2_sdmc_minus_lcdm"]=float(out["sdmc"]["score"]-out["lcdm"]["score"])
    return out

def run(path,output,stage):
    df=pd.read_csv(path)
    need=set(BASE_KEYS)|{"status","chi2_planck","chi2_bao"}|{"chi2_"+x for x in DATASETS}
    if not need.issubset(df.columns):raise ValueError("Missing fields "+str(sorted(need-set(df.columns))))
    good=df[df.status.eq("OK")].copy()
    output.mkdir(parents=True,exist_ok=True)
    result={"stage":stage,"model_families":list(FAMILIES),
            "scoring":"covariance-weighted uncalibrated SN distance moduli with analytic nuisance/intercept profiling",
            "interpretation":"Separate SN distances do not infer calibrated H0; joint is Planck-lite + DESI BAO + one SN, not DESI full shape",
            "datasets":{}}
    for short,display in DATASETS.items():
        field="chi2_"+short
        metric=field if stage=="sn" else "joint_"+short
        if stage=="joint":
            # No triple SN product. Only ONE SN likelihood per objective.
            good[metric]=good["chi2_planck"]+good["chi2_bao"]+good[field]
        result["datasets"][short]={"label":display,"metric":metric,
              "results":choose(good,metric,stage)}
        rank=good.sort_values(["family",metric],kind="stable")
        cols=list(BASE_KEYS)+[metric,field,"chi2_planck","chi2_bao"]
        cols=list(dict.fromkeys(cols))
        rank[cols].to_csv(output/f"{stage}_{short}_rank.csv",index=False)
    result["success_count_by_family"]={f:int((good.family==f).sum()) for f in FAMILIES}
    (output/f"{stage}_summary.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("STANDALONE_SN_FIRST_RESULT" if stage=="sn" else "ONLY_AFTER_SN_JOINT_RESULT",
          json.dumps(result,sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--stage",choices=["sn","joint"],required=True)
    v=p.parse_args()
    run(v.input,v.output,v.stage)
