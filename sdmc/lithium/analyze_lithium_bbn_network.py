#!/usr/bin/env python3
"""Summarize independent (legacy) AlterBBN network and screening calculations.

This is not a complete SDMC field-driven BBN nuclear-network evolution: the
network represents the small early tracker by a radiationlike extra density.
The source library is a pinned, historical AlterBBN fork, not PRIMAT.
"""
import argparse,csv,json,math,re
from pathlib import Path
import numpy as np

OMEGA_LAT=.022083219194622913
OMEGA_LCDM=.022406369378007947
LAMBDA=18.40625
FR=4/LAMBDA**2
OBS={"Yp":{"mean":.2453,"sigma":.0034},
     "D_H":{"mean":2.508e-5,"sigma":.029e-5},
     "Li7_total":{"mean":1.45e-10,"sigma":.25e-10}}

def screen():
    previous={"Yp":.246990,"D_H":2.43955e-5,"Li7_total":5.54707e-10,
              "omega_ref":.02239952,"f_ref":.01}
    f=FR; r=f/previous["f_ref"]; baryon=OMEGA_LAT/previous["omega_ref"]
    return {
      "scope":"LOCAL-SENSITIVITY ONLY; baseline at older omega_b and lambda=20",
      "lambda_e":LAMBDA,"early_tracker_fraction":f,
      "speedup":1/math.sqrt(1-f),"delta_Neff_radiationlike_proxy":
                43/7*f/(1-f),
      "predicted_Yp_approx":previous["Yp"]+
        .0016*(6.115*(OMEGA_LAT-previous["omega_ref"])/.02233)+
        previous["Yp"]*.00336*r,
      "predicted_D_H_approx":previous["D_H"]*baryon**(-1.6)*(1+.00835*r),
      "predicted_Li7_total_approx":previous["Li7_total"]*baryon**2*(1-.00568*r),
      "assumed_baryon_sensitivities":{"D_H":"eta^-1.6","Li7_total":"eta^2"},
      "observational_pdg_2025":OBS
    }

def parse(path):
    rows=[]
    for line in Path(path).read_text().splitlines():
        if "BBN_RESULT " not in line:continue
        raw=line.split("BBN_RESULT ",1)[1].strip()
        d=dict(x.split("=",1) for x in raw.split() if "=" in x)
        label=d.pop("label")
        row={"label":label,**{key:float(val) for key,val in d.items()}}
        rows.append(row)
    if not rows:raise ValueError("Missing direct BBN network results")
    seen={x["label"] for x in rows}
    required={"local021_standard","late300_same_baryon_standard","late300_early_tracker",
              "late300_lambda20","late300_lambda15","late300_lambda10","late300_lambda8"}
    if not required.issubset(seen):raise ValueError(f"missing runs: {required-seen}")
    unique={k:next(x for x in rows if x["label"]==k) for k in seen}
    control=unique["late300_same_baryon_standard"]
    sdmc=unique["late300_early_tracker"]
    ratios={k:sdmc[k]/control[k] for k in ("Yp","D_H","He3_H","Li7_direct","Be7","Li7_total")}
    assert abs(sdmc["fr"]-FR)<1e-9
    assert all(0.9<ratios[k]<1.1 for k in ("Yp","D_H","Li7_total")),ratios
    assert all(abs(x["Li7_total"]-x["Li7_direct"]-x["Be7"])<1e-18 for x in rows)
    vsobs={k:{"network_model_to_observed":sdmc[k]/v["mean"],
               "central_difference":sdmc[k]-v["mean"]}
          for k,v in OBS.items()}
    return {"model":"SDMC Late300 early radiation-tracker equivalent in pinned legacy AlterBBN",
       "scientific_status":"INDEPENDENT NUCLEAR-NETWORK CROSSCHECK with radiationlike density proxy; NOT full PRIMAT or full microscopic H(T)",
       "all_runs":rows,"same_baryon_SDMC_vs_standard_percent":
                     {k:(v-1)*100 for k,v in ratios.items()},
       "observational_comparison":vsobs,"local_sensitivity_screen":screen(),
       "source_pin":"espensem/AlterBBN_master@8219dd4b578244e36d56ebda31541705ae382093",
       "prior_stage":"Part VIII/IX PRIMAT-calibrated local-sensitivity reported -0.568% lithium for lambda_e=20",
       "interpretation":"Neither a small radiation-like tracker nor a small baryon-density change demonstrates a lithium solution. Check D and Yp. No new matter creation and no structural Jacobian-to-nuclear coupling assumed."}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=parse(args.input)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"summary.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    (args.output/"runs.csv").write_text("") # csv module with consistent columns
    with (args.output/"runs.csv").open("w",newline="") as f:
       w=csv.DictWriter(f,fieldnames=list(result["all_runs"][0]))
       w.writeheader();w.writerows(result["all_runs"])
    print("DIRECT_BBN_NETWORK_7BE_LI7_RESULT",json.dumps({
      "same_baryon_percent":result["same_baryon_SDMC_vs_standard_percent"],
      "lithium_network_total":next(x for x in result["all_runs"]
              if x["label"]=="late300_early_tracker")["Li7_total"],
      "observed_lithium":OBS["Li7_total"]["mean"],
      "relative_to_observed":result["observational_comparison"]["Li7_total"]["network_model_to_observed"]
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
