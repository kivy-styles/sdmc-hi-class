#!/usr/bin/env python3
"""
Cut-stability audit for the late300 native-Weyl DES-Y3 tomography pilot.

Uses the already generated late300/local021 Weyl+matter spectra from the
successful native-Weyl run and scans progressively more conservative angular
cuts.  The purpose is to identify the regime in which the linear native-Weyl
comparison is stable and to diagnose how strongly the crude HALOFIT-envelope
variant is driven by small scales.

This is still not the official DES-Y3 likelihood: photo-z shifts, shear
calibration, baryonic feedback and eta_IA are not profiled here.
"""
from pathlib import Path
import importlib.util, json, numpy as np, csv

spec=importlib.util.spec_from_file_location("weylbase","sdmc/run_late300_desy3_weyl_tomography.py")
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)

ROOT=Path("output/weyl_cutscan")
ROOT.mkdir(parents=True,exist_ok=True)
DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")

h_l=69.71482083084993/100.0
Om_l=(0.022083219194622913+0.12299536722293603)/h_l**2
h_c=68.56859/100.0
Om_c=(0.02240637+0.11824151)/h_c**2
late=base.load_model("late300","output/weyl_spectra",h_l,Om_l)
lcdm=base.load_model("local021","output/weyl_spectra",h_c,Om_c)
official=base.load_des(DATA,CUTS,extra_conservative=False)

def subset(des,xip_min,xim_min):
    keep=[]
    for i,r in enumerate(des["rows"]):
        threshold=xip_min if r["typ"]=="xip" else xim_min
        if r["ang"]>=threshold: keep.append(i)
    ix=np.array(keep,dtype=int)
    return dict(rows=[des["rows"][i] for i in keep],
                data=des["data"][ix],
                cov=des["cov"][np.ix_(ix,ix)],
                zsrc=des["zsrc"],nz=des["nz"],
                ndata_full=des.get("ndata_full",len(des["data"])))

cuts=[
  ("official",0.0,0.0),
  ("c20_100",20.0,100.0),
  ("c30_120",30.0,120.0),
  ("c40_150",40.0,150.0),
  ("c60_180",60.0,180.0),
  ("c80_200",80.0,200.0),
  ("c100_250",100.0,250.0),
]
out={"status":"native-Weyl DES-Y3 angular-cut stability audit",
     "gr_normalization":base.gr_identity(lcdm),
     "cuts":{}}
rows=[]
for label,xp,xm in cuts:
    des=subset(official,xp,xm)
    if len(des["data"])<20: continue
    rec={"ndata":len(des["data"]),"xip_min_arcmin":xp,"xim_min_arcmin":xm,
         "cov_condition":float(np.linalg.cond(des["cov"]))}
    for variant in ["linear","halofit_envelope"]:
        rr={}
        for m in [late,lcdm]:
            T=base.model_theory(m,des,variant)
            pr=base.chi_profile(des,T[:3])
            pr["valid_grid_fraction"]=T[3]["valid_grid_fraction"]
            rr[m["name"]]=pr
        rr["delta_chi2_late300_minus_local021"]=rr["late300"]["chi2_profile"]-rr["local021"]["chi2_profile"]
        rec[variant]=rr
        rows.append([label,xp,xm,len(des["data"]),variant,
                     rr["late300"]["chi2_profile"],rr["local021"]["chi2_profile"],
                     rr["delta_chi2_late300_minus_local021"],
                     rr["late300"]["AIA_best"],rr["local021"]["AIA_best"]])
    out["cuts"][label]=rec
    print("WEYL_CUTSCAN",label,json.dumps(rec,sort_keys=True),flush=True)

with (ROOT/"late300_weyl_cutscan.csv").open("w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["label","xip_min_arcmin","xim_min_arcmin","ndata","variant",
                "chi2_late300","chi2_local021","delta_chi2","AIA_late300","AIA_local021"])
    w.writerows(rows)
(ROOT/"late300_weyl_cutscan.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
