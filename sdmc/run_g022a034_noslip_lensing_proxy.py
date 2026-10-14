#!/usr/bin/env python3
"""
Exact-background No-Slip lensing-amplitude correction for g022-A034.

This is NOT a full DES-Y3/HSC/KiDS shear likelihood.  It corrects the earlier
compressed S8-only diagnostic by the exact SDMC Weyl coupling
  Sigma(z) = M_Pl^2 / M_*^2(z) = 1 / Mstar2(z)
for the exact g022-A034 background.  Since cosmic-shear power scales
approximately as Sigma^2 P_delta, the corresponding amplitude-equivalent
compressed variable is S8_eff(z) = Sigma(z) S8.

We report a redshift bracket instead of pretending that one survey has a
single source redshift.
"""
from pathlib import Path
import json,re
import numpy as np,pandas as pd

IN=Path("inputs/a034/desi_cov_00_background.dat")
OUT=Path("output/g022a034_noslip_lensing"); OUT.mkdir(parents=True,exist_ok=True)
S8=0.8383384799534958
LOCAL_S8=0.8074580907496921
WL={
 "DES_Y3_3x2":dict(mu=0.776,lo=0.017,hi=0.017),
 "HSC_Y3_DESIcal":dict(mu=0.805,lo=0.018,hi=0.018),
 "KiDS_Legacy":dict(mu=0.815,lo=0.021,hi=0.016),
}
def wl_chi(x,d):
    sig=d["hi"] if x>=d["mu"] else d["lo"]
    return ((x-d["mu"])/sig)**2
local={k:wl_chi(LOCAL_S8,d) for k,d in WL.items()}

lines=IN.read_text().splitlines()
hdr=[x for x in lines if x.startswith("#") and re.search(r"1:z",x)][-1].lstrip("#").strip()
ms=list(re.finditer(r"(\d+):",hdr)); names=[]
for i,m in enumerate(ms):
    e=ms[i+1].start() if i+1<len(ms) else len(hdr)
    names.append(hdr[m.end():e].strip())
df=pd.DataFrame(np.loadtxt(IN),columns=names).sort_values("z")

zs=[0.0,0.25,0.5,0.75,1.0,1.5,2.0,3.0]
rows=[]
for z in zs:
    m2=float(np.interp(z,df.z,df["M*^2_smg"]))
    geff=float(np.interp(z,df.z,df["G_eff_smg"]))
    slip=float(np.interp(z,df.z,df["slip_eff_smg"]))
    sigma=1.0/m2
    s8eff=S8*sigma
    r={"z":z,"Mstar2":m2,"Sigma":sigma,"Sigma2":sigma*sigma,
       "G_eff":geff,"slip_eta":slip,"S8_raw":S8,"S8_eff":s8eff}
    for k,d in WL.items():
        chi=wl_chi(s8eff,d)
        r[f"{k}_chi2_proxy"]=chi
        r[f"{k}_delta_vs_local021_proxy"]=chi-local[k]
    rows.append(r)
pd.DataFrame(rows).to_csv(OUT/"g022a034_noslip_lensing_redshift_bracket.csv",index=False)

# Main reporting bracket: z=0.5..1.5, roughly covering the bulk lensing kernel
sel=[r for r in rows if 0.5<=r["z"]<=1.5]
summary={
 "scientific_status":"No-Slip-corrected compressed amplitude proxy; not a full shear likelihood",
 "relation":"Sigma=1/Mstar2, eta=1; shear power approximately scales as Sigma^2 P_delta",
 "raw_S8":S8,
 "Sigma_range_z0p5_to_1p5":[min(r["Sigma"] for r in sel),max(r["Sigma"] for r in sel)],
 "Sigma2_range_z0p5_to_1p5":[min(r["Sigma2"] for r in sel),max(r["Sigma2"] for r in sel)],
 "S8_eff_range_z0p5_to_1p5":[min(r["S8_eff"] for r in sel),max(r["S8_eff"] for r in sel)],
 "redshift_rows":rows
}
(OUT/"g022a034_noslip_lensing_summary.json").write_text(json.dumps(summary,indent=2))
print("G022A034_NOSLIP_LENSING",json.dumps(summary,sort_keys=True),flush=True)
