#!/usr/bin/env python3
"""
Candidate-specific JWST/FRESCO likelihood-equivalent stress test for fj028.

This is not a catalogue-level JWST MCMC. It propagates the standalone
FRESCO efficiency-vs-collapse-amplification calibration onto the exact fj028
linear spectrum. The nonlinear mapper coefficient g3 remains independent of
the linear hi_class sector.
"""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd

ROOT=Path("inputs"); OUT=Path("output/fj028_jwst"); OUT.mkdir(parents=True,exist_ok=True)
old=np.loadtxt(ROOT/"d044/desi_cov_00_z1_pk.dat")
new=np.loadtxt(ROOT/"fj028/desi_cov_00_z1_pk.dat")
kold,pold=old[:,0],old[:,1]; knew,pnew=new[:,0],new[:,1]
ks=np.array([0.1,0.2,0.5,1.0,2.0])
rat=[]
for k in ks:
    po=float(np.interp(k,kold,pold)); pn=float(np.interp(k,knew,pnew))
    rat.append(math.sqrt(pn/po))
r_sigma=float(np.mean(rat))

GREF=0.01573642223; QREF=0.054
AMGRID=np.array([1.000,1.100,1.163,1.200,1.300,1.360,1.457])
E95=np.array([0.630,0.455,0.365,0.325,0.245,0.205,0.165])
E68=np.array([0.790,0.575,0.465,0.415,0.310,0.270,0.210])

def eps_for(g3):
    q=QREF*g3/GREF; am=1+3*q; aeff=am*r_sigma
    return dict(g3=float(g3),Qeff=float(q),A_M=float(am),A_eff_with_linear_spectrum=float(aeff),
                epsilon95_equiv=float(np.interp(aeff,AMGRID,E95,left=E95[0],right=E95[-1])),
                epsilon68_equiv=float(np.interp(aeff,AMGRID,E68,left=E68[0],right=E68[-1])))

def g_for_target(target):
    am=target/r_sigma; q=(am-1)/3
    return GREF*q/QREF

# Diagnostic only: g3 = A_F is NOT treated as a fundamental identity.
g_af=0.02111394613161683
g95=g_for_target(1.360); g68=g_for_target(1.457)
rows=[]
for lab,g in [("g3_equals_AF_diagnostic",g_af),
              ("epsilon20_95_equiv_target",g95),
              ("epsilon20_68_equiv_target",g68)]:
    r=eps_for(g); r["label"]=lab; rows.append(r)
pd.DataFrame(rows).to_csv(OUT/"fj028_jwst_candidates.csv",index=False)
summary={"sigma_ratio_by_k":{str(k):float(v) for k,v in zip(ks,rat)},
         "mean_sigma_ratio":r_sigma,
         "g3_targets":{"diagnostic_g3_equals_AF":g_af,
                       "epsilon20_95_equiv":g95,
                       "epsilon20_68_equiv":g68},
         "candidates":rows,
         "scientific_status":"conditional likelihood-equivalent stress test; not catalogue-level JWST MCMC; g3 remains independent of the linear sector"}
(OUT/"fj028_jwst_summary.json").write_text(json.dumps(summary,indent=2))
print("FJ028_JWST",json.dumps(summary,sort_keys=True),flush=True)
