#!/usr/bin/env python3
"""
Candidate-specific JWST/FRESCO likelihood-equivalent stress test for
late300 accepted-action.

This is deliberately not an exact catalogue-level JWST MCMC. It propagates the
already-audited standalone FRESCO efficiency-vs-collapse-amplification curve
onto the exact late300 accepted-action linear spectrum. The exact old-D044 and new-A034
P(k) files are used to correct for the small primordial-amplitude/tilt change.

The nonlinear UV coefficient g3 is kept independent of the linear hi_class
sector. That is essential: the candidate local mapper operator is constructed
to vanish in the linear regime, so it must not be silently identified with A_F.
"""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd

ROOT=Path("inputs")
OUT=Path("output/late300_locked_jwst")
OUT.mkdir(parents=True,exist_ok=True)

old=np.loadtxt(ROOT/"s117/desi_cov_00_z1_pk.dat")
new=np.loadtxt(ROOT/"late300_action/linear_cov_free_00_z1_pk.dat")
kold,pold=old[:,0],old[:,1]
knew,pnew=new[:,0],new[:,1]

# Halo-relevant diagnostic k range. The exact nonlinear HMF would require a
# full sigma(M) integration; here we preserve the previous likelihood-equivalent
# construction and quantify the linear-spectrum correction directly.
ks=np.array([0.1,0.2,0.5,1.0,2.0])
rat=[]
for k in ks:
    po=float(np.interp(k,kold,pold)); pn=float(np.interp(k,knew,pnew))
    rat.append(math.sqrt(pn/po))
r_sigma=float(np.mean(rat))

# Standalone continuation calibration.
GREF=0.01573642223
QREF=0.054
AMGRID=np.array([1.000,1.100,1.163,1.200,1.300,1.360,1.457])
E95=np.array([0.630,0.455,0.365,0.325,0.245,0.205,0.165])
E68=np.array([0.790,0.575,0.465,0.415,0.310,0.270,0.210])

def eps_for(g3):
    q=QREF*g3/GREF
    am=1+3*q
    # Because sigma is slightly lower for A034, the rare-tail quantity depends
    # approximately on A_M * sigma_new/sigma_old.
    aeff=am*r_sigma
    e95=float(np.interp(aeff,AMGRID,E95,left=E95[0],right=E95[-1]))
    e68=float(np.interp(aeff,AMGRID,E68,left=E68[0],right=E68[-1]))
    return dict(g3=float(g3),Qeff=float(q),A_M=float(am),A_eff_with_linear_spectrum=float(aeff),
                epsilon95_equiv=e95,epsilon68_equiv=e68)

# Solve g3 required to restore the old A_eff thresholds after the A034 spectrum shift.
def g_for_target_aeff(target):
    am=target/r_sigma
    q=(am-1)/3
    return GREF*q/QREF

g_af=0.02048146490100771
g95=g_for_target_aeff(1.360)
g68=g_for_target_aeff(1.457)
rows=[eps_for(g) for g in [g_af,g95,g68]]
labels=["g3_equals_AF_diagnostic","epsilon20_95_equiv_target","epsilon20_68_equiv_target"]
for r,l in zip(rows,labels): r["label"]=l
pd.DataFrame(rows).to_csv(OUT/"late300_locked_jwst_candidates.csv",index=False)
summary={
  "sigma_ratio_by_k":{str(k):float(v) for k,v in zip(ks,rat)},
  "mean_sigma_ratio":r_sigma,
  "g3_targets":{"diagnostic_g3_equals_AF":g_af,"epsilon20_95_equiv":g95,"epsilon20_68_equiv":g68},
  "candidates":rows,
  "scientific_status":"conditional likelihood-equivalent stress test; not catalogue-level JWST MCMC and g3 is not yet action-derived"
}
(OUT/"late300_locked_jwst_summary.json").write_text(json.dumps(summary,indent=2))
print("LATE300_LOCKED_JWST",json.dumps(summary,sort_keys=True),flush=True)
