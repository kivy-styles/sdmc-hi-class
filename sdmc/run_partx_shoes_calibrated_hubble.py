#!/usr/bin/env python3
"""
Part X — calibrated Pantheon+SH0ES local-distance-ladder test.

This is deliberately different from the Part IX intercept-profiled Hubble-flow
SN test. Here the Cepheid-calibrator rows are retained. Their CEPH_DIST values
anchor the SN absolute magnitude, while the rows flagged USED_IN_SH0ES_HF
determine the Hubble-flow normalization. The fit therefore contains genuine
absolute-distance / H0 information.

For each supplied exact background we report:
  * the fixed-background calibrated-ladder chi2 (M_B profiled only);
  * the best local H0 obtained by allowing one Hubble-flow distance-scale
    shift while keeping the shape E(z) of that background fixed;
  * the 1-sigma local H0 uncertainty from the full covariance;
  * Delta chi2 between the acoustic-locked fixed H0 and the locally profiled
    H0 (one extra normalization degree of freedom).

No intercept is profiled in a way that can erase the calibrator/H0 information:
M_B is anchored by the Cepheid calibrators, and the Hubble-flow normalization
is a separate fitted column.
"""
from pathlib import Path
import json, math, re, sys
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve

C_LIGHT = 299792.458
DATA = Path("sn_data/PantheonPlus")
OUT = Path("output/partx_shoes_calibrated")
OUT.mkdir(parents=True, exist_ok=True)

def class_table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:
        a=a[1:]
    if len(a)!=n*n:
        raise RuntimeError(f"covariance size mismatch: {len(a)} vs {n*n}")
    return a.reshape(n,n)

pp=pd.read_csv(DATA/"Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
Cfull=read_cov(DATA/"Pantheon+SH0ES_STAT+SYS.cov",len(pp))

is_cal=pp["IS_CALIBRATOR"].to_numpy(int)==1
is_hf=pp["USED_IN_SH0ES_HF"].to_numpy(int)==1
# A calibrator is used as an absolute-distance anchor, never as a Hubble-flow datum.
mask=is_cal | is_hf
sel=np.where(mask)[0]
d=pp.iloc[sel].reset_index(drop=True)
C=Cfull[np.ix_(sel,sel)]
cf=cho_factor(C,lower=True,check_finite=False)

cal=d["IS_CALIBRATOR"].to_numpy(int)==1
hf=(d["USED_IN_SH0ES_HF"].to_numpy(int)==1) & (~cal)
if not cal.any() or not hf.any():
    raise RuntimeError("SH0ES calibrator/Hubble-flow selection failed")
if np.any(d.loc[cal,"CEPH_DIST"].to_numpy(float)<=0):
    raise RuntimeError("invalid CEPH_DIST in calibrator rows")

mobs=d["m_b_corr"].to_numpy(float)
z=d["zHD"].to_numpy(float)
zh=d["zHEL"].to_numpy(float)
ceph=d["CEPH_DIST"].to_numpy(float)

def cinv(v):
    return cho_solve(cf,v,check_finite=False)

def h0_of_bg(bg):
    zarr=bg["z"].to_numpy(float)
    Harr=bg["H [1/Mpc]"].to_numpy(float)
    return float(np.interp(0.0,zarr,Harr)*C_LIGHT)

def mu_hf(bg):
    DM=np.interp(z,bg["z"],bg["comov. dist."])
    DA=DM/(1.+z)
    DL=(1.+zh)*(1.+z)*DA
    if np.any(DL[hf]<=0):
        raise RuntimeError("non-positive luminosity distance in Hubble-flow rows")
    return 5.*np.log10(DL)+25.

def solve_gls(y,A):
    CiA=np.column_stack([cinv(A[:,j]) for j in range(A.shape[1])])
    N=A.T@CiA
    cov=np.linalg.inv(N)
    beta=cov@(A.T@cinv(y))
    r=y-A@beta
    chi=float(r@cinv(r))
    return beta,cov,chi,r

def score(name,path):
    bg=class_table(path)
    H0base=h0_of_bg(bg)
    mucos=mu_hf(bg)
    mu0=np.where(cal,ceph,mucos)

    # y = M_B + delta_mu * I_HF
    y=mobs-mu0
    A=np.column_stack([np.ones(len(y)),hf.astype(float)])
    beta,cov2,chi_prof,_=solve_gls(y,A)
    MB=float(beta[0]); dmu=float(beta[1])
    s_dmu=float(math.sqrt(cov2[1,1]))
    H0_local=float(H0base*10.0**(-dmu/5.0))
    s_H0=float(H0_local*math.log(10.)/5.*s_dmu)

    # Fixed cosmological normalization: only M_B is profiled.
    A1=np.ones((len(y),1))
    b1,cov1,chi_fixed,_=solve_gls(y,A1)
    delta=float(chi_fixed-chi_prof)
    sig=float(math.sqrt(max(delta,0.0)))

    # Calibrator and HF residual RMS after the profiled local fit. These are
    # descriptive only; the covariance chi2 above is the actual statistic.
    rprof=y-A@beta
    rec={
      "model":name,
      "n_total_selected":int(len(y)),
      "n_calibrators":int(cal.sum()),
      "n_hubble_flow":int(hf.sum()),
      "H0_background":H0base,
      "M_B_fixedH0_profile":float(b1[0]),
      "chi2_fixed_H0":chi_fixed,
      "M_B_local_profile":MB,
      "delta_mu_HF":dmu,
      "sigma_delta_mu_HF":s_dmu,
      "H0_local_profile":H0_local,
      "sigma_H0_local":s_H0,
      "chi2_local_profile":chi_prof,
      "delta_chi2_fixed_vs_local_profile":delta,
      "sqrt_delta_chi2":sig,
      "calibrator_rms_mag":float(np.sqrt(np.mean(rprof[cal]**2))),
      "hubble_flow_rms_mag":float(np.sqrt(np.mean(rprof[hf]**2))),
    }
    print("PARTX_SHOES_RESULT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

models=[]
for arg in sys.argv[1:]:
    name,path=arg.split("=",1)
    models.append(score(name,path))

# Cross-background consistency: the local ladder should return nearly the same
# H0 for modestly different viable background shapes if the implementation is sound.
out={"selection":{
        "n_rows_full":int(len(pp)),
        "n_selected":int(mask.sum()),
        "n_calibrators":int(cal.sum()),
        "n_hubble_flow":int(hf.sum())
     },
     "models":models}
if len(models)>=2:
    h=np.array([x["H0_local_profile"] for x in models])
    out["local_H0_spread"]=float(h.max()-h.min())
(OUT/"partx_shoes_calibrated_summary.json").write_text(json.dumps(out,indent=2))
pd.DataFrame(models).to_csv(OUT/"partx_shoes_calibrated_summary.csv",index=False)
print("PARTX_SHOES_SUMMARY",json.dumps(out,sort_keys=True),flush=True)
