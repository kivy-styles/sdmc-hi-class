#!/usr/bin/env python3
"""
Part X — official SH0ES-2022 compact full-ladder verification.

Uses the high-level SH0ES y-L-C products released with Riess et al. (2022).
The ladder is a linear Gaussian system:
    y = theta @ L + noise,  noise ~ N(0,C).
We solve the generalized least-squares problem exactly, marginalizing all
other linear ladder parameters, then evaluate the profile penalty for fixed
cosmological H0 values from the exact SDMC backgrounds.

This is complementary to run_partx_shoes_calibrated_hubble.py:
  * that script keeps the exact model D_L(z) shape and uses the Pantheon+
    calibrator + Hubble-flow covariance;
  * this script reproduces the official high-level SH0ES ladder compression
    with all released nuisance directions.
"""
from pathlib import Path
import json, math, sys
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from astropy.io import fits

ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path("official_shoes/SH0ES_Data")
OUT=Path("output/partx_shoes_full_ladder")
OUT.mkdir(parents=True,exist_ok=True)

Yp=ROOT/"ally_shoes_ceph_topantheonwt6.0_112221.fits"
Lp=ROOT/"alll_shoes_ceph_topantheonwt6.0_112221.fits"
Cp=ROOT/"allc_shoes_ceph_topantheonwt6.0_112221.fits"
Rp=ROOT/"lstsq_results.txt"

def load_fits_array(path):
    with fits.open(path,memmap=False) as hdul:
        for h in hdul:
            if h.data is not None:
                return np.asarray(h.data,dtype=float)
    raise RuntimeError(f"no data in {path}")

def load_cov(path):
    try:
        a=load_fits_array(path)
        if a.ndim==2 and a.shape[0]==a.shape[1]:
            return a
    except Exception:
        pass
    vals=np.asarray(np.loadtxt(path),float)
    if vals.ndim==2 and vals.shape[0]==vals.shape[1]:
        return vals
    vals=vals.reshape(-1)
    n=int(round(math.sqrt(vals.size)))
    if n*n==vals.size:
        return vals.reshape(n,n)
    if vals.size>1:
        n0=int(round(vals[0]))
        if n0*n0==vals.size-1:
            return vals[1:].reshape(n0,n0)
    raise RuntimeError(f"cannot parse covariance {path}; shape={vals.shape}, size={vals.size}")

Y=load_fits_array(Yp).reshape(-1)
L=load_fits_array(Lp)
C=load_cov(Cp)

# Official MCMC code uses np.dot(theta,L), so L must be (npar,nobs).
if L.shape[1] != Y.size and L.shape[0] == Y.size:
    L=L.T
if L.shape[1] != Y.size:
    raise RuntimeError(f"L/Y mismatch: L={L.shape}, Y={Y.shape}")
if C.shape != (Y.size,Y.size):
    raise RuntimeError(f"C/Y mismatch: C={C.shape}, Y={Y.shape}")

# Stable generalized least-squares solve.  Avoid forming the normal
# equations directly because the distance-ladder design matrix spans very
# different nuisance scales and can be ill-conditioned.
from scipy.linalg import cholesky, solve_triangular
Lc=cholesky(C,lower=True,check_finite=False)
yw=solve_triangular(Lc,Y,lower=True,check_finite=False)
Aw=solve_triangular(Lc,L.T,lower=True,check_finite=False)  # (nobs,npar)
Q,R=np.linalg.qr(Aw,mode="reduced")
theta=np.linalg.solve(R,Q.T@yw)
Rinv=solve_triangular(R,np.eye(R.shape[0]),lower=False,check_finite=False)
cov_theta=Rinv@Rinv.T
rw=yw-Aw@theta
chi2=float(rw@rw)
res=Y-theta@L

p=float(theta[-1])           # official code: final parameter = 5 log10 H0
sp=float(math.sqrt(cov_theta[-1,-1]))
H0=float(10.0**(p/5.0))
sH0=float(H0*math.log(10.)/5.0*sp)

published_lstsq=np.loadtxt(Rp)
p_pub=float(published_lstsq[-1,0])
sp_pub=float(published_lstsq[-1,1])
H0_pub=float(10.0**(p_pub/5.0))
sH0_pub=float(H0_pub*math.log(10.)/5.0*sp_pub)

fixed={
 "O4D030":71.21323206347779,
 "EDGE004_QBEST":70.78157952616118,
 "EDGE024_QBEST":70.79188564969033,
 "LCDM_LOCAL021":68.56858744695782,
 "LCDM_CONTROL":67.36000006269403,
}
tests=[]
for name,h in fixed.items():
    pp=5.0*math.log10(h)
    dp=pp-p
    dchi=float(dp*dp/(sp*sp))
    tests.append({
      "model":name,
      "H0_fixed":h,
      "delta_5logH0":dp,
      "delta_chi2_profile":dchi,
      "sqrt_delta_chi2":math.sqrt(dchi),
      "signed_sigma":dp/sp,
    })

out={
 "n_data":int(Y.size),
 "n_parameters":int(L.shape[0]),
 "chi2_gls":chi2,
 "dof_nominal":int(Y.size-L.shape[0]),
 "five_log10_H0_gls":p,
 "sigma_five_log10_H0_gls":sp,
 "H0_gls":H0,
 "sigma_H0_gls":sH0,
 "release_lstsq_last_parameter":p_pub,
 "release_lstsq_last_sigma":sp_pub,
 "release_lstsq_H0":H0_pub,
 "release_lstsq_sigma_H0_linearized":sH0_pub,
 "fixed_H0_profile_tests":tests,
}
(OUT/"partx_shoes_full_ladder_summary.json").write_text(json.dumps(out,indent=2))
print("PARTX_SHOES_FULL_LADDER",json.dumps(out,sort_keys=True),flush=True)
for t in tests:
    print("PARTX_SHOES_FULL_FIXED",json.dumps(t,sort_keys=True),flush=True)
