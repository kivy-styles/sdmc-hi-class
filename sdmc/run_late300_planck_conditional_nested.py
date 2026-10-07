#!/usr/bin/env python3
"""
Conditional Planck Bayesian evidence for fixed cosmological spectra.

This integrates the full 21-dimensional Planck nuisance sector with dynesty
for a supplied lensed C_ell spectrum. It matches the nuisance likelihood and
priors used by the repository's full-Plik profiles.

Important: cosmological and late300 structural coordinates are FIXED.
Therefore this is a genuine nuisance-marginalized conditional evidence,
not the final model-family evidence for SDMC versus LCDM.

The omitted normalization constants of the shared Gaussian nuisance priors
cancel in the late300/local021 Bayes factor because the priors/bounds are
identical between the two fixed-spectrum models.
"""
from pathlib import Path
import argparse, json, math
import numpy as np
from dynesty import NestedSampler
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE import TTTEEE
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

TCMB=2.7255
fixed={
  "cib_index":-1.3,
  "galf_TE_index":-2.4,
  "galf_EE_index":-2.4,
  "A_sbpx_100_100_TT":1.0,"A_sbpx_143_143_TT":1.0,
  "A_sbpx_143_217_TT":1.0,"A_sbpx_217_217_TT":1.0,
  "galf_EE_A_100":0.055,"galf_EE_A_100_143":0.040,
  "galf_EE_A_100_217":0.094,"galf_EE_A_143":0.086,
  "galf_EE_A_143_217":0.21,"galf_EE_A_217":0.70,
  "A_cnoise_e2e_100_100_EE":1.0,"A_cnoise_e2e_143_143_EE":1.0,
  "A_cnoise_e2e_217_217_EE":1.0,
  "A_sbpx_100_100_EE":1.0,"A_sbpx_100_143_EE":1.0,
  "A_sbpx_100_217_EE":1.0,"A_sbpx_143_143_EE":1.0,
  "A_sbpx_143_217_EE":1.0,"A_sbpx_217_217_EE":1.0,
  "A_pol":1.0,"calib_100P":1.021,"calib_143P":0.966,"calib_217P":1.040,
}
names=[
  "A_planck","calib_100T","calib_217T",
  "A_cib_217","xi_sz_cib","A_sz","ksz_norm",
  "gal545_A_100","gal545_A_143","gal545_A_143_217","gal545_A_217",
  "ps_A_100_100","ps_A_143_143","ps_A_143_217","ps_A_217_217",
  "galf_TE_A_100","galf_TE_A_100_143","galf_TE_A_100_217",
  "galf_TE_A_143","galf_TE_A_143_217","galf_TE_A_217"
]
bounds=[
  (.97,1.03),
  (.9946,1.0058),(.99285,1.00325),
  (0,200),(0,1),(0,10),(0,10),
  (0,24.6),(0,26.6),(0,91.5),(0,251.9),
  (0,400),(0,400),(0,400),(0,400),
  (0,0.466),(0,0.418),(0,1.18),(0,0.783),(0,1.41),(0,6.258)
]
gauss={
  "A_planck":(1.,.0025),
  "calib_100T":(1.0002,.0007),"calib_217T":(.99805,.00065),
  "gal545_A_100":(8.6,2.),"gal545_A_143":(10.6,2.),
  "gal545_A_143_217":(23.5,8.5),"gal545_A_217":(91.9,20.),
  "galf_TE_A_100":(.130,.042),"galf_TE_A_100_143":(.130,.036),
  "galf_TE_A_100_217":(.46,.09),"galf_TE_A_143":(.207,.072),
  "galf_TE_A_143_217":(.69,.09),"galf_TE_A_217":(1.938,.54),
}

def load(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    dl={k:np.zeros(n) for k in ["tt","ee","bb","te","pp","tp","ep"]}
    dl["tt"][ell]=a[:,1]*conv
    dl["ee"][ell]=a[:,2]*conv
    dl["te"][ell]=a[:,3]*conv
    dl["bb"][ell]=a[:,4]*conv
    ll=ell.astype(float)*(ell.astype(float)+1.)
    dl["pp"][ell]=a[:,5]*ll
    dl["tp"][ell]=a[:,6]*np.sqrt(ll)*(TCMB*1e6)
    dl["ep"][ell]=a[:,7]*np.sqrt(ll)*(TCMB*1e6)
    cl={k:np.zeros(n) for k in ["tt","ee","bb","te"]}
    fac=np.zeros(n); q=ell>=2; fac[ell[q]]=2*np.pi/ll[q]
    for k in cl: cl[k][ell]=dl[k][ell]*fac[ell]
    return cl,dl

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--model",required=True)
    ap.add_argument("--spectrum",required=True)
    ap.add_argument("--planck-full",default="planck_full")
    ap.add_argument("--planck-packages",default="planck_packages")
    ap.add_argument("--out",required=True)
    ap.add_argument("--nlive",type=int,default=120)
    ap.add_argument("--dlogz",type=float,default=0.5)
    ap.add_argument("--maxcall",type=int,default=80000)
    args=ap.parse_args()

    high=TTTEEE(packages_path=args.planck_full)
    lowT=TT(packages_path=args.planck_packages)
    lowE=EE(packages_path=args.planck_packages)
    lens=LensingNative(packages_path=args.planck_packages)
    expected=set(high.expected_params); supplied=set(fixed)|set(names)
    if expected!=supplied:
        raise RuntimeError(f"param mismatch missing={sorted(expected-supplied)} extra={sorted(supplied-expected)}")

    cl,dl=load(Path(args.spectrum))
    calls=0

    def prior_chi2(p):
        c=0.
        for n,(mu,sig) in gauss.items():
            c+=((p[n]-mu)/sig)**2
        c+=((p["ksz_norm"]+1.6*p["A_sz"]-9.5)/3.0)**2
        return float(c)

    def loglike(x):
        nonlocal calls
        calls+=1
        p=dict(fixed); p.update(dict(zip(names,map(float,x))))
        try:
            lh=float(high.log_likelihood(cl,**p))
            A=p["A_planck"]
            lt=float(lowT.log_likelihood(dl["tt"],calib=A))
            le=float(lowE.log_likelihood(dl["ee"],calib=A))
            lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
            ll=float(lens.log_likelihood(dl,**lp))
            vals=[lh,lt,le,ll]
            if not all(np.isfinite(vals)): return -1e300
            ans=lh+lt+le+ll-0.5*prior_chi2(p)
            if calls%1000==0:
                print("NESTED_PROGRESS",args.model,calls,ans,flush=True)
            return float(ans)
        except Exception as e:
            if calls<20: print("NESTED_LIKE_EXCEPTION",repr(e),flush=True)
            return -1e300

    lo=np.array([b[0] for b in bounds],float)
    hi=np.array([b[1] for b in bounds],float)
    def ptform(u):
        u=np.asarray(u,float)
        return lo+(hi-lo)*u

    sampler=NestedSampler(loglike,ptform,len(names),nlive=args.nlive,
                          bound="multi",sample="rwalk",
                          rstate=np.random.default_rng(300 if args.model=="late300" else 21))
    sampler.run_nested(dlogz=args.dlogz,maxcall=args.maxcall,print_progress=True)
    r=sampler.results
    logz=float(r.logz[-1]); logzerr=float(r.logzerr[-1])
    ncall=int(np.sum(r.ncall))
    out={
      "status":"conditional fixed-spectrum full-Planck nuisance nested evidence pilot",
      "model":args.model,"spectrum":str(args.spectrum),
      "ndim":len(names),"nlive":args.nlive,"requested_dlogz":args.dlogz,
      "maxcall":args.maxcall,"ncall":ncall,
      "logZ_shared_prior_normalization_omitted":logz,
      "logZerr":logzerr,
      "qualification":(
        "This marginalizes Planck nuisance parameters only. Cosmological and SDMC structural "
        "coordinates are fixed, so it is not the final SDMC-vs-LCDM model-family evidence. "
        "Shared omitted nuisance-prior normalization constants cancel in the Bayes factor."
      )
    }
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    Path(args.out).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("PLANCK_CONDITIONAL_NESTED_EVIDENCE",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
