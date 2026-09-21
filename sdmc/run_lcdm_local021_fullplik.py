#!/usr/bin/env python3
from pathlib import Path
import csv, json, numpy as np
from iminuit import Minuit
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE import TTTEEE
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output")
TCMB=2.7255
high=TTTEEE(packages_path="planck_full")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

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
x0=np.array([
  1.0,1.0002,0.99805,
  67.,0.05,7.,3.,
  8.6,10.6,23.5,91.9,
  257.,47.,40.,104.,
  .130,.130,.46,.207,.69,1.938
],float)

# Previously validated full-Plik optimum for the frozen exact leader, used only
# as a deterministic nuisance-sector seed.
old_best=np.array([
  1.0009365762341846,0.9997467106847509,0.9981560739731901,
  43.193162287179746,0.9998056674508373,6.895039361411885,0.00015794192579733727,
  8.785878818836895,10.919557463345576,20.3264677879129,96.2535342066465,
  243.54911332511628,53.56831616579629,60.740988423524556,126.16191877781048,
  0.11423183829752921,0.13473694444670203,0.48236745411099513,
  0.22545632687086326,0.6667960509398663,2.0893124002260377
],float)

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
expected=set(high.expected_params)
supplied=set(fixed)|set(names)
if expected!=supplied:
    raise RuntimeError(f"param mismatch missing={sorted(expected-supplied)} extra={sorted(supplied-expected)}")

def load(path):
    a=np.loadtxt(path)
    ell=a[:,0].astype(int)
    n=int(ell.max())+1
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
    fac=np.zeros(n)
    q=ell>=2
    fac[ell[q]]=2*np.pi/ll[q]
    for k in cl:
        cl[k][ell]=dl[k][ell]*fac[ell]
    return cl,dl

def prior_chi2(p):
    c=0.
    for n,(mu,sig) in gauss.items():
        c+=((p[n]-mu)/sig)**2
    c+=((p["ksz_norm"]+1.6*p["A_sz"]-9.5)/3.0)**2
    return float(c)

def make_objective(path):
    cl,dl=load(path)
    def objective(*x):
        p=dict(fixed)
        p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p))
        A=p["A_planck"]
        lt=float(lowT.log_likelihood(dl["tt"],calib=A))
        le=float(lowE.log_likelihood(dl["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        if not all(np.isfinite([lh,lt,le,ll])):
            return 1e100
        return float(-2*(lh+lt+le+ll)+prior_chi2(p))
    def pieces(x):
        p=dict(fixed)
        p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p))
        A=p["A_planck"]
        lt=float(lowT.log_likelihood(dl["tt"],calib=A))
        le=float(lowE.log_likelihood(dl["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        return p,lh,lt,le,ll
    return objective,pieces

def run_one(label,path,seeds):
    fun,pieces=make_objective(path)
    best=None
    for iseed,seed in enumerate(seeds):
        m=Minuit(fun,*seed,name=names)
        m.errordef=1.0
        m.strategy=1
        m.tol=0.1
        for n,b in zip(names,bounds):
            m.limits[n]=b
        for n in names:
            m.errors[n]=max(1e-5,0.02*max(abs(m.values[n]),1.0))
        m.migrad(ncall=6500)
        if not m.fmin.is_valid:
            m.simplex(ncall=3000)
            m.migrad(ncall=6500)
        row=(float(m.fval),np.array([m.values[n] for n in names],float),
             bool(m.fmin.is_valid),int(m.nfcn))
        print("LCDM_CONTROL_FULLPLIK_SEED",label,iseed,*[row[0],row[2],row[3]],flush=True)
        if best is None or row[0]<best[0]:
            best=row
    fval,x,valid,nfcn=best
    p,lh,lt,le,ll=pieces(x)
    row={
      "model":label,
      "chi2_profile_total":fval,
      "chi2_high_full":-2*lh,
      "chi2_lowT":-2*lt,
      "chi2_lowE":-2*le,
      "chi2_lensing":-2*ll,
      "chi2_nuisance_priors":prior_chi2(p),
      "A_planck":p["A_planck"],
      "valid":valid,
      "nfcn":nfcn,
    }
    row.update({n:p[n] for n in names if n!="A_planck"})
    print("LCDM_CONTROL_FULLPLIK_RESULT",json.dumps(row,sort_keys=True),flush=True)
    return row,x

models=[
  ("LCDM_REFERENCE",OUT/"lcdm_reference_00_cl_lensed.dat"),
  ("LCDM_LOCAL021",OUT/"lcdm_local021_00_cl_lensed.dat"),
]
for label,path in models:
    if not path.exists():
        raise FileNotFoundError(path)

results=[]
seed_prev=None
for label,path in models:
    seeds=[x0,old_best]
    if seed_prev is not None:
        seeds.append(seed_prev)
    displaced=x0.copy()
    displaced[3]=80.; displaced[5]=5.; displaced[6]=2.
    seeds.append(displaced)
    row,seed_prev=run_one(label,path,seeds)
    results.append(row)

r={x["model"]:x for x in results}
d=r["LCDM_LOCAL021"]["chi2_profile_total"]-r["LCDM_REFERENCE"]["chi2_profile_total"]
print("LCDM_LOCAL021_FULLPLIK_DELTA",d,flush=True)

fields=list(results[0].keys())
with (OUT/"lcdm_local021_fullplik_profile.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields)
    w.writeheader()
    w.writerows(results)
