#!/usr/bin/env python3
"""Profile official full Planck nuisance sector for named C_l candidate files."""
from pathlib import Path
import argparse,csv,json
import numpy as np
from iminuit import Minuit
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE import TTTEEE
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

ap=argparse.ArgumentParser()
ap.add_argument("--out",required=True)
ap.add_argument("cases",nargs="+")
args=ap.parse_args()
cases=[]
for item in args.cases:
    if "=" not in item: raise ValueError(item)
    label,path=item.split("=",1); cases.append((label,Path(path)))

TCMB=2.7255
high=TTTEEE(packages_path="planck_full")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

fixed={
'cib_index':-1.3,'galf_TE_index':-2.4,'galf_EE_index':-2.4,
'A_sbpx_100_100_TT':1.,'A_sbpx_143_143_TT':1.,'A_sbpx_143_217_TT':1.,'A_sbpx_217_217_TT':1.,
'galf_EE_A_100':.055,'galf_EE_A_100_143':.040,'galf_EE_A_100_217':.094,
'galf_EE_A_143':.086,'galf_EE_A_143_217':.21,'galf_EE_A_217':.70,
'A_cnoise_e2e_100_100_EE':1.,'A_cnoise_e2e_143_143_EE':1.,'A_cnoise_e2e_217_217_EE':1.,
'A_sbpx_100_100_EE':1.,'A_sbpx_100_143_EE':1.,'A_sbpx_100_217_EE':1.,
'A_sbpx_143_143_EE':1.,'A_sbpx_143_217_EE':1.,'A_sbpx_217_217_EE':1.,
'A_pol':1.,'calib_100P':1.021,'calib_143P':.966,'calib_217P':1.040}
names=['A_planck','calib_100T','calib_217T','A_cib_217','xi_sz_cib','A_sz','ksz_norm',
'gal545_A_100','gal545_A_143','gal545_A_143_217','gal545_A_217',
'ps_A_100_100','ps_A_143_143','ps_A_143_217','ps_A_217_217',
'galf_TE_A_100','galf_TE_A_100_143','galf_TE_A_100_217','galf_TE_A_143',
'galf_TE_A_143_217','galf_TE_A_217']
x0=np.array([1.,1.0002,.99805,67.,.05,7.,3.,8.6,10.6,23.5,91.9,
257.,47.,40.,104.,.130,.130,.46,.207,.69,1.938])
current=np.array([1.008245535114707,.9998750001280258,.9982119122897034,
44.22091109885158,.9997030504321224,6.692287855624427,.00016370049661516226,
8.45309812343869,10.60269070310368,20.195272036637586,97.0693107700243,
252.7505102917758,57.67427029377452,63.72262893550619,129.8627383677383,
.11563079474456282,.13456076987301357,.4804478052831501,.22690668463414007,
.6690126964433923,2.1213120205986318])
bounds=[(.97,1.03),(.9946,1.0058),(.99285,1.00325),(0,200),(0,1),(0,10),(0,10),
(0,24.6),(0,26.6),(0,91.5),(0,251.9),(0,400),(0,400),(0,400),(0,400),
(0,.466),(0,.418),(0,1.18),(0,.783),(0,1.41),(0,6.258)]
gauss={'A_planck':(1.,.0025),'calib_100T':(1.0002,.0007),'calib_217T':(.99805,.00065),
'gal545_A_100':(8.6,2.),'gal545_A_143':(10.6,2.),'gal545_A_143_217':(23.5,8.5),
'gal545_A_217':(91.9,20.),'galf_TE_A_100':(.130,.042),'galf_TE_A_100_143':(.130,.036),
'galf_TE_A_100_217':(.46,.09),'galf_TE_A_143':(.207,.072),
'galf_TE_A_143_217':(.69,.09),'galf_TE_A_217':(1.938,.54)}

def load(path):
 a=np.loadtxt(path); ell=a[:,0].astype(int); n=ell.max()+1; conv=(TCMB*1e6)**2
 dl={k:np.zeros(n) for k in ['tt','ee','bb','te','pp','tp','ep']}
 dl['tt'][ell]=a[:,1]*conv; dl['ee'][ell]=a[:,2]*conv; dl['te'][ell]=a[:,3]*conv; dl['bb'][ell]=a[:,4]*conv
 ll=ell.astype(float)*(ell+1.); dl['pp'][ell]=a[:,5]*ll
 dl['tp'][ell]=a[:,6]*np.sqrt(ll)*(TCMB*1e6); dl['ep'][ell]=a[:,7]*np.sqrt(ll)*(TCMB*1e6)
 cl={k:np.zeros(n) for k in ['tt','ee','bb','te']}; fac=np.zeros(n); q=ell>=2; fac[ell[q]]=2*np.pi/ll[q]
 for k in cl: cl[k][ell]=dl[k][ell]*fac[ell]
 return cl,dl

def pri(p):
 d={n:((p[n]-mu)/sig)**2 for n,(mu,sig) in gauss.items()}
 d['ksz_sz_combo']=((p['ksz_norm']+1.6*p['A_sz']-9.5)/3.)**2
 return d

def run(label,path):
 cl,dl=load(path)
 def objective(*x):
  p=dict(fixed); p.update(dict(zip(names,map(float,x))))
  lh=float(high.log_likelihood(cl,**p)); A=p['A_planck']
  lt=float(lowT.log_likelihood(dl['tt'],calib=A)); le=float(lowE.log_likelihood(dl['ee'],calib=A))
  lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
  ll=float(lens.log_likelihood(dl,**lp))
  return float(-2*(lh+lt+le+ll)+sum(pri(p).values())) if np.all(np.isfinite([lh,lt,le,ll])) else 1e100
 best=None
 for seed in [current,x0]:
  m=Minuit(objective,*seed,name=names); m.errordef=1.; m.strategy=1; m.tol=.1
  for n,b in zip(names,bounds): m.limits[n]=b
  for n in names: m.errors[n]=max(1e-5,.02*max(abs(m.values[n]),1.))
  m.migrad(ncall=6000)
  if not m.fmin.is_valid: m.simplex(ncall=2500); m.migrad(ncall=6000)
  row=(float(m.fval),np.array([m.values[n] for n in names]),bool(m.fmin.is_valid),int(m.nfcn))
  if best is None or row[0]<best[0]: best=row
 fval,x,valid,nfcn=best; p=dict(fixed); p.update(dict(zip(names,map(float,x))))
 lh=float(high.log_likelihood(cl,**p)); A=p['A_planck']; lt=float(lowT.log_likelihood(dl['tt'],calib=A)); le=float(lowE.log_likelihood(dl['ee'],calib=A))
 lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}; ll=float(lens.log_likelihood(dl,**lp)); pr=pri(p)
 rec={'model':label,'chi2_profile_total':fval,'chi2_high_full':-2*lh,'chi2_lowT':-2*lt,
 'chi2_lowE':-2*le,'chi2_lensing':-2*ll,'chi2_nuisance_priors':sum(pr.values()),
 'A_planck':A,'A_planck_prior_chi2':pr['A_planck'],'valid':valid,'nfcn':nfcn}
 rec.update({n:p[n] for n in names if n!='A_planck'})
 print("COLD_FULLPLIK_RESULT",json.dumps(rec,sort_keys=True),flush=True); return rec

rows=[run(label,path) for label,path in cases]
base=next(r for r in rows if r['model']=='CURRENT')
for r in rows: r['delta_vs_current']=r['chi2_profile_total']-base['chi2_profile_total']
out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("COLD_FULLPLIK_RANKED",json.dumps(sorted(rows,key=lambda r:r['chi2_profile_total']),sort_keys=True),flush=True)
