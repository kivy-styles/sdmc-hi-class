#!/usr/bin/env python3
"""Focused tau-Q-ns refinement on the stable physical No-Slip lensing ridge."""
from pathlib import Path
import json,math,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative
OUT=Path("output/phys_tauqns"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL=.0025
H0=69.7095092787552; OB=.02211778594478965; OC=.1258596659367904
OX=.6953970160452582
AF=.0661676420550793; ZC=3.0664408076554537; WIDTH=.5638187950477004
D0=.34231919445927034; DFLOOR=.045; LAM=17.925; ZT=17.775
NS0=.9637720508351921; TAU0=.06456347454525531; Q0=2.948208928124793
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")
def cls(path):
 a=np.loadtxt(path); ell=a[:,0].astype(int); n=ell.max()+1; conv=(TCMB*1e6)**2
 tt=np.zeros(n);ee=np.zeros(n);te=np.zeros(n);pp=np.zeros(n)
 tt[ell]=a[:,1]*conv;ee[ell]=a[:,2]*conv;te[ell]=a[:,3]*conv
 L=ell.astype(float);pp[ell]=a[:,5]*L*(L+1.)
 return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}
def score(path):
 h,d=cls(path)
 def p(A):
  ch=float(high.chi_squared(h,A_planck=A));ct=float(-2*lowT.log_likelihood(d["tt"],calib=A));ce=float(-2*lowE.log_likelihood(d["ee"],calib=A))
  lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
  cl=float(-2*lens.log_likelihood(d,**lp));cp=((A-1)/CAL)**2
  return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
 o=minimize_scalar(lambda A:p(float(A))[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10}); A=float(o.x); x=p(A)
 return dict(A_planck=A,chi2_high=x[0],chi2_lowT=x[1],chi2_lowE=x[2],chi2_lensing=x[3],chi2_cal=x[4],chi2_total=x[5])
def ini(root,ns,tau,Q):
 lnAs=Q+2*tau; As=math.exp(lnAs)/1e10
 return textwrap.dedent(f"""\
 H0 = {H0}
 omega_b = {OB}
 omega_cdm = {OC}
 N_ncdm = 0
 N_ur = 3.046
 T_cmb = 2.7255
 YHe = 0.2453
 A_s = {As:.17e}
 n_s = {ns:.16g}
 tau_reio = {tau:.16g}
 Omega_Lambda = 0
 Omega_fld = 3.1443554e-8
 fluid_equation_of_state = SDMC_TRACKER
 cs2_fld = 0.003
 use_ppf = no
 Omega_smg = -1
 gravity_model = sdmc_v3_independent_kinetic
 parameters_smg = {AF},{ZC},{WIDTH},{D0},1.0,{DFLOOR}
 expansion_model = sdmc_full
 expansion_smg = {OX},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
 pert_initial_conditions_smg = zero
 method_qs_smg = fully_dynamic
 modes = s
 output = tCl,pCl,lCl
 lensing = yes
 l_max_scalars = 3000
 root = {root}
 format = class
 input_verbose = 0
 background_verbose = 0
 thermodynamics_verbose = 0
 perturbations_verbose = 0
 spectra_verbose = 0
 lensing_verbose = 0
 output_verbose = 0
 """)
anchors=[
 ("current",NS0,TAU0,Q0),
 ("tau0544",NS0,.0544,Q0),
 ("tau056",NS0,.056,Q0),
 ("tau0584",NS0,.0584,Q0),
 ("tau060",NS0,.060,Q0),
 ("tau0584_ns966",.966,.0584,Q0),
 ("tau0584_ns968",.968,.0584,Q0),
]
LOW=np.array([.958,.052,Q0-.008]); HIGH=np.array([.970,.066,Q0+.008])
pts=qmc.scale(qmc.Sobol(3,scramble=True,seed=20260920).random_base2(m=6),LOW,HIGH)
design=anchors+[(f"sobol{i+1:03d}",*map(float,x)) for i,x in enumerate(pts)]
rows=[]
for i,(tag,ns,tau,Q) in enumerate(design,1):
 root=str(OUT/f"{i:03d}_{tag}_"); ip=OUT/f"{i:03d}_{tag}.ini"; ip.write_text(ini(root,ns,tau,Q))
 cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
 rec=dict(id=tag,index=i,n_s=ns,tau_reio=tau,Q=Q,ln10As=Q+2*tau,A_s=math.exp(Q+2*tau)/1e10,returncode=cp.returncode,status="FAIL")
 p=Path(root+"00_cl_lensed.dat")
 if cp.returncode==0 and p.exists(): rec.update(score(p));rec["status"]="OK"
 else: rec["error"]=cp.stdout[-800:].replace("\n"," | ")
 rows.append(rec); print("PHYSTAUQNS_POINT",json.dumps(rec,sort_keys=True),flush=True)
df=pd.DataFrame(rows); base=df[df.id=="current"].iloc[0]
for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]: df["delta_"+k+"_vs_current"]=df[k]-base[k]
df.to_csv(OUT/"phys_tau_q_ns_refine.csv",index=False)
ok=df[df.status=="OK"]
print("PHYSTAUQNS_CURRENT",base.to_dict(),flush=True)
print("PHYSTAUQNS_BEST_TOTAL",ok.nsmallest(15,"chi2_total").to_dict("records"),flush=True)
