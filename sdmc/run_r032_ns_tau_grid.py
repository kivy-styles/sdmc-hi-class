#!/usr/bin/env python3
from pathlib import Path
import json, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/ns_tau"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def p(A):
        A=float(A)
        hi=float(high.chi_squared(harr,A_planck=A))
        lt=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        le=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        ln=float(-2*lens.log_likelihood(dls,**lp))
        cal=((A-1.)/CAL_SIGMA)**2
        return hi,lt,le,ln,cal,hi+lt+le+ln+cal
    o=minimize_scalar(lambda A:p(A)[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(o.x); return (A,)+p(A)

def sdmc_ini(ns,tau,root):
    return textwrap.dedent(f"""\
    H0 = 70.8514
    omega_b = 0.02239952
    omega_cdm = 0.12444227328918850
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = 2.16715426219737e-9
    n_s = {ns}
    tau_reio = {tau}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = 0.0715,5.0,1.426,0.34231919445927034,1.0,0.045
    expansion_model = sdmc_full
    expansion_smg = 0.7073985893,17.925,17.775,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
    write background = yes
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

def lcdm_ini(root):
    return textwrap.dedent(f"""\
    H0 = 67.36
    omega_b = 0.02237
    omega_cdm = 0.1200
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = 2.10e-9
    n_s = 0.9649
    tau_reio = 0.0544
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
    write background = yes
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

lr=str(OUT/"lcdm_"); li=OUT/"lcdm.ini"; li.write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(li)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
l=score(lr+"00_cl_lensed.dat")
print("NSTAUGRID_LCDM",l,flush=True)

NS=[0.960,0.964,0.968,0.972,0.976]
TAU=[0.050,0.0544,0.058,0.062,0.066,0.070]
rows=[]
for ns in NS:
  for tau in TAU:
    tag=f"n{ns:.3f}_t{tau:.4f}".replace(".","p")
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(sdmc_ini(ns,tau,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"n_s":ns,"tau":tau,"returncode":cp.returncode,"status":"FAIL"}
    if cp.returncode==0 and Path(root+"00_cl_lensed.dat").exists():
        p=score(root+"00_cl_lensed.dat")
        rec.update(status="OK",A_planck=p[0],chi2_high=p[1],chi2_lowT=p[2],chi2_lowE=p[3],
                   chi2_lensing=p[4],chi2_cal=p[5],chi2_total=p[6],delta_planck=p[6]-l[-1])
    else: rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    print("NSTAUGRID_POINT",json.dumps(rec,sort_keys=True),flush=True); rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"ns_tau_grid.csv",index=False)
ok=df[df.status=="OK"]
print("NSTAUGRID_BEST",ok.nsmallest(15,"delta_planck").to_dict("records"),flush=True)
print("NSTAUGRID_BEST_CAL",ok.nsmallest(10,"chi2_cal").to_dict("records"),flush=True)
print("NSTAUGRID_BEST_HIGH",ok.nsmallest(10,"chi2_high").to_dict("records"),flush=True)
