#!/usr/bin/env python3
from pathlib import Path
import json, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/sharp_tau"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; SIG=.0025
AF=0.036947266645729546; ZC=4.2729806587100025; WIDTH=0.3756495761351945
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def load(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),dict(tt=tt,te=te,ee=ee,pp=pp)
def score(path):
    h,d=load(path)
    def p(A):
        A=float(A)
        ch=float(high.chi_squared(h,A_planck=A)); ct=float(-2*lowT.log_likelihood(d["tt"],calib=A)); ce=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(d,**lp)); cp=((A-1)/SIG)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    o=minimize_scalar(lambda A:p(A)[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(o.x); return (A,)+p(A)

def sdmc(tau,root):
    return textwrap.dedent(f"""\
H0 = 70.8514
omega_b = 0.02239952
omega_cdm = 0.12444227328918850
N_ncdm = 0
N_ur = 3.046
T_cmb = 2.7255
YHe = 0.2453
A_s = 2.16715426219737e-9
n_s = 0.964
tau_reio = {tau:.12g}
Omega_Lambda = 0
Omega_fld = 3.1443554e-8
fluid_equation_of_state = SDMC_TRACKER
cs2_fld = 0.003
use_ppf = no
Omega_smg = -1
gravity_model = sdmc_v3_independent_kinetic
parameters_smg = {AF:.17g},{ZC:.17g},{WIDTH:.17g},0.34231919445927034,1.0,0.045
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
def lcdm(root):
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
lr=str(OUT/"lcdm_"); p=OUT/"lcdm.ini"; p.write_text(lcdm(lr)); cp=subprocess.run(["./class",str(p)],capture_output=True,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout+cp.stderr)
l=score(lr+"00_cl_lensed.dat"); print("SHARPTAU_LCDM",l,flush=True)

taus=[0.0544,0.056,0.058,0.060,0.061,0.062,0.063,0.064,0.065,0.066,0.068,0.070]
rows=[]
for i,tau in enumerate(taus):
    root=str(OUT/f"t{i:02d}_"); ip=OUT/f"t{i:02d}.ini"; ip.write_text(sdmc(tau,root))
    cp=subprocess.run(["./class",str(ip)],capture_output=True,text=True,timeout=240)
    rec={"tau":tau,"status":"FAIL","returncode":cp.returncode}
    if cp.returncode==0:
        s=score(root+"00_cl_lensed.dat")
        rec.update(status="OK",A_planck=s[0],chi2_high=s[1],chi2_lowT=s[2],chi2_lowE=s[3],chi2_lensing=s[4],chi2_cal=s[5],chi2_total=s[6],delta_planck=s[6]-l[-1])
    else: rec["error"]=(cp.stdout+cp.stderr)[-1000:].replace("\n"," | ")
    print("SHARPTAU_POINT",json.dumps(rec,sort_keys=True),flush=True); rows.append(rec)
df=pd.DataFrame(rows); df.to_csv(OUT/"sharp_tau.csv",index=False)
print("SHARPTAU_BEST",df[df.status=="OK"].nsmallest(len(df),"delta_planck").to_dict("records"),flush=True)
