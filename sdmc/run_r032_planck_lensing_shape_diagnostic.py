#!/usr/bin/env python3
"""
Diagnostic only: scan CLASS A_L at the frozen exact-leader SDMC point.

A_L rescales the CMB lensing-potential transfer used to lens TT/TE/EE. This
test asks whether the remaining Planck mismatch has the shape of excess/deficit
lensing smoothing. A_L is NOT promoted to an SDMC model parameter.
"""
from pathlib import Path
import json, math, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/planck_lensing_shape")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; SIG=.0025
AL=[0.65,0.75,0.85,0.90,0.95,1.00,1.05,1.10,1.15,1.25,1.35]
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def load(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load(path)
    def pcs(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/SIG)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    o=minimize_scalar(lambda x:pcs(x)[-1],bounds=(.97,1.03),method="bounded",
                      options={"xatol":1e-10})
    A=float(o.x); return (A,)+pcs(A)

def ini(al,root):
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
    tau_reio = 0.0544

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = 0.0715,5.00,1.426,0.34231919445927034,1.0,0.045
    expansion_model = sdmc_full
    expansion_smg = 0.7073985893,17.925,17.775,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3

    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    A_L = {al:.12g}
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

rows=[]
for i,al in enumerate(AL,1):
    root=str(OUT/(f"al{i:02d}_"))
    p=OUT/(f"al{i:02d}.ini"); p.write_text(ini(al,root))
    cp=subprocess.run(["./class",str(p)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec={"A_L":al,"returncode":cp.returncode,"status":"FAIL"}
    cl=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and cl.exists():
        v=score(cl)
        rec.update({"status":"OK","A_planck":v[0],"chi2_high":v[1],
                    "chi2_lowT":v[2],"chi2_lowE":v[3],"chi2_lensing":v[4],
                    "chi2_cal":v[5],"chi2_total":v[6]})
    else:
        rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    print("PLANCK_LENSING_SHAPE_POINT",json.dumps(rec,sort_keys=True),flush=True)
    rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"planck_lensing_shape_scan.csv",index=False)
ok=df[df.status.eq("OK")].copy()
base=ok.iloc[(ok.A_L-1.).abs().argsort()[:1]].iloc[0]
for c in ["chi2_high","chi2_lensing","chi2_lowT","chi2_lowE","chi2_cal","chi2_total"]:
    ok["d_"+c]=ok[c]-float(base[c])
ok=ok.sort_values("chi2_total")
ok.to_csv(OUT/"planck_lensing_shape_ranked.csv",index=False)
print("PLANCK_LENSING_SHAPE_BASE",base.to_dict(),flush=True)
print("PLANCK_LENSING_SHAPE_RANKED",ok.to_dict("records"),flush=True)
