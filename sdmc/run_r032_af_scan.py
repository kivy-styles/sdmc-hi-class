#!/usr/bin/env python3
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/af_scan"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.8514; OB=0.02239952; OC=0.12444227328918850
NS=0.964; AS=2.16715426219737e-9; TAU=0.0544
OX=0.7073985893
ZC=5.00; WIDTH=1.426; D0=0.34231919445927034; DFLOOR=0.045
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def parts(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    o=minimize_scalar(lambda A:parts(A)[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(o.x); return (A,)+parts(A)

def ini(AF,root):
    return textwrap.dedent(f"""\
    H0 = {H0}
    omega_b = {OB}
    omega_cdm = {OC}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {AS:.17e}
    n_s = {NS}
    tau_reio = {TAU}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF:.15g}, {ZC}, {WIDTH}, {D0:.17g}, 1.0, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX}, 17.925, 17.775, 0.5, 0.01105624999, 0.25, 0.01951933685, 1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
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
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
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

lr=str(OUT/"lcdm_"); lip=OUT/"lcdm.ini"; lip.write_text(lcdm(lr))
cp=subprocess.run(["./class",str(lip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
ls=score(lr+"00_cl_lensed.dat")
print("AFSCAN_LCDM",ls,flush=True)

AFS=[0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.0715,0.08,0.09,0.10]
rows=[]
for i,af in enumerate(AFS):
    root=str(OUT/f"a{i:02d}_"); ip=OUT/f"a{i:02d}.ini"; ip.write_text(ini(af,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"AF":af,"returncode":cp.returncode,"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
        bg=table(bgp); z=bg["z"].to_numpy(); iz=np.argmin(np.abs(z))
        rec.update(F0=float(bg["M*^2_smg"].to_numpy()[iz]),
                   min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    if cp.returncode==0 and clp.exists():
        ps=score(clp)
        rec.update(status="OK",A_planck=ps[0],chi2_high=ps[1],chi2_lowT=ps[2],chi2_lowE=ps[3],
                   chi2_lensing=ps[4],chi2_cal=ps[5],chi2_planck=ps[6],delta_planck=ps[6]-ls[-1])
    else:
        rec["error"]=cp.stdout[-1500:].replace("\n"," | ")
    print("AFSCAN_POINT",json.dumps(rec,sort_keys=True),flush=True); rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"af_scan.csv",index=False)
ok=df[df.status=="OK"].copy()
print("AFSCAN_BEST_TOTAL",ok.nsmallest(len(ok),"delta_planck").to_dict("records"),flush=True)
print("AFSCAN_BEST_LENS",ok.nsmallest(len(ok),"chi2_lensing").to_dict("records"),flush=True)
print("AFSCAN_BEST_CAL",ok.nsmallest(len(ok),"chi2_cal").to_dict("records"),flush=True)
