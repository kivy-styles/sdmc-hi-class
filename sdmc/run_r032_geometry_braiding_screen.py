#!/usr/bin/env python3
"""
Local positive-braiding screen on the current geometry-compromise leader.

Background/ordinary cosmology is held fixed at t=0.25:
H0=70.6661848725751, wb=0.022329086486197414,
wc=0.12479662145108897, ns=0.963943012708798,
tau=0.0544, ln(1e10 As)=3.060.
Structural No-Slip envelope: AF=0.05709640212077647,
zc=3.0902165911160413, width=0.5979282500222326.

Only Delta_B(z)=A exp[-(z-zcB)^2/(2 sigma_B^2)] is varied.
D is held fixed by the sdmc_v3_lens_braiding implementation.
This is an EFT diagnostic screen, not yet a covariant G3 result.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/sobol049_braid_screen"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.5659273196888
OB=0.02213832986768335
OC=0.12419066331610083
NS=0.9614044621046632
TAU=0.05464880801877007
LNAS=3.051477308589965
AS=math.exp(LNAS)/1e10
OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2
AF=0.0573794336710125
ZC=2.8532181778922676
WIDTH=0.575869657304138
D0=0.34231919445927034
DFLOOR=0.045

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
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def pieces(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pieces(A)[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def ini(root,A,zl,sl):
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {OB:.15g}
    omega_cdm = {OC:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {AS:.17e}
    n_s = {NS:.15g}
    tau_reio = {TAU:.15g}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_lens_braiding
    parameters_smg = {AF:.17g},{ZC:.17g},{WIDTH:.17g},{D0:.17g},1.0,{DFLOOR:.17g},{A:.17g},{zl:.17g},{sl:.17g}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},17.925,17.775,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    a_ini_over_a_today_default = 1.e-8
    a_ini_test_qs_smg = 1.e-8
    pert_ic_ini_z_ref_smg = 1.e7
    a_min_stability_test_smg = 1.e-8
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

design=[
    (0.000,1.50,0.80),
    (0.032,1.50,0.80),
    (0.036,1.50,0.80),
    (0.040,1.50,0.80),
    (0.044,1.50,0.80),
    (0.048,1.50,0.80),
    (0.052,1.50,0.80),
    (0.056,1.50,0.80),
    (0.032,1.50,1.00),
    (0.040,1.50,1.00),
]

rows=[]
for i,(A,zl,sl) in enumerate(design,1):
    tag=f"a{A:.3f}_z{zl:.2f}_s{sl:.2f}".replace(".","p")
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini")
    ip.write_text(ini(root,A,zl,sl))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,index=i,A_lens=A,z_lens=zl,sigma_lens=sl,returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_slip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    if cp.returncode==0 and bgp.exists() and clp.exists():
        rec["stable_subluminal"]=bool(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.)
        if rec["stable_subluminal"]:
            rec.update(score(clp)); rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    rows.append(rec)
    print("S49_BRAID_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
zero=df[(df.status=="OK")&(df.A_lens==0.)].iloc[0]
for k in ["chi2_total","chi2_high","chi2_lensing","chi2_lowT","chi2_lowE","chi2_cal"]:
    df["delta_"+k+"_vs_zero"]=df[k]-zero[k]
df.to_csv(OUT/"sobol049_braiding_screen.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
best=ok.nsmallest(15,"chi2_total")
best.to_csv(OUT/"sobol049_braiding_shortlist.csv",index=False)
print("S49_BRAID_ZERO",json.dumps(zero.to_dict(),default=float,sort_keys=True),flush=True)
print("S49_BRAID_BEST",json.dumps(best.to_dict("records"),default=float,sort_keys=True),flush=True)
