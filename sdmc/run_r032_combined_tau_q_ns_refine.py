#!/usr/bin/env python3
"""
Local primordial/reionization refinement around the combined SDMC candidate.
Background, acoustic ruler, and localized braiding are fixed. We scan
(n_s, tau, Q), Q = ln(1e10 A_s) - 2 tau, to separate the primary-CMB
amplitude direction from reionization and tilt.
"""
from pathlib import Path
import json, math, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/tau_q_ns"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025
H0=69.7095092787552; OB=.02211778594478965; OC=.1258596659367904
Q0=2.948208928124793; NS0=.9637720508351921; TAU0=.06456347454525531
OX=.6953970160452582
AF=.0715; ZC=5.; WIDTH=1.426; D0=.34231919445927034; DFLOOR=.045
LAM=17.925; ZT=17.775
ALENS=.024; ZLENS=1.5; SLENS=.8

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=ell.max()+1; conv=(TCMB*1e6)**2
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
    op=minimize_scalar(lambda x:pieces(x)[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

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
    gravity_model = sdmc_v3_lens_braiding
    parameters_smg = {AF},{ZC},{WIDTH},{D0:.17g},1.0,{DFLOOR},{ALENS},{ZLENS},{SLENS}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
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

# Current point plus physically motivated anchors and deterministic Sobol design.
anchors=[
 ("current",NS0,TAU0,Q0),
 ("tau0544_constQ",NS0,.0544,Q0),
 ("tau056_constQ",NS0,.056,Q0),
 ("tau058_constQ",NS0,.058,Q0),
 ("tau060_constQ",NS0,.060,Q0),
 ("tau058_ns966",.966,.058,Q0),
 ("tau058_ns968",.968,.058,Q0),
]
LOW=np.array([.958,.052,Q0-.010])
HIGH=np.array([.970,.066,Q0+.010])
sob=qmc.Sobol(d=3,scramble=True,seed=20260920)
pts=qmc.scale(sob.random_base2(m=6),LOW,HIGH)
design=anchors+[(f"sobol{i+1:03d}",float(x[0]),float(x[1]),float(x[2])) for i,x in enumerate(pts)]

rows=[]
for i,(tag,ns,tau,Q) in enumerate(design,1):
    lnAs=Q+2*tau
    root=str(OUT/f"{i:03d}_{tag}_"); ip=OUT/f"{i:03d}_{tag}.ini"
    ip.write_text(ini(root,ns,tau,Q))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,index=i,n_s=ns,tau_reio=tau,Q=Q,ln10As=lnAs,
             A_s=math.exp(lnAs)/1e10,returncode=cp.returncode,status="FAIL")
    clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and clp.exists():
        rec.update(score(clp)); rec["status"]="OK"
    else:
        rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    rows.append(rec)
    print("TAUQNS_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
base=df[df.id=="current"].iloc[0]
for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
    df["delta_"+k+"_vs_current"]=df[k]-base[k]
df.to_csv(OUT/"tau_q_ns_refine.csv",index=False)
ok=df[df.status=="OK"].copy()
print("TAUQNS_CURRENT",base.to_dict(),flush=True)
print("TAUQNS_BEST_TOTAL",ok.nsmallest(15,"chi2_total").to_dict("records"),flush=True)
print("TAUQNS_BEST_LOWE",ok.nsmallest(15,"chi2_lowE").to_dict("records"),flush=True)
print("TAUQNS_BEST_HIGH",ok.nsmallest(15,"chi2_high").to_dict("records"),flush=True)
