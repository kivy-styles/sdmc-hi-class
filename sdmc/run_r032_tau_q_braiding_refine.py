#!/usr/bin/env python3
from pathlib import Path
import json, math, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/tau_q_braid"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=69.7095092787552
OB=0.02211778594478965
OC=0.1258596659367904
OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2

AF=0.0715; ZC=5.; WIDTH=1.426; D0=0.34231919445927034; DFLOOR=.045
LAM=17.925; ZT=17.775
ALENS=.024; ZLENS=1.5; SIGLENS=.8

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=ell.max()+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); te=np.zeros(n); ee=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def pcs(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pcs(A)[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]; q=(k>1e-4)&(p>0)&np.isfinite(p)
    k=k[q]; p=p[q]; x=8*k; W=np.ones_like(x); m=x!=0
    W[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    return float(np.sqrt(np.trapezoid(k**3*p*W**2/(2*np.pi**2),x=np.log(k))))

def ini(root,ns,tau,Q):
    lnAs=Q+2*tau; As=math.exp(lnAs)/1e10
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {OB:.15g}
    omega_cdm = {OC:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {As:.17e}
    n_s = {ns:.15g}
    tau_reio = {tau:.15g}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_lens_braiding
    parameters_smg = {AF},{ZC},{WIDTH},{D0:.17g},1.0,{DFLOOR},{ALENS},{ZLENS},{SIGLENS}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
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

lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
base=score(Path(lr+"00_cl_lensed.dat"))
print("TAUQ_LCDM",base,flush=True)

taus=[.052,.055,.058,.061,.064]
Qs=[2.944,2.948,2.952,2.956]
nss=[.960,.964,.968]
design=[(0.06456347454525531,2.948208928124793,0.9637720508351921,"acoustic_best")]
for tau in taus:
    for Q in Qs:
        for ns in nss:
            design.append((tau,Q,ns,"grid"))

rows=[]
for i,(tau,Q,ns,kind) in enumerate(design,1):
    lnAs=Q+2*tau
    root=str(OUT/f"p{i:03d}_"); ip=OUT/f"p{i:03d}.ini"
    ip.write_text(ini(root,ns,tau,Q))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(index=i,kind=kind,tau_reio=tau,Q=Q,n_s=ns,ln10As=lnAs,
             A_s=math.exp(lnAs)/1e10,returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat"); pkp=Path(root+"00_pk.dat")
    if cp.returncode==0 and clp.exists() and pkp.exists():
        sc=score(clp); rec.update(sc); rec["sigma8"]=sigma8(pkp); rec["status"]="OK"
        rec["delta_planck"]=rec["chi2_total"]-base["chi2_total"]
    else:
        rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    rows.append(rec)
    print("TAUQ_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"tau_q_braiding_refine.csv",index=False)
ok=df[df.status=="OK"].copy()
print("TAUQ_BEST_TOTAL",ok.nsmallest(20,"chi2_total").to_dict("records"),flush=True)
print("TAUQ_BEST_LOWE",ok.nsmallest(20,"chi2_lowE").to_dict("records"),flush=True)
print("TAUQ_BEST_LENS",ok.nsmallest(20,"chi2_lensing").to_dict("records"),flush=True)
