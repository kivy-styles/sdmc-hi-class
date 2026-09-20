#!/usr/bin/env python3
"""
Focused primordial/reionization refinement around the exact physical No-Slip
structure that reached Delta chi2_full-Plik = +0.60195.

Structural sector is fixed:
  AF=0.0525, zc=4.6, width=0.680590909090909
  Dfloor=0.045, D0=0.34231919445927034
  lambda_e=17.925, z_t=17.775

Ordinary background is fixed to the current physical seed; only
(n_s, tau, Q=ln(1e10 A_s)-2 tau) are scanned.  Q is used because primary
high-l CMB amplitude approximately follows A_s exp(-2 tau).
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/physical_best_primordial"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.66173072; OB=0.02224; OC=0.12444227328918850
OX=0.7061451696512845
AF=0.0525; ZC=4.6; WIDTH=0.680590909090909
D0=0.34231919445927034; DFLOOR=0.045
LAMBDA=17.925; ZT=17.775

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
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def ini(ns,tau,Q,root):
    lnAs=Q+2*tau
    As=math.exp(lnAs)/1e10
    return textwrap.dedent(f"""\
    H0 = {H0}
    omega_b = {OB}
    omega_cdm = {OC}
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
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF},{ZC},{WIDTH},{D0:.17g},1.0,{DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX},{LAMBDA},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic

    modes = s
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
    output_background_smg = 3
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

lr=str(OUT/"lcdm_"); ip=OUT/"lcdm.ini"; ip.write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2500:])
lcdm=score(Path(lr+"00_cl_lensed.dat"))
print("PRIM_REFINE_LCDM",lcdm,flush=True)

NS=[0.960,0.962,0.964,0.966,0.968]
TAUS=[0.052,0.055,0.058,0.061,0.064]
QS=[2.944,2.948,2.952,2.956,2.960]
design=[]
# exact current physical seed first
design.append(("seed",0.964,0.061,3.074-2*0.061))
for ns in NS:
  for tau in TAUS:
    for Q in QS:
      design.append((f"n{ns:.3f}_t{tau:.3f}_q{Q:.3f}",ns,tau,Q))

rows=[]
for i,(tag,ns,tau,Q) in enumerate(design,1):
    lnAs=Q+2*tau
    root=str(OUT/(f"{i:03d}_{tag}_")); ip=OUT/(f"{i:03d}_{tag}.ini")
    ip.write_text(ini(ns,tau,Q,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,n_s=ns,tau_reio=tau,Q=Q,ln10As=lnAs,A_s=math.exp(lnAs)/1e10,
             returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()))
    if cp.returncode==0 and clp.exists():
        stable=rec.get("min_D",-1)>0 and rec.get("min_cs2",-1)>0 and rec.get("max_cs2",2)<=1.
        rec["stable_subluminal"]=bool(stable)
        if stable:
            sc=score(clp); rec.update(sc); rec["status"]="OK"
            for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
                rec["delta_"+k]=sc[k]-lcdm[k]
    if rec["status"]!="OK": rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    rows.append(rec)
    print("PRIM_REFINE_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"physical_best_primordial_refine.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
print("PRIM_REFINE_BEST_TOTAL",ok.nsmallest(15,"chi2_total").to_dict("records"),flush=True)
print("PRIM_REFINE_BEST_HIGHELL",ok.nsmallest(15,"chi2_high").to_dict("records"),flush=True)
print("PRIM_REFINE_BEST_LOWE",ok.nsmallest(15,"chi2_lowE").to_dict("records"),flush=True)
print("PRIM_REFINE_BEST_LENSING",ok.nsmallest(15,"chi2_lensing").to_dict("records"),flush=True)
