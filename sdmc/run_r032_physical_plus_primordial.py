#!/usr/bin/env python3
"""
Focused Planck closure screen at the best stable final-r032 physical lensing point.

Structural point held fixed to the best stable candidate found by the first
physical-lensing ridge screen:
  A_F = 0.066168... , z_c = 3.066441... , width = 0.563819...

The scan varies only (n_s, ln(1e10 A_s), tau_reio).  No phenomenological A_L
is used.  This stage asks whether the physical lensing repair and the
primordial/reionization calibration repair reinforce each other.

Screen only: finalists still require exact covariant reconstruction, full
Planck Plik nuisance profiling, raw DESI DR1 full-shape, and separate SN tests.
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

OUT=Path("output/physical_plus_primordial"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025

# Final-r032 expansion and ordinary cosmology.
H0=70.8514
OB=0.02239952
OC=0.12444227328918850
OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2
LAM=17.925
ZT=17.775
D0=0.34231919445927034
DFLOOR=0.045

# Best stable physical-lensing screen point (run 35535117956, sobol045).
AF=0.0661676420550793
ZC=3.0664408076554537
WIDTH=0.5638187950477004

NS_VALUES=[0.960,0.964,0.968,0.972]
TAU_VALUES=[0.0564,0.0584,0.0604,0.0624,0.0644]
LNAS_VALUES=[3.058,3.063,3.068,3.073]

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
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp},a

def score(path):
    harr,dls,a=load_cls(path)
    def pcs(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda A:pcs(A)[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5]),a

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]
    q=(k>1e-4)&(p>0)&np.isfinite(p); k=k[q]; p=p[q]
    x=8*k; W=np.ones_like(x); m=x!=0
    W[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    return float(np.sqrt(np.trapezoid(k**3*p*W**2/(2*np.pi**2),x=np.log(k))))

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

def ini(root,ns,lnAs,tau):
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
    n_s = {ns:.12g}
    tau_reio = {tau:.12g}

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF:.17g}, {ZC:.17g}, {WIDTH:.17g}, {D0:.17g}, 1.0, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g}, {LAM}, {ZT}, 0.5, 0.01105624999, 0.25, 0.01951933685, 1.5
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

lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,
                  stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
base,lcl=score(Path(lr+"00_cl_lensed.dat"))
lell=lcl[:,0].astype(int)
print("PHYSPRIM_LCDM",base,flush=True)

def ratio(a,col,L):
    e=a[:,0].astype(int)
    i=np.where(e==L)[0]; j=np.where(lell==L)[0]
    if not len(i) or not len(j): return np.nan
    return float(a[i[0],col]/lcl[j[0],col])

design=[
 ("physical_baseline",0.964,3.076,0.0544),
 ("diag_target",0.964,3.066,0.0584),
 ("diag_target_tilt",0.968,3.066,0.0584),
]
for ns in NS_VALUES:
    for tau in TAU_VALUES:
        for lnAs in LNAS_VALUES:
            design.append((f"n{ns:.3f}_t{tau:.4f}_A{lnAs:.3f}".replace(".","p"),ns,lnAs,tau))

rows=[]
for i,(tag,ns,lnAs,tau) in enumerate(design,1):
    root=str(OUT/(f"{i:03d}_{tag}_")); ip=OUT/(f"{i:03d}_{tag}.ini")
    ip.write_text(ini(root,ns,lnAs,tau))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec=dict(id=tag,index=i,A_F=AF,z_c=ZC,width=WIDTH,n_s=ns,ln10As=lnAs,
             A_s=math.exp(lnAs)/1e10,tau_reio=tau,Q=lnAs-2*tau,
             returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat"); pkp=Path(root+"00_pk.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()))
    if cp.returncode==0 and clp.exists() and pkp.exists():
        stable=rec.get("min_D",0)>0 and rec.get("min_cs2",0)>0 and rec.get("max_cs2",2)<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            sc,a=score(clp); rec.update(sc); rec["status"]="OK"; rec["sigma8"]=sigma8(pkp)
            for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
                rec["delta_"+k]=rec[k]-base[k]
            for L in [40,100,400,1000]: rec[f"pp_ratio_L{L}"]=ratio(a,5,L)
            for L in [200,500,1000,2000]: rec[f"tt_ratio_l{L}"]=ratio(a,1,L)
    if rec["status"]!="OK": rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    rows.append(rec)
    print("PHYSPRIM_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"physical_plus_primordial_profile.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
if ok.empty: raise SystemExit("No stable physical+primordial candidates")
keep=set(["physical_baseline","diag_target","diag_target_tilt"])
for col in ["chi2_total","chi2_high","chi2_lensing","chi2_cal","chi2_lowE"]:
    keep.update(ok.nsmallest(15,col)["id"].tolist())
short=ok[ok.id.isin(keep)].sort_values("chi2_total")
short.to_csv(OUT/"physical_plus_primordial_shortlist.csv",index=False)
print("PHYSPRIM_BEST_TOTAL",ok.nsmallest(20,"chi2_total").to_dict("records"),flush=True)
print("PHYSPRIM_BEST_HIGH",ok.nsmallest(15,"chi2_high").to_dict("records"),flush=True)
print("PHYSPRIM_BEST_LENS",ok.nsmallest(15,"chi2_lensing").to_dict("records"),flush=True)
print("PHYSPRIM_SHORTLIST",short.to_dict("records"),flush=True)
