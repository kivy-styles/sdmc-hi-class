#!/usr/bin/env python3
"""
Final-r032 physical lensing-ridge scan.

Vary only the No-Slip Planck-mass transition (A_F, z_c, width) at the exact
r032 expansion/background point and current ordinary cosmology. No artificial
A_L is used. The goal is to test whether a stable trajectory can keep F(a)
closer to unity through the CMB-lensing kernel while preserving the primary
CMB fit.

This is a Planck screening diagnostic. Any finalist must later be rebuilt as
an exact covariant replay and tested with full Plik, raw DESI DR1 full shape,
and each SN compilation separately.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/physical_lensing_ridge"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025

# Exact r032 expansion / current ordinary cosmology.
H0=70.8514
OB=0.02239952
OC=0.12444227328918850
NS=0.964
LNAS=3.076
TAU=0.0544
AS=math.exp(LNAS)/1e10
OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2
LAM=17.925
ZT=17.775
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
    a=np.loadtxt(path)
    ell=a[:,0].astype(int); n=int(ell.max())+1
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
    op=minimize_scalar(lambda x:pcs(x)[-1],bounds=(0.97,1.03),
                       method="bounded",options={"xatol":1e-10})
    A=float(op.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5]),a

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]
    q=(k>1e-4)&(p>0)&np.isfinite(p); k=k[q]; p=p[q]
    x=8*k
    W=np.ones_like(x); m=x!=0
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

def sdmc_ini(root,AF,zc,width):
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
    parameters_smg = {AF:.17g}, {zc:.17g}, {width:.17g}, {D0:.17g}, 1.0, {DFLOOR}
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

# LCDM baseline.
lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,
                  stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
base,lcl=score(Path(lr+"00_cl_lensed.dat"))
lell=lcl[:,0].astype(int)
print("PHYS_LENS_LCDM",base,flush=True)

def ratio(a,col,L):
    e=a[:,0].astype(int)
    i=np.where(e==L)[0]; j=np.where(lell==L)[0]
    if not len(i) or not len(j): return np.nan
    return float(a[i[0],col]/lcl[j[0],col])

# Target current + physically motivated anchors + Sobol.
anchors=[
    ("current",(0.0715,5.0,1.426)),
    ("later_sameA",(0.0715,2.5,1.20)),
    ("later_narrow",(0.0715,2.0,0.80)),
    ("lowerA_later",(0.0650,2.5,1.00)),
    ("lowerA_laterwide",(0.0600,2.0,1.40)),
    ("mildlower_later",(0.0680,3.0,1.20)),
    ("higherA_later",(0.0750,2.5,1.20)),
    ("historical_joint",(0.07310,4.90,1.430)),
    ("historical_micro",(0.073128,4.84,1.426)),
    ("historical_edge",(0.07360,4.25,1.430)),
]
LOW=np.array([0.050,1.0,0.45])
HIGH=np.array([0.080,8.0,2.00])
sob=qmc.Sobol(d=3,scramble=True,seed=20260920)
pts=qmc.scale(sob.random_base2(m=6),LOW,HIGH)
design=anchors+[(f"sobol{i+1:03d}",tuple(x)) for i,x in enumerate(pts)]

rows=[]
for i,(tag,(AF,zc,width)) in enumerate(design,1):
    root=str(OUT/(f"{i:03d}_{tag}_")); ip=OUT/(f"{i:03d}_{tag}.ini")
    ip.write_text(sdmc_ini(root,float(AF),float(zc),float(width)))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec=dict(id=tag,index=i,A_F=float(AF),z_c=float(zc),width=float(width),
             returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat")
    clp=Path(root+"00_cl_lensed.dat")
    pkp=Path(root+"00_pk.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
                   F0=float(np.interp(0,bg.z,bg["M*^2_smg"])))
        for z in [0.5,1,2,3,5]:
            rec[f"F_z{str(z).replace('.','p')}"]=float(np.interp(z,bg.z,bg["M*^2_smg"]))
    if cp.returncode==0 and clp.exists() and pkp.exists():
        stable=rec.get("min_D",0)>0 and rec.get("min_cs2",0)>0 and rec.get("max_cs2",2)<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            sc,a=score(clp); rec.update(sc); rec["status"]="OK"; rec["sigma8"]=sigma8(pkp)
            for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
                rec["delta_"+k]=rec[k]-base[k]
            for L in [40,100,400,1000]:
                rec[f"pp_ratio_L{L}"]=ratio(a,5,L)
            for L in [200,500,1000,2000]:
                rec[f"tt_ratio_l{L}"]=ratio(a,1,L)
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    rows.append(rec)
    print("PHYS_LENS_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"physical_lensing_ridge.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
if ok.empty: raise SystemExit("No stable physical lensing candidates")
# Keep multiple objectives because the best Planck total need not be the best
# physical lensing candidate after DESI is included.
keep=set(["current"])
for col in ["chi2_total","chi2_high","chi2_lensing","chi2_cal"]:
    keep.update(ok.nsmallest(12,col)["id"].tolist())
# Prefer candidates that restore pp while staying near LCDM high-l.
ok["pp_mismatch"]=abs(ok["pp_ratio_L100"]-1)+abs(ok["pp_ratio_L400"]-1)+abs(ok["pp_ratio_L1000"]-1)
keep.update(ok.nsmallest(12,"pp_mismatch")["id"].tolist())
short=ok[ok.id.isin(keep)].copy()
short["screen_rank"]=short[["delta_chi2_total","delta_chi2_high","delta_chi2_lensing"]].sum(axis=1)
short=short.sort_values("delta_chi2_total")
short.to_csv(OUT/"physical_lensing_shortlist.csv",index=False)

print("PHYS_LENS_BEST_TOTAL",ok.nsmallest(15,"chi2_total").to_dict("records"),flush=True)
print("PHYS_LENS_BEST_LENS",ok.nsmallest(15,"chi2_lensing").to_dict("records"),flush=True)
print("PHYS_LENS_BEST_PP",ok.nsmallest(15,"pp_mismatch").to_dict("records"),flush=True)
print("PHYS_LENS_SHORTLIST",short.to_dict("records"),flush=True)
