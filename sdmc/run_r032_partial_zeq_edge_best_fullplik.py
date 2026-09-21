#!/usr/bin/env python3
from pathlib import Path
import json, math, subprocess, textwrap
import numpy as np
import pandas as pd
from iminuit import Minuit
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE import TTTEEE
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/partial_zeq_edge_best_fullplik")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255

# Improved acoustic geometry.
H0=69.7095092787552
OB=0.02211778594478965
OC=0.1258596659367904
NS=0.9637720508351921

# Best stable physical No-Slip Planck-mass ridge.
AF=0.0573794336710125
ZC=3.0902165911160413
WIDTH=0.5979282500222326
D0=0.34231919445927034
DFLOOR=0.045
LAM=17.7
ZT=17.1
OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2

CANDIDATES=[
    dict(id="partial_zeq_edge_best",H0=70.66175429589923,ob=0.02206374637655471,
         oc=0.12328183150909171,ns=0.9625227132590487,tau=0.055202901571989066,
         lnAs=3.0598696043919773,AF=0.0529375,
         zc=2.89375,width=0.54625,dfloor=0.03825,
         lam=18.475,zt=16.962500000000002),
]
def table(path):
    import re
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def ini(root,c):
    H0=float(c["H0"]); ob=float(c["ob"]); oc=float(c["oc"]); ns=float(c["ns"])
    tau=float(c["tau"]); lnAs=float(c["lnAs"])
    As=math.exp(lnAs)/1e10
    h=H0/100.; OX=1.-(ob+oc+OR)/(h*h)
    return textwrap.dedent(f"""\\
    H0 = {H0:.15g}
    omega_b = {ob:.15g}
    omega_cdm = {oc:.15g}
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
    parameters_smg = {float(c["AF"]):.17g},{float(c["zc"]):.17g},{float(c["width"]):.17g},{D0:.17g},1.0,{float(c["dfloor"]):.17g}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{float(c["lam"]):.17g},{float(c["zt"]):.17g},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3

    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
    write background = yes
    write thermodynamics = yes
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

meta={}
for c in CANDIDATES:
    root=str(OUT/(c["id"]+"_"))
    ip=OUT/(c["id"]+".ini")
    ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    if cp.returncode:
        print("PZEQEDGEBEST_CLASS_FAIL",c["id"],cp.stdout[-3000:],flush=True)
        raise SystemExit(f"CLASS failed {c['id']}")
    bg=table(root+"00_background.dat")
    rec=dict(c)
    rec.update(
        min_D=float(bg["kin (D)"].min()),
        min_cs2=float(bg["c_s^2"].min()),
        max_cs2=float(bg["c_s^2"].max()),
        max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
        F0=float(bg.iloc[np.argmin(np.abs(bg.z.to_numpy()))]["M*^2_smg"]),
    )
    rec["stable_subluminal"]=bool(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.)
    print("PZEQEDGEBEST_STABILITY",json.dumps(rec,sort_keys=True),flush=True)
    if not rec["stable_subluminal"]:
        raise SystemExit(f"stability failed {c['id']}")
    meta[c["id"]]=rec

lcdm_root=str(OUT/"lcdm_")
lcdm_ini=OUT/"lcdm.ini"
lcdm_ini.write_text(textwrap.dedent(f"""\
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
output = tCl,pCl,lCl
modes = s
lensing = yes
l_max_scalars = 3000
root = {lcdm_root}
format = class
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
spectra_verbose = 0
lensing_verbose = 0
output_verbose = 0
"""))
cp=subprocess.run(["./class",str(lcdm_ini)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
if cp.returncode:
    raise RuntimeError(cp.stdout[-3000:])

high=TTTEEE(packages_path="planck_full")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

fixed={
  "cib_index":-1.3,"galf_TE_index":-2.4,"galf_EE_index":-2.4,
  "A_sbpx_100_100_TT":1.0,"A_sbpx_143_143_TT":1.0,"A_sbpx_143_217_TT":1.0,"A_sbpx_217_217_TT":1.0,
  "galf_EE_A_100":0.055,"galf_EE_A_100_143":0.040,"galf_EE_A_100_217":0.094,
  "galf_EE_A_143":0.086,"galf_EE_A_143_217":0.21,"galf_EE_A_217":0.70,
  "A_cnoise_e2e_100_100_EE":1.0,"A_cnoise_e2e_143_143_EE":1.0,"A_cnoise_e2e_217_217_EE":1.0,
  "A_sbpx_100_100_EE":1.0,"A_sbpx_100_143_EE":1.0,"A_sbpx_100_217_EE":1.0,
  "A_sbpx_143_143_EE":1.0,"A_sbpx_143_217_EE":1.0,"A_sbpx_217_217_EE":1.0,
  "A_pol":1.0,"calib_100P":1.021,"calib_143P":0.966,"calib_217P":1.040,
}
names=[
  "A_planck","calib_100T","calib_217T","A_cib_217","xi_sz_cib","A_sz","ksz_norm",
  "gal545_A_100","gal545_A_143","gal545_A_143_217","gal545_A_217",
  "ps_A_100_100","ps_A_143_143","ps_A_143_217","ps_A_217_217",
  "galf_TE_A_100","galf_TE_A_100_143","galf_TE_A_100_217","galf_TE_A_143","galf_TE_A_143_217","galf_TE_A_217"
]
default=np.array([1.0,1.0002,0.99805,67.,0.05,7.,3.,8.6,10.6,23.5,91.9,257.,47.,40.,104.,.130,.130,.46,.207,.69,1.938],float)
old_cov=np.array([
1.0045213554451686,0.9998424064177515,0.9982224709218035,45.72423851559507,0.7248916197434251,6.976523105384968,5.500029149071257e-06,
8.535071574604354,10.780908169898211,20.060917984547913,96.15277043474971,254.14725387485439,53.74852187444169,57.054180993875775,125.69753555422906,
0.11540298701803288,0.13448306014088762,0.48055099721034644,0.2266942865821831,0.66858903191407,2.1183252356550293],float)
bounds=[
(.97,1.03),(.9946,1.0058),(.99285,1.00325),(0,200),(0,1),(0,10),(0,10),
(0,24.6),(0,26.6),(0,91.5),(0,251.9),(0,400),(0,400),(0,400),(0,400),
(0,0.466),(0,0.418),(0,1.18),(0,0.783),(0,1.41),(0,6.258)]
gauss={
"A_planck":(1.,.0025),"calib_100T":(1.0002,.0007),"calib_217T":(.99805,.00065),
"gal545_A_100":(8.6,2.),"gal545_A_143":(10.6,2.),"gal545_A_143_217":(23.5,8.5),"gal545_A_217":(91.9,20.),
"galf_TE_A_100":(.130,.042),"galf_TE_A_100_143":(.130,.036),"galf_TE_A_100_217":(.46,.09),
"galf_TE_A_143":(.207,.072),"galf_TE_A_143_217":(.69,.09),"galf_TE_A_217":(1.938,.54)}

if set(high.expected_params)!=(set(fixed)|set(names)):
    raise RuntimeError("full Plik parameter set mismatch")

def load(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    dl={k:np.zeros(n) for k in ["tt","ee","bb","te","pp","tp","ep"]}
    dl["tt"][ell]=a[:,1]*conv; dl["ee"][ell]=a[:,2]*conv; dl["te"][ell]=a[:,3]*conv; dl["bb"][ell]=a[:,4]*conv
    ll=ell.astype(float)*(ell.astype(float)+1.)
    dl["pp"][ell]=a[:,5]*ll
    dl["tp"][ell]=a[:,6]*np.sqrt(ll)*(TCMB*1e6)
    dl["ep"][ell]=a[:,7]*np.sqrt(ll)*(TCMB*1e6)
    cl={k:np.zeros(n) for k in ["tt","ee","bb","te"]}
    fac=np.zeros(n); q=ell>=2; fac[ell[q]]=2*np.pi/ll[q]
    for k in cl: cl[k][ell]=dl[k][ell]*fac[ell]
    return cl,dl

def prior_chi2(p):
    v=sum(((p[n]-mu)/sig)**2 for n,(mu,sig) in gauss.items())
    v+=((p["ksz_norm"]+1.6*p["A_sz"]-9.5)/3.)**2
    return float(v)

def make_objective(path):
    cl,dl=load(path)
    def fun(*x):
        p=dict(fixed); p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p)); A=p["A_planck"]
        lt=float(lowT.log_likelihood(dl["tt"],calib=A))
        le=float(lowE.log_likelihood(dl["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        if not all(np.isfinite([lh,lt,le,ll])): return 1e100
        return float(-2*(lh+lt+le+ll)+prior_chi2(p))
    def pieces(x):
        p=dict(fixed); p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p)); A=p["A_planck"]
        lt=float(lowT.log_likelihood(dl["tt"],calib=A))
        le=float(lowE.log_likelihood(dl["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        return p,lh,lt,le,ll
    return fun,pieces

def run_one(label,path,seeds):
    fun,pieces=make_objective(path)
    best=None
    for i,seed in enumerate(seeds):
        m=Minuit(fun,*seed,name=names); m.errordef=1.; m.strategy=1; m.tol=.1
        for n,b in zip(names,bounds): m.limits[n]=b
        for n in names: m.errors[n]=max(1e-5,.02*max(abs(m.values[n]),1.))
        m.migrad(ncall=7000)
        if not m.fmin.is_valid:
            m.simplex(ncall=3000); m.migrad(ncall=7000)
        row=(float(m.fval),np.array([m.values[n] for n in names]),bool(m.fmin.is_valid),int(m.nfcn))
        print("PZEQEDGEBEST_FULLPLIK_SEED",label,i,row[0],row[2],row[3],flush=True)
        if best is None or row[0]<best[0]: best=row
    fval,x,valid,nfcn=best
    p,lh,lt,le,ll=pieces(x)
    row=dict(model=label,chi2_profile_total=fval,chi2_high_full=-2*lh,chi2_lowT=-2*lt,
             chi2_lowE=-2*le,chi2_lensing=-2*ll,chi2_nuisance_priors=prior_chi2(p),
             A_planck=p["A_planck"],valid=valid,nfcn=nfcn)
    row.update({n:p[n] for n in names if n!="A_planck"})
    print("PZEQEDGEBEST_FULLPLIK_RESULT",json.dumps(row,sort_keys=True),flush=True)
    return row,x

lcdm,lx=run_one("LCDM",OUT/"lcdm_00_cl_lensed.dat",[default,old_cov])
rows=[lcdm]
prev=lx
for c in CANDIDATES:
    row,cx=run_one(c["id"],OUT/(c["id"]+"_00_cl_lensed.dat"),[prev,lx,default,old_cov])
    prev=cx
    row.update(meta[c["id"]])
    row["delta_chi2_vs_lcdm"]=row["chi2_profile_total"]-lcdm["chi2_profile_total"]
    rows.append(row)
    print("PZEQEDGEBEST_FULLPLIK_DELTA",c["id"],row["delta_chi2_vs_lcdm"],flush=True)

pd.DataFrame(rows).to_csv(OUT/"partial_zeq_edge_best_fullplik.csv",index=False)
best=min(rows[1:],key=lambda r:r["chi2_profile_total"])
summary=dict(lcdm_chi2=lcdm["chi2_profile_total"],best_model=best["model"],
             best_delta_vs_lcdm=best["chi2_profile_total"]-lcdm["chi2_profile_total"])
(OUT/"partial_zeq_edge_best_fullplik_summary.json").write_text(json.dumps(summary,indent=2))
print("PZEQEDGEBEST_FULLPLIK_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
