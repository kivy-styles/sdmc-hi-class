#!/usr/bin/env python3
"""
Primordial closure around the exact-promoted lrref008 structural solution.

Hold the acoustically locked geometry and SDMC structural sector fixed.
Re-open only the primordial CMB trio using:
  n_s,
  tau,
  Q = ln(10^10 A_s) - 2 tau.

Q isolates the primary-amplitude combination while tau controls the amount
of A_s restored for lensing. This is a compact closure test, not a new broad
cosmology search.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/lrref008_primordial_closure"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0=70.7927624481251
OB=0.022083219194622913
OC=0.12299536722293603
NS0=0.9625227132590487
TAU0=0.055202901571989066
Q0=2.944463801247999
LNAS0=Q0+2.*TAU0

AF=0.035991880638059234
ZC=3.7172300264239313
WIDTH=0.4018391790613532
D0=0.34231919445927034
DF=0.045995755137875675
LAM=18.40625
ZT=16.24029149528593

# n_s, tau, Q
LOW=np.array([0.9560,0.0450,Q0-0.014],float)
HIGH=np.array([0.9700,0.0660,Q0+0.014],float)

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

def pscore(path):
    harr,dls=load_cls(path)
    def pc(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,ns,tau,lnAs):
    As=math.exp(lnAs)/1e10
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={As:.17e}
n_s={ns:.17g}
tau_reio={tau:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
output_background_smg=3
modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
write background=yes
root={root}
format=class
input_verbose=0
background_verbose=0
thermodynamics_verbose=0
perturbations_verbose=0
spectra_verbose=0
lensing_verbose=0
output_verbose=0
""")

def run(tag,ns,tau,Q):
    lnAs=Q+2.*tau
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini")
    ip.write_text(ini(root,ns,tau,lnAs))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=tag,n_s=ns,tau_reio=tau,Q=Q,ln10As=lnAs,A_s=math.exp(lnAs)/1e10,
             status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1)
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(pscore(clp)); rec["status"]="OK"
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    print("EDGE_PRIM_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[]
base=run("lowridge_center",NS0,TAU0,Q0)
if base["status"]!="OK": raise RuntimeError("zt-best center failed")
BASE=float(base["chi2_planck"])
rows.append(base)

# Axis checks guarantee direct turnover information.
for i,ns in enumerate([0.958,0.960,0.9615,0.9635,0.965,0.9675,0.969]):
    rows.append(run(f"ns_axis_{i:02d}",ns,TAU0,Q0))
for i,tau in enumerate([0.047,0.050,0.0525,0.0575,0.060,0.0635]):
    rows.append(run(f"tau_axis_{i:02d}",NS0,tau,Q0))
for i,dq in enumerate([-0.014,-0.009,-0.0045,0.0045,0.009,0.014]):
    rows.append(run(f"Q_axis_{i:02d}",NS0,TAU0,Q0+dq))

# 32 deterministic joint Sobol points.
sam=qmc.Sobol(d=3,scramble=False)
pts=qmc.scale(sam.random_base2(m=5),LOW,HIGH)
for i,p in enumerate(pts):
    ns,tau,Q=map(float,p)
    rows.append(run(f"sobol{i:03d}",ns,tau,Q))

for rec in rows:
    if rec.get("status")=="OK":
        rec["delta_planck_vs_lowridge"]=rec["chi2_planck"]-BASE
        rec["delta_high_vs_lowridge"]=rec["chi2_high"]-base["chi2_high"]
        rec["delta_lowT_vs_lowridge"]=rec["chi2_lowT"]-base["chi2_lowT"]
        rec["delta_lowE_vs_lowridge"]=rec["chi2_lowE"]-base["chi2_lowE"]
        rec["delta_lensing_vs_lowridge"]=rec["chi2_lensing"]-base["chi2_lensing"]
        rec["delta_cal_vs_lowridge"]=rec["chi2_cal"]-base["chi2_cal"]

df=pd.DataFrame(rows)
df.to_csv(OUT/"lrref008_primordial_closure.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
best=ok.nsmallest(12,"chi2_planck").to_dict("records")
summary={
  "center_chi2":BASE,
  "n_total":int(len(df)),
  "n_stable":int(len(ok)),
  "best":best[0] if best else None,
  "best_delta_vs_lowridge":float(ok.chi2_planck.min()-BASE) if len(ok) else None
}
(OUT/"lrref008_primordial_closure_summary.json").write_text(json.dumps(summary,indent=2))
print("EDGE_PRIM_BEST",json.dumps(best,sort_keys=True),flush=True)
print("EDGE_PRIM_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
