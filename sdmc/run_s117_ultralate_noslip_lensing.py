#!/usr/bin/env python3
"""
S117-background ultra-late narrow No-Slip lensing ridge.

Purpose:
  Search a much wider Planck-mass transition domain than the local S117 scans,
  while freezing the exact S117 ordinary cosmology and late background that
  already close Planck+DESI+all three SN alternatives.

This tests a concrete hypothesis suggested by the earlier physical-lensing
workflows: weak-lensing pressure may be relieved by a delayed/broader No-Slip
transition rather than by lowering A_s. Only A_F, z_c and width vary here.

Any point is only a screening candidate until exact full Plik + official DESI
full shape are promoted. Because the background is frozen to S117, SN, BAO,
age and calibrated-ladder geometry are inherited at screening level.
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

OUT=Path("output/s117_ultralate_lensing"); OUT.mkdir(parents=True,exist_ok=True)
T=2.7255; SIG=.0025; OR=4.17998772e-5
H0=69.71482083084993
OB=0.022083219194622913
OC=0.12299536722293603
AS=2.117602412937069e-9
NS=0.9625227132590487
TAU=0.055202901571989066
D0=0.34231919445927034
DF=0.052829781055450435
POWER=1.0
LAM=18.40625
ZT=16.19317622240633
DNT=0.5
ALATE=0.024822032167576252
BLATE=0.01264704344701022
S117_LITE=1020.917005539309
S117_S8=0.8408471470723361

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(T*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    ha,s=load_cls(path)
    def f(A):
        v=float(high.chi_squared(ha,A_planck=A))
        v+=float(-2*lowT.log_likelihood(s["tt"],calib=A))
        v+=float(-2*lowE.log_likelihood(s["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        v+=float(-2*lens.log_likelihood(s,**lp))
        v+=((A-1)/SIG)**2
        return v
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def sig8(path):
    a=np.loadtxt(path); k=a[:,0]; P=a[:,1]; x=8*k
    W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
    return float(np.sqrt(np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))))

def ini(root,AF,zc,width):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={AS:.17e}
n_s={NS:.17g}
tau_reio={TAU:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={AF:.17g},{zc:.17g},{width:.17g},{D0:.17g},{POWER:.17g},{DF:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM:.17g},{ZT:.17g},{DNT:.17g},{ALATE:.17g},0.25,{BLATE:.17g},1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
modes=s
output=tCl,pCl,lCl,mPk
lensing=yes
l_max_scalars=3000
P_k_max_h/Mpc=2.5
z_pk=0
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

def run(label,AF,zc,width):
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,AF,zc,width))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    r=dict(id=label,A_F=float(AF),z_c=float(zc),width=float(width),status="FAIL",returncode=cp.returncode)
    try:
        bg=table(root+"00_background.dat")
        r["D_min"]=float(bg["kin (D)"].min()); r["cs2_min"]=float(bg["c_s^2"].min()); r["cs2_max"]=float(bg["c_s^2"].max())
        r["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        r["F0"]=float(np.interp(0,bg.z,bg["M*^2_smg"]))
        stable=cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1
        r["stable_subluminal"]=bool(stable)
        if stable:
            pc,Ap=pscore(root+"00_cl_lensed.dat")
            pks=sorted(OUT.glob(label+"_*pk.dat")); s8=sig8(pks[0])
            i0=int(np.argmin(np.abs(bg.z.to_numpy()))); Om=float(bg.iloc[i0]["Omega_m(z)"]); S8=s8*math.sqrt(Om/0.3)
            r.update(status="OK",chi2_planck_lite=pc,A_planck=Ap,delta_planck=pc-S117_LITE,
                     sigma8=s8,Omega_m0=Om,S8=S8,delta_S8=S8-S117_S8,
                     DESY3_chi2=((S8-0.776)/0.017)**2,
                     HSC_chi2=((S8-0.805)/0.018)**2,
                     KiDS_chi2=((S8-0.815)/0.020)**2)
    except Exception as e:r["error"]=repr(e)
    if r["status"]!="OK":r["error"]=r.get("error",cp.stdout[-700:].replace("\n"," | "))
    print("S117_ULTRALATE_LENS_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

anchors=[
 ("S117",0.02114389337040484,3.715736017236486,0.35696151830255984),
 ("ultraA",0.060,0.50,0.12),
 ("ultraB",0.100,0.50,0.12),
 ("ultraC",0.080,0.30,0.08),
 ("ultraD",0.120,0.75,0.18),
 ("ultraE",0.150,0.35,0.10),
]
rows=[run(*x) for x in anchors]
LOW=np.array([0.020,0.15,0.05]); HIGH=np.array([0.160,1.20,0.35])
sob=qmc.Sobol(d=3,scramble=True,seed=117920)
for i,x in enumerate(qmc.scale(sob.random_base2(m=6),LOW,HIGH)):
    rows.append(run(f"l{i:03d}",*map(float,x)))

df=pd.DataFrame(rows); df.to_csv(OUT/"s117_ultralate_lensing.csv",index=False)
ok=df[df.status=="OK"].copy()
front=[]
for i,r in ok.iterrows():
    dom=((ok.chi2_planck_lite<=r.chi2_planck_lite)&(ok.S8<=r.S8)&
         ((ok.chi2_planck_lite<r.chi2_planck_lite)|(ok.S8<r.S8))).any()
    if not dom:front.append(i)
pareto=ok.loc[front].sort_values(["S8","chi2_planck_lite"])
pareto.to_csv(OUT/"s117_ultralate_lensing_pareto.csv",index=False)
# Three promotion channels: Planck-safe, WL-improved, and compromise.
safe=ok[ok.chi2_planck_lite<=S117_LITE+0.5].sort_values(["S8","chi2_planck_lite"]).head(8)
mid=ok[ok.chi2_planck_lite<=S117_LITE+2.0].sort_values(["S8","chi2_planck_lite"]).head(8)
lowS=ok.sort_values(["S8","chi2_planck_lite"]).head(8)
prom=pd.concat([safe,mid,lowS]).drop_duplicates("id")
prom.to_csv(OUT/"s117_ultralate_lensing_promote.csv",index=False)
summary={"n_ok":int(len(ok)),"pareto":pareto.head(30).to_dict("records"),
         "promotion":prom.to_dict("records"),
         "rule":"Exact full Plik + official DESI are required before a promoted point can replace S117."}
(OUT/"s117_ultralate_lensing_summary.json").write_text(json.dumps(summary,indent=2))
print("S117_ULTRALATE_LENS_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
