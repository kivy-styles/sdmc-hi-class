#!/usr/bin/env python3
"""
A034/g022 fixed-background full-strength perturbative search.

The exact a034 distance background is frozen so Pantheon+, Union3, DES-Y5,
BAO, age and calibrated-ladder geometry are preserved. Search only the SDMC
structural/kinetic + primordial sector:
  As, ns, tau, A_F, z_c, width, D_floor, D0, power.
The aim is to deepen the already successful a034 Planck+DESI closure while
improving the scalar stability margin and, where possible, S8.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/a034_fullstrength_9d"); OUT.mkdir(parents=True,exist_ok=True)
T=2.7255; SIG=.0025; OR=4.17998772e-5

# a034/g022 background - freeze distance sector.
H0=69.71482083084993
OB=0.022083219194622913
OC=0.12299536722293603
LAM=18.40625
ZT=16.19317622240633
ALATE=0.024822032167576252
BLATE=0.01264704344701022

BASE=dict(As=2.1047019076951746e-09,ns=0.9628416801318527,tau=0.053167105823755265,
          AF=0.02048146490100771,zc=3.927876388467848,width=0.33114133956842123,
          Dfloor=0.05181542627513409,D0=0.34231919445927034,power=1.0)
KEYS=["As","ns","tau","AF","zc","width","Dfloor","D0","power"]
LOW=np.array([2.075e-9,0.9605,0.0470,0.0186,3.76,0.295,0.0470,0.30,0.72])
HIGH=np.array([2.125e-9,0.9660,0.0570,0.0224,4.08,0.375,0.0600,0.405,1.28])

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(T*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    h,s=load_cls(path)
    def f(A):
        v=float(high.chi_squared(h,A_planck=A))
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

def ini(root,c):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0}
omega_b={OB}
omega_cdm={OC}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={c['As']:.17e}
n_s={c['ns']:.17g}
tau_reio={c['tau']:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={c['AF']:.17g},{c['zc']:.17g},{c['width']:.17g},{c['D0']:.17g},{c['power']:.17g},{c['Dfloor']:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},0.5,{ALATE},0.25,{BLATE},1.5
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

def run(label,c):
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    r=dict(id=label,**{k:float(c[k]) for k in KEYS},status="FAIL")
    try:
        bg=table(root+"00_background.dat")
        r["D_min"]=float(bg["kin (D)"].min()); r["cs2_min"]=float(bg["c_s^2"].min()); r["cs2_max"]=float(bg["c_s^2"].max())
        r["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1
        if stable:
            pc,Ap=pscore(root+"00_cl_lensed.dat"); pks=sorted(OUT.glob(label+"_*pk.dat"))
            s=sig8(pks[0]); i0=int(np.argmin(np.abs(bg["z"].to_numpy()))); Om=float(bg.iloc[i0]["Omega_m(z)"])
            r.update(status="OK",planck_lite=pc,A_planck=Ap,sigma8=s,Omega_m0=Om,S8=s*math.sqrt(Om/0.3))
    except Exception as e:r["error"]=repr(e)
    print("A034_9D_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

rows=[run("a034",BASE)]
sob=qmc.Sobol(d=len(KEYS),scramble=True,seed=3402026)
for i,x in enumerate(qmc.scale(sob.random_base2(m=7),LOW,HIGH)):
    rows.append(run(f"q{i:03d}",{k:float(v) for k,v in zip(KEYS,x)}))

df=pd.DataFrame(rows)
base=df[df.id=="a034"].iloc[0]
df["d_planck"]=df["planck_lite"]-float(base.planck_lite)
df["d_S8"]=df["S8"]-float(base.S8)
df.to_csv(OUT/"a034_9d_scan.csv",index=False)
ok=df[(df.status=="OK")&(df.cs2_min>=0.006)].copy()
# Pareto rank Planck, S8 and stability; do not overfit a scalar score.
front=[]
for i,r in ok.iterrows():
    dom=((ok.planck_lite<=r.planck_lite)&(ok.S8<=r.S8)&(ok.cs2_min>=r.cs2_min)&
         ((ok.planck_lite<r.planck_lite)|(ok.S8<r.S8)|(ok.cs2_min>r.cs2_min))).any()
    if not dom:front.append(i)
pareto=ok.loc[front].sort_values(["planck_lite","S8"],ascending=[True,True])
prom=pd.concat([
    ok.sort_values(["planck_lite","S8"]).head(8),
    ok.sort_values(["S8","planck_lite"]).head(8),
    pareto.head(12)
]).drop_duplicates("id").head(20)
prom.to_csv(OUT/"a034_9d_promote.csv",index=False)
pareto.to_csv(OUT/"a034_9d_pareto.csv",index=False)
summary={"base":base.to_dict(),"n_ok":int(len(ok)),
         "promotion":prom.to_dict("records"),
         "rule":"Background fixed to a034; SN/BAO/age/SH0ES geometry preserved. Exact full Plik + official DESI required for promotion."}
(OUT/"a034_9d_summary.json").write_text(json.dumps(summary,indent=2))
print("A034_9D_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
