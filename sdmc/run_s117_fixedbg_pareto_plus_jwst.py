#!/usr/bin/env python3
"""
Second-stage S117 fixed-background perturbation search.

The S117 ordinary cosmological and late-background geometry is frozen so that
its exact SN/BAO/age/SH0ES behavior is preserved. Only the SDMC perturbation /
independent-kinetic structural coordinates are varied. The aim is to reduce
S8/weak-lensing pressure without giving away the Planck closure.

JWST is deliberately separated from the linear sector: the nonlinear mapper
coefficient g3 is an independent UV coefficient at this stage. The workflow
therefore reports the g3 needed for the FRESCO epsilon~0.20 target but does not
use g3 to alter CMB, DESI-linear, SN, BAO or weak-lensing observables.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/s117_fixedbg_pareto"); OUT.mkdir(parents=True,exist_ok=True)
T=2.7255; SIG=.0025
H0=69.71482083084993
OB=0.022083219194622913
OC=0.12299536722293603
AS=2.117602412937069e-9
NS=0.9625227132590487
TAU=0.055202901571989066
LAM=18.40625
ZT=16.19317622240633
DNT=0.5
ALATE=0.024822032167576252
BLATE=0.01264704344701022
OR=4.17998772e-5
S117_LITE=1020.8312884546846
S117_S8=0.8401099567317601

# S117 original and the later g073 point.
SEEDS=[
 ("S117",0.02114389337040484,3.715736017236486,0.35696151830255984,0.34231919445927034,1.0,0.052829781055450435),
 ("g073",0.020284758737310768,3.6509963892400266,0.35059972247108817,0.34231919445927034,0.45,0.05064729745686054),
]

# Wider than the old 4D scan; adds D0 and power, which later Lagrangian-guided
# work showed can change the kinetic response while leaving the S117 geometry fixed.
LOW=np.array([0.0185,3.45,0.285,0.300,0.30,0.044])
HIGH=np.array([0.0245,4.10,0.425,0.430,0.85,0.062])

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

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

def ini(root,AF,ZC,W,D0,POWER,DF):
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
parameters_smg={AF:.17g},{ZC:.17g},{W:.17g},{D0:.17g},{POWER:.17g},{DF:.17g}
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

def run(label,*pars):
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,*pars))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    AF,ZC,W,D0,POWER,DF=pars
    r=dict(id=label,A_F=AF,z_c=ZC,width=W,D0=D0,power=POWER,D_floor=DF,status="FAIL")
    try:
        bg=table(root+"00_background.dat")
        r["D_min"]=float(bg["kin (D)"].min()); r["cs2_min"]=float(bg["c_s^2"].min()); r["cs2_max"]=float(bg["c_s^2"].max())
        r["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1
        r["stable_subluminal"]=bool(stable)
        if stable:
            pc,Ap=pscore(root+"00_cl_lensed.dat")
            pks=sorted(OUT.glob(label+"_*pk.dat"))
            s8=sig8(pks[0]); i0=int(np.argmin(np.abs(bg["z"].to_numpy()))); Om=float(bg.iloc[i0]["Omega_m(z)"])
            S8=s8*math.sqrt(Om/0.3)
            r.update(status="OK",chi2_planck_lite=pc,A_planck=Ap,delta_planck=pc-S117_LITE,
                     sigma8=s8,Omega_m0=Om,S8=S8,delta_S8=S8-S117_S8)
    except Exception as e:r["error"]=repr(e)
    if r["status"]!="OK":r["error"]=r.get("error",cp.stdout[-500:].replace("\n"," | "))
    print("S117_FIXEDBG_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

rows=[run(*s) for s in SEEDS]
sob=qmc.Sobol(d=6,scramble=True,seed=1172026)
for i,x in enumerate(qmc.scale(sob.random_base2(m=7),LOW,HIGH)):
    rows.append(run(f"p{i:03d}",*map(float,x)))

df=pd.DataFrame(rows); df.to_csv(OUT/"s117_fixedbg_pareto.csv",index=False)
ok=df[df.status=="OK"].copy()
front=[]
for i,r in ok.iterrows():
    dom=((ok.chi2_planck_lite<=r.chi2_planck_lite)&(ok.S8<=r.S8)&
         ((ok.chi2_planck_lite<r.chi2_planck_lite)|(ok.S8<r.S8))).any()
    if not dom:front.append(i)
pareto=ok.loc[front].sort_values(["S8","chi2_planck_lite"])
pareto.to_csv(OUT/"s117_fixedbg_pareto_front.csv",index=False)
near=ok[ok.chi2_planck_lite<=S117_LITE+0.50].sort_values(["S8","chi2_planck_lite"]).head(12)
near.to_csv(OUT/"s117_fixedbg_promote.csv",index=False)

# Nonlinear JWST target is orthogonal at this stage.
g3_95=0.01573642223*(0.120/0.054)
g3_68=0.01573642223*(0.152/0.054)
summary={
 "n_ok":int(len(ok)),
 "S117":{"chi2_planck_lite":S117_LITE,"S8":S117_S8},
 "pareto":pareto.head(30).to_dict("records"),
 "promotion":near.to_dict("records"),
 "JWST_independent_nonlinear_target":{"g3_for_epsilon20_95_equiv":g3_95,
                                     "g3_for_epsilon20_68_equiv":g3_68,
                                     "note":"candidate UV coefficient; no linear-CMB effect assumed"}
}
(OUT/"s117_fixedbg_summary.json").write_text(json.dumps(summary,indent=2))
print("S117_FIXEDBG_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
