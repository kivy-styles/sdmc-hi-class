#!/usr/bin/env python3
"""
S117 amplitude-compensated growth search.

Purpose
-------
The exact S117 background is frozen. This preserves the already-tested SN,
compressed BAO, BBN, cosmic-age and calibrated-distance-ladder geometry.
We then search the primordial-amplitude / reionization / structural-kinetic
subspace for lower S8 while retaining the S117 Planck quality.

The scan varies:
  A_s, n_s, tau_reio,
  A_F, z_c, width, D0, power, D_floor.

This directly targets the remaining downstream weakness identified by Part IX:
high S8 / compressed weak-lensing pressure.

The nonlinear JWST mapper coefficient g3 remains a separate UV coefficient.
It is not used in this linear Boltzmann search, so no JWST tuning is allowed to
back-react on CMB, DESI-linear or weak-lensing observables.
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

OUT=Path("output/s117_ampcomp")
OUT.mkdir(parents=True,exist_ok=True)

TCMB=2.7255
CAL_SIG=.0025
OR=4.17998772e-5

# Frozen exact S117 background/ordinary geometry.
H0=69.71482083084993
OB=0.022083219194622913
OC=0.12299536722293603
LAM=18.40625
ZT=16.19317622240633
DNT=0.5
ALATE=0.024822032167576252
BLATE=0.01264704344701022

# Recomputed native-lite center in the current workflow.
S117_PLANCK=1020.917005539309
S117_S8=0.8404896714875673

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
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
    a=np.loadtxt(path)
    ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float)
    pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    ha,s=load_cls(path)
    def f(A):
        v=float(high.chi_squared(ha,A_planck=A))
        v+=float(-2*lowT.log_likelihood(s["tt"],calib=A))
        v+=float(-2*lowE.log_likelihood(s["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        v+=float(-2*lens.log_likelihood(s,**lp))
        v+=((A-1.)/CAL_SIG)**2
        return v
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; P=a[:,1]
    x=8*k
    W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
    return float(np.sqrt(np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))))

def ini(root,p):
    h=H0/100.
    ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={p['As']:.17e}
n_s={p['ns']:.17g}
tau_reio={p['tau']:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={p['AF']:.17g},{p['zc']:.17g},{p['width']:.17g},{p['D0']:.17g},{p['power']:.17g},{p['Dfloor']:.17g}
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

def run(label,p):
    root=str(OUT/(label+"_"))
    ip=OUT/(label+".ini")
    ip.write_text(ini(root,p))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=420)
    r=dict(id=label,**p,status="FAIL",returncode=cp.returncode)
    try:
        bg=table(root+"00_background.dat")
        r["D_min"]=float(bg["kin (D)"].min())
        r["cs2_min"]=float(bg["c_s^2"].min())
        r["cs2_max"]=float(bg["c_s^2"].max())
        r["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=(cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1)
        r["stable_subluminal"]=bool(stable)
        if stable:
            pc,Ap=pscore(root+"00_cl_lensed.dat")
            pks=sorted(OUT.glob(label+"_*pk.dat"))
            s8=sigma8(pks[0])
            i0=int(np.argmin(np.abs(bg["z"].to_numpy())))
            Om=float(bg.iloc[i0]["Omega_m(z)"])
            S8=s8*math.sqrt(Om/0.3)
            r.update(status="OK",chi2_planck_lite=pc,A_planck=Ap,
                     delta_planck=pc-S117_PLANCK,sigma8=s8,Omega_m0=Om,S8=S8,
                     delta_S8=S8-S117_S8)
            # Objective is only a screen: exact full-Plik decides promotion.
            # Penalize Planck loss strongly but reward meaningful S8 relief.
            r["objective"]=float(max(0,pc-S117_PLANCK-0.15)+
                                 0.35*((S8-0.815)/0.015)**2)
    except Exception as e:
        r["error"]=repr(e)
    if r["status"]!="OK":
        r["error"]=r.get("error",cp.stdout[-700:].replace("\n"," | "))
    print("S117_AMPCOMP_POINT",json.dumps(r,sort_keys=True),flush=True)
    for q in OUT.glob(label+"_*"):
        try:q.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

# Degeneracy-guided seeds: preserve approximately A_s exp(-2 tau).
As0=2.117602412937069e-9
tau0=0.055202901571989066
def tau_amp(As):
    return tau0+0.5*math.log(As/As0)

base_struct=dict(AF=0.02114389337040484,zc=3.715736017236486,width=0.35696151830255984,
                 D0=0.34231919445927034,power=1.0,Dfloor=0.052829781055450435)
seeds=[
 ("S117",dict(As=As0,ns=0.9625227132590487,tau=tau0,**base_struct)),
 ("amp206",dict(As=2.06e-9,ns=0.9630,tau=tau_amp(2.06e-9),**base_struct)),
 ("amp204",dict(As=2.04e-9,ns=0.9640,tau=tau_amp(2.04e-9),**base_struct)),
 ("amp202",dict(As=2.02e-9,ns=0.9650,tau=tau_amp(2.02e-9),**base_struct)),
 ("amp200",dict(As=2.00e-9,ns=0.9660,tau=max(0.025,tau_amp(2.00e-9)),**base_struct)),
 ("amp202g073",dict(As=2.02e-9,ns=0.9650,tau=tau_amp(2.02e-9),
                      AF=0.020284758737310768,zc=3.6509963892400266,width=0.35059972247108817,
                      D0=0.34231919445927034,power=0.45,Dfloor=0.05064729745686054)),
]
rows=[run(label,p) for label,p in seeds]

# 256-point 9D Sobol search around the amplitude/reionization degeneracy and
# the viable S117 kinetic basin.
LOW=np.array([1.985e-9,0.9580,0.0250,0.0185,3.45,0.285,0.300,0.30,0.0440])
HIGH=np.array([2.125e-9,0.9700,0.0585,0.0245,4.10,0.425,0.430,0.90,0.0620])
KEYS=["As","ns","tau","AF","zc","width","D0","power","Dfloor"]
sob=qmc.Sobol(d=9,scramble=True,seed=202610061)
cloud=qmc.scale(sob.random_base2(m=8),LOW,HIGH)
for i,x in enumerate(cloud):
    p={k:float(v) for k,v in zip(KEYS,x)}
    rows.append(run(f"a{i:03d}",p))

df=pd.DataFrame(rows)
df.to_csv(OUT/"s117_ampcomp_scan.csv",index=False)
ok=df[df.status=="OK"].copy()

front=[]
for i,r in ok.iterrows():
    dominated=((ok.chi2_planck_lite<=r.chi2_planck_lite)&(ok.S8<=r.S8)&
               ((ok.chi2_planck_lite<r.chi2_planck_lite)|(ok.S8<r.S8))).any()
    if not dominated:front.append(i)
pareto=ok.loc[front].sort_values(["S8","chi2_planck_lite"])
pareto.to_csv(OUT/"s117_ampcomp_pareto.csv",index=False)

# Three promotion tiers.
tight=ok[ok.chi2_planck_lite<=S117_PLANCK+0.25].sort_values(["S8","chi2_planck_lite"]).head(10)
moderate=ok[ok.chi2_planck_lite<=S117_PLANCK+0.75].sort_values(["S8","chi2_planck_lite"]).head(10)
allbest=ok.sort_values(["objective","S8","chi2_planck_lite"]).head(12)
for name,d in [("tight",tight),("moderate",moderate),("objective",allbest)]:
    d.to_csv(OUT/f"s117_ampcomp_{name}_promote.csv",index=False)

# Independent nonlinear JWST UV target from the preceding derivation.
g3_95=0.01573642223*(0.120/0.054)
g3_68=0.01573642223*(0.152/0.054)
summary={
 "n_total":len(rows),"n_ok":int(len(ok)),
 "S117":{"planck_native_lite":S117_PLANCK,"S8":S117_S8},
 "best_tight":tight.to_dict("records"),
 "best_moderate":moderate.to_dict("records"),
 "best_objective":allbest.to_dict("records"),
 "pareto":pareto.head(30).to_dict("records"),
 "JWST_nonlinear_completion":{"g3_95_equiv_epsilon20":g3_95,
                               "g3_68_equiv_epsilon20":g3_68,
                               "linear_sector_effect":"none by construction of nonlinear gate",
                               "status":"conditional UV coefficient until action-level normalization is derived"},
 "promotion_rule":"No point is a final winner until full Plik, official DESI DR1 full shape, exact 3-SN, calibrated SH0ES and downstream battery are rerun."
}
(OUT/"s117_ampcomp_summary.json").write_text(json.dumps(summary,indent=2))
print("S117_AMPCOMP_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
