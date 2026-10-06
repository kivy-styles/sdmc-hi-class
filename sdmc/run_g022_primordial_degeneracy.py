#!/usr/bin/env python3
"""
g022 primordial-degeneracy search.

Keep the exact S117/g022 late background and structural sector frozen. Vary the
primordial amplitude and optical depth along the near-constant
A_s exp(-2 tau) direction, with a small n_s freedom. Goal: retain the g022
Planck-quality basin while reducing sigma8/S8 and therefore the compressed
weak-lensing penalty. Because the background is unchanged, SN, BAO, age and
calibrated-SH0ES background geometry are inherited from g022/S117.
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

OUT=Path("output/g022_primordial_degeneracy"); OUT.mkdir(parents=True,exist_ok=True)
T=2.7255; SIG=.0025
H0=69.71482083084993; OB=0.022083219194622913; OC=0.12299536722293603
AS0=2.117602412937069e-9; NS0=0.9625227132590487; TAU0=0.055202901571989066
AEFF0=AS0*math.exp(-2*TAU0)
AF=0.02048146490100771; ZC=3.927876388467848; WIDTH=0.33114133956842123
D0=0.34231919445927034; POWER=1.0; DF=0.05181542627513409
LAM=18.40625; ZT=16.19317622240633; ALATE=0.024822032167576252; BLATE=0.01264704344701022
OR=4.17998772e-5
BASE_CHI=None

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

def ini(root,As,ns,tau):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0}
omega_b={OB}
omega_cdm={OC}
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
parameters_smg={AF},{ZC},{WIDTH},{D0},{POWER},{DF}
expansion_model=sdmc_full
expansion_smg={ox},{LAM},{ZT},0.5,{ALATE},0.25,{BLATE},1.5
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

def run(label,tau,ns,amp_tilt):
    # amp_tilt permits ±0.4% departure from exact constant Aeff.
    As=AEFF0*math.exp(2*tau)*(1+amp_tilt)
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,As,ns,tau))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    r=dict(id=label,As=As,ns=ns,tau=tau,amp_tilt=amp_tilt,status="FAIL")
    try:
        bg=table(root+"00_background.dat")
        r["D_min"]=float(bg["kin (D)"].min()); r["cs2_min"]=float(bg["c_s^2"].min()); r["cs2_max"]=float(bg["c_s^2"].max())
        r["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1
        if stable:
            pc,Ap=pscore(root+"00_cl_lensed.dat"); pks=sorted(OUT.glob(label+"_*pk.dat"))
            s8=sig8(pks[0]); i0=int(np.argmin(np.abs(bg["z"].to_numpy()))); Om=float(bg.iloc[i0]["Omega_m(z)"])
            S8=s8*math.sqrt(Om/0.3)
            r.update(status="OK",chi2_planck_lite=pc,A_planck=Ap,sigma8=s8,Omega_m0=Om,S8=S8,
                     Aeff=As*math.exp(-2*tau))
    except Exception as e:r["error"]=repr(e)
    if r["status"]!="OK": r["error"]=r.get("error",cp.stdout[-500:].replace("\n"," | "))
    print("G022_PRIM_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

rows=[]
rows.append(run("center",TAU0,NS0,0.0))
sob=qmc.Sobol(d=3,scramble=True,seed=22026)
# tau 0.041-0.056, ns 0.9605-0.9685, Aeff tilt +/-0.4%
lo=np.array([0.0410,0.9605,-0.004]); hi=np.array([0.0560,0.9685,0.004])
for i,x in enumerate(qmc.scale(sob.random_base2(m=6),lo,hi)):
    rows.append(run(f"a{i:03d}",*map(float,x)))

df=pd.DataFrame(rows)
center=df[df.id=="center"].iloc[0]
df["delta_planck"]=df["chi2_planck_lite"]-float(center.chi2_planck_lite)
df["delta_S8"]=df["S8"]-float(center.S8)
df.to_csv(OUT/"g022_primordial_scan.csv",index=False)
ok=df[df.status=="OK"].copy()
# prefer at least 0.005 S8 improvement and <=0.35 native-lite cost; if empty relax.
prom=ok[(ok.delta_planck<=0.35)&(ok.delta_S8<=-0.004)].sort_values(["S8","delta_planck"]).head(12)
if len(prom)<4:
    prom=ok[ok.delta_planck<=0.6].sort_values(["S8","delta_planck"]).head(12)
prom.to_csv(OUT/"g022_primordial_promote.csv",index=False)
summary={"center":center.to_dict(),"promotion":prom.to_dict("records"),
         "rule":"Candidates require exact full-Plik and official DESI promotion before replacing exact g022."}
(OUT/"g022_primordial_summary.json").write_text(json.dumps(summary,indent=2))
print("G022_PRIM_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
