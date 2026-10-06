#!/usr/bin/env python3
"""
Nine-dimensional fixed-background CMB/perturbation refinement around fj014.

The entire fj014 distance background is frozen:
H0, omega_b, omega_cdm, z_t, A_late, B_late.
Therefore its Pantheon+, Union3, DES-Y5, background BAO geometry, age and
calibrated-ladder geometry are preserved at screening level. Only primordial
and structural/kinetic coordinates are varied to remove the remaining Planck
gap without paying another SN penalty.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/fj014_fixedbg_cmb_refine"); OUT.mkdir(parents=True,exist_ok=True)
T=2.7255; SIG=.0025; OR=4.17998772e-5
H0=69.51135345487855; OB=0.0221336738740094; OC=0.12258112205294892
ZT=16.06170386472717; ALATE=0.011254649444762618; BLATE=0.00994506031619385
LAM=18.40625

CENTER=np.array([2.1096255943551658e-9,0.9639264293648303,0.05461631424725056,
                 0.01948851037286222,3.644913007412106,0.35893273874651643,
                 0.3865188211016357,1.0877414136193693,0.05569390742853284])
LOW=np.array([2.095e-9,0.9620,0.0485,0.0184,3.52,0.325,0.345,0.92,0.0500])
HIGH=np.array([2.122e-9,0.9665,0.0580,0.0204,3.78,0.400,0.435,1.22,0.0620])
KEYS=["As","ns","tau","AF","zc","width","D0","power","Dfloor"]

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)
def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(T*1e6)**2
    tt=np.zeros(n);ee=np.zeros(n);te=np.zeros(n);pp=np.zeros(n)
    tt[ell]=a[:,1]*conv;ee[ell]=a[:,2]*conv;te[ell]=a[:,3]*conv;L=ell.astype(float);pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}
def pscore(path):
    h,s=load_cls(path)
    def f(Ap):
        v=float(high.chi_squared(h,A_planck=Ap))
        v+=float(-2*lowT.log_likelihood(s["tt"],calib=Ap))
        v+=float(-2*lowE.log_likelihood(s["ee"],calib=Ap))
        lp={lens.calibration_param:Ap} if getattr(lens,"calibration_param",None) else {}
        v+=float(-2*lens.log_likelihood(s,**lp)); v+=((Ap-1)/SIG)**2
        return v
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)
def sig8(path):
    a=np.loadtxt(path);k=a[:,0];P=a[:,1];x=8*k;W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
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
def run(label,x):
    c={k:float(v) for k,v in zip(KEYS,x)};root=str(OUT/(label+"_"));ip=OUT/(label+".ini");ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    r=dict(id=label,**c,status="FAIL")
    try:
        bg=tab(root+"00_background.dat");r["D_min"]=float(bg["kin (D)"].min());r["cs2_min"]=float(bg["c_s^2"].min());r["cs2_max"]=float(bg["c_s^2"].max())
        r["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1
        if stable:
            pc,Ap=pscore(root+"00_cl_lensed.dat");pks=sorted(OUT.glob(label+"_*pk.dat"));s=sig8(pks[0]);i0=int(np.argmin(np.abs(bg["z"].to_numpy())));Om=float(bg.iloc[i0]["Omega_m(z)"])
            r.update(status="OK",planck_lite=pc,A_planck=Ap,sigma8=s,Omega_m0=Om,S8=s*math.sqrt(Om/0.3))
    except Exception as e:r["error"]=repr(e)
    print("FJ014_FIXEDBG_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

rows=[run("center",CENTER)]
sob=qmc.Sobol(d=len(KEYS),scramble=True,seed=14014)
for i,x in enumerate(qmc.scale(sob.random_base2(m=7),LOW,HIGH)):
    rows.append(run(f"r{i:03d}",x))
df=pd.DataFrame(rows)
base=float(df[df.id=="center"].iloc[0].planck_lite)
df["d_planck_vs_fj014"]=df["planck_lite"]-base
df.to_csv(OUT/"fj014_fixedbg_scan.csv",index=False)
ok=df[(df.status=="OK")&(df.cs2_min>=0.006)].sort_values(["planck_lite","S8"]).head(12)
ok.to_csv(OUT/"fj014_fixedbg_promote.csv",index=False)
summary={"center_planck_lite":base,"promotion":ok.to_dict("records"),
         "background":"fj014 fixed; SN/BAO-distance/age/SH0ES geometry unchanged at screen level",
         "rule":"Promote exact full Plik first; then official DESI full shape and exact three-SN replay."}
(OUT/"fj014_fixedbg_summary.json").write_text(json.dumps(summary,indent=2))
print("FJ014_FIXEDBG_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
