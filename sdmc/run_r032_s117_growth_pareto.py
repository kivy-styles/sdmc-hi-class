#!/usr/bin/env python3
"""
S117 structural growth Pareto scan.

The exact S117 late background and ordinary cosmological sector are frozen.
Only four perturbative / independent-kinetic structural coordinates vary:
  A_F, z_c, width, D_floor.

For each point:
  * exact hi_class stability / subluminality / No-Slip diagnostics,
  * Planck 2018 native-lite TTTEEE+lowT+lowE+lensing score,
  * linear z=0 matter spectrum,
  * sigma8 and S8.

No weak-lensing likelihood is folded into the Planck score.  The output is a
Pareto map: points that lower S8 while retaining the S117 Planck quality are
shortlisted for exact full-Plik + raw-DESI promotion.
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

OUT=Path("output/s117_g073_structural_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Exact S117 cosmological/background sector.
H0=69.71482083084993
OB=0.022083219194622913
OC=0.12299536722293603
AS=2.117602412937069e-9
NS=0.9625227132590487
TAU=0.055202901571989066
D0=0.34231919445927034
LAM=18.40625
ZT=16.19317622240633
DNT=0.5
A_LATE=0.024822032167576252
TAUA=0.25
B_LATE=0.01264704344701022
TAUB=1.5

# Exact S117 structural center.
AF0=0.020284758737310768
ZC0=3.6509963892400266
W0=0.35059972247108817
DF0=0.05064729745686054

S117_LITE=1020.658201861972
S117_S8=0.840968957

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
        x=float(high.chi_squared(harr,A_planck=A))
        x+=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        x+=float(-2*lens.log_likelihood(dls,**lp))
        x+=((A-1.)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(pc,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def sigma8_from_pk(path):
    a=np.loadtxt(path)
    k=a[:,0]; P=a[:,1]
    x=8.0*k
    W=np.where(np.abs(x)<1e-5,1-x*x/10.0,3*(np.sin(x)-x*np.cos(x))/x**3)
    sig2=np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))
    return float(np.sqrt(sig2))

def ini(root,AF,ZC,W,DF):
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
parameters_smg={AF:.17g},{ZC:.17g},{W:.17g},{D0:.17g},1.0,{DF:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM:.17g},{ZT:.17g},{DNT},{A_LATE:.17g},{TAUA},{B_LATE:.17g},{TAUB}
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

def run(label,AF,ZC,W,DF):
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini")
    ip.write_text(ini(root,AF,ZC,W,DF))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    rec=dict(id=label,A_F=float(AF),z_c=float(ZC),width=float(W),D_floor=float(DF),
             status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    pks=sorted(OUT.glob(label+"_*pk.dat"))
    if cp.returncode==0 and bgp.exists() and clp.exists() and pks:
        try:
            bg=table(bgp)
            rec["min_D"]=float(bg["kin (D)"].min())
            rec["min_cs2"]=float(bg["c_s^2"].min())
            rec["max_cs2"]=float(bg["c_s^2"].max())
            rec["max_abs_noslip"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
            stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
            rec["stable_subluminal"]=bool(stable)
            if stable:
                pc,Ap=pscore(clp)
                sig8=sigma8_from_pk(pks[0])
                i0=int(np.argmin(np.abs(bg["z"].to_numpy())))
                Om0=float(bg.iloc[i0]["Omega_m(z)"])
                S8=sig8*math.sqrt(Om0/0.3)
                rec.update(status="OK",chi2_planck=pc,A_planck=Ap,sigma8=sig8,
                           Omega_m0=Om0,S8=S8,
                           delta_planck_vs_s117=pc-S117_LITE,
                           delta_S8_vs_s117=S8-S117_S8)
        except Exception as e:
            rec["error"]=repr(e)
    if rec["status"]!="OK":
        rec["error"]=rec.get("error",cp.stdout[-700:].replace("\n"," | "))
    print("S117_G073_REFINE_POINT",json.dumps(rec,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return rec

# Center + the strongest neighboring structures from the preceding 128-point rescue.
known=[
 ("g073",AF0,ZC0,W0,DF0),
 ("g109",0.020244,3.748068,0.344728,0.048076),
 ("g022",0.02048146490100771,3.927876388467848,0.33114133956842123,0.05181542627513409)
]
# Replace abbreviated two entries with exact values only if read from deterministic scan is unnecessary;
# the Sobol cloud below covers them. Keep the first three exact known seeds.
rows=[]
for item in known:
    rows.append(run(*item))

# 128-point local Sobol map around S117.
low=np.array([0.0185,3.48,0.325,0.045])
highb=np.array([0.0212,3.84,0.378,0.056])
sob=qmc.Sobol(d=4,scramble=True,seed=25092026)
for i,x in enumerate(qmc.scale(sob.random_base2(m=6),low,highb)):
    rows.append(run(f"g{i:03d}",*map(float,x)))

df=pd.DataFrame(rows)
df.to_csv(OUT/"s117_g073_refine.csv",index=False)
ok=df[df.status=="OK"].copy()

# Non-dominated Pareto frontier in (Planck chi2, S8).
front=[]
for idx,r in ok.iterrows():
    dominated=((ok.chi2_planck<=r.chi2_planck)&(ok.S8<=r.S8)&
               ((ok.chi2_planck<r.chi2_planck)|(ok.S8<r.S8))).any()
    if not dominated: front.append(idx)
pareto=ok.loc[front].sort_values(["chi2_planck","S8"])
pareto.to_csv(OUT/"s117_g073_refine_front.csv",index=False)

# Promotion sets: preserve Planck quality first, then look for growth relief.
near=ok[ok.chi2_planck<=S117_LITE+0.08].sort_values(["S8","chi2_planck"]).head(15)
buffer=ok[ok.chi2_planck<S117_LITE].sort_values(["S8","chi2_planck"]).head(15)
summary={
 "g073_center":{"chi2_planck_lite":S117_LITE,"S8":S117_S8,
                "A_F":AF0,"z_c":ZC0,"width":W0,"D_floor":DF0},
 "n_ok":int(len(ok)),
 "pareto":pareto.head(30).to_dict("records"),
 "best_S8_within_plus_0p08_planck":near.to_dict("records"),
 "planck_better_than_s117":buffer.to_dict("records"),
 "promotion_rule":"No point supersedes S117 until exact full-Plik and raw-DESI are rerun; fixed S117 background preserves the SN/SH0ES geometry."
}
(OUT/"s117_g073_refine_summary.json").write_text(json.dumps(summary,indent=2))
print("S117_G073_REFINE_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
