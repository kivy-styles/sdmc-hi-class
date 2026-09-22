#!/usr/bin/env python3
"""
Bayesian structural refinement around SN-aware s110.

Late geometry is frozen at sobol223, so the three SN shape likelihoods are
unchanged. We optimize only the structural/No-Slip coordinates for Planck-lite,
with explicit stability and subluminality gates. The search box extends below
the previous Sobol boundaries in A_F and z_t because s110 sat close to both.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from skopt import gp_minimize
from skopt.space import Real
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure223_s110_bo_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# SN-improved late background: frozen.
H0=69.85635133907199
A_LATE=0.018589897081255913
B_LATE=0.013815921754576268
OB=0.022083219194622913
OC=0.12299536722293603
AS=2.117602412937069e-9
NS=0.9625227132590487
TAU=0.055202901571989066
D0=0.34231919445927034
LAM=18.40625
DNT=0.5
TAUA=0.25
TAUB=1.5

# Current exact-promotion point.
S110=[0.026252536563202738,3.79610941009596,0.3229390993900597,
      0.046828035046346486,15.982883202601224]

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
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,AF,ZC,W,DF,ZT):
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

rows=[]; counter=0
def evaluate(x,tag=None):
    global counter
    AF,ZC,W,DF,ZT=map(float,x)
    label=tag or f"bo{counter:03d}"; counter+=1
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,AF,ZC,W,DF,ZT))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=label,A_F=AF,z_c=ZC,width=W,D_floor=DF,z_t=ZT,status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    score=1.0e6
    if cp.returncode==0 and bgp.exists() and clp.exists():
        try:
            bg=table(bgp)
            rec["min_D"]=float(bg["kin (D)"].min())
            rec["min_cs2"]=float(bg["c_s^2"].min())
            rec["max_cs2"]=float(bg["c_s^2"].max())
            rec["max_abs_noslip"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
            stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1)
            rec["stable_subluminal"]=bool(stable)
            if stable:
                rec.update(pscore(clp)); rec["status"]="OK"; score=float(rec["chi2_planck"])
            else:
                # soft penalty carries useful boundary information to the GP
                score=5000.+1000.*max(0.,-rec["min_cs2"])+1000.*max(0.,rec["max_cs2"]-1.)
        except Exception as e:
            rec["error"]=repr(e)
    if rec["status"]!="OK":
        rec["error"]=rec.get("error",cp.stdout[-700:].replace("\n"," | "))
    rows.append(rec)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    print("S110BO_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return score

# Evaluate s110 exactly in this workflow first.
s110_score=evaluate(S110,"s110")

space=[
    Real(0.0180,0.0300,name="A_F"),
    Real(3.40,4.20,name="z_c"),
    Real(0.245,0.390,name="width"),
    Real(0.0390,0.0550,name="D_floor"),
    Real(15.35,16.20,name="z_t"),
]

# GP Bayesian optimization, seeded at s110.
res=gp_minimize(
    evaluate,space,
    x0=[S110],y0=[s110_score],
    n_calls=72,n_initial_points=18,
    acq_func="gp_hedge",noise=1e-8,
    random_state=94231,verbose=False,
)

df=pd.DataFrame(rows)
if "chi2_planck" in df:
    df["delta_vs_s110"]=df["chi2_planck"]-s110_score
df.to_csv(OUT/"snclosure223_s110_bo_refine.csv",index=False)
ok=df[df.status=="OK"].sort_values("chi2_planck")
best=ok.head(16).to_dict("records")
summary={
 "s110_score":float(s110_score),
 "s110":dict(zip(["A_F","z_c","width","D_floor","z_t"],S110)),
 "n_evaluated":int(len(df)),
 "n_ok":int(len(ok)),
 "best":best[0],
 "best_delta_vs_s110":float(best[0]["chi2_planck"]-s110_score),
 "bo_x":list(map(float,res.x)),
 "bo_fun":float(res.fun),
 "top16":best,
}
(OUT/"snclosure223_s110_bo_refine_summary.json").write_text(json.dumps(summary,indent=2))
print("S110BO_BEST",json.dumps(best,sort_keys=True),flush=True)
print("S110BO_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
