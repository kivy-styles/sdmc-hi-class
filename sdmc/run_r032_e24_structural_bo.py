#!/usr/bin/env python3
"""
Bayesian structural re-closure on the original edge024-Q late background.

Goal:
  find a stable perturbative/kinetic completion that improves the CMB enough
  to preserve the strong edge024-Q raw-DESI score while retaining its fixed
  SN distance curve.

Search variables:
  AF, zc, width, D_floor, Q
where Q = ln(1e10 A_s) - 2 tau.

This is a screening optimizer using Planck native-lite + low-l TT/EE + lensing.
Top points must be promoted to exact covariant full-Plik + raw DESI before any
claim of closure.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from skopt import Optimizer
from skopt.space import Real
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/e24_structural_bo"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Frozen edge024-Q late background and ordinary densities.
H0=70.79187943335671
OB=0.022083219194622913
OC=0.12299536722293603
NS=0.9625227132590487
TAU=0.055202901571989066
AL=0.01105624999
BL=0.01951933685
OX=0.7104246827081668
LAM=18.40625
ZT=16.173189924377947
DN=.5; TAUA=.25; TAUB=1.5
D0=0.34231919445927034

# Original exact successful edge024-Q center.
CENTER=[0.03200255395658314,4.007449422683567,0.34799921004101636,
        0.04607192634791136,2.942463801247999]

space=[
 Real(0.024,0.040,name="AF"),
 Real(3.65,4.35,name="ZC"),
 Real(0.285,0.430,name="width"),
 Real(0.042,0.064,name="D_floor"),
 Real(2.9395,2.9465,name="Q"),
]

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
    a=np.loadtxt(path); e=a[:,0].astype(int); n=int(e.max())+1
    cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[e]=a[:,1]*cv; ee[e]=a[:,2]*cv; te[e]=a[:,3]*cv
    L=e.astype(float); pp[e]=a[:,5]*L*(L+1.)
    return np.column_stack([e,tt[e],te[e],ee[e]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    harr,d=load_cls(path)
    def pc(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(d,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    p=pc(float(op.x))
    return dict(A_planck=float(op.x),chi2_high=p[0],chi2_lowT=p[1],
                chi2_lowE=p[2],chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,AF,ZC,W,DF,Q):
    AS=math.exp(Q+2*TAU)/1e10
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
expansion_smg={OX:.17g},{LAM},{ZT},{DN},{AL},{TAUA},{BL},{TAUB}
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

counter=0
rows=[]
def evaluate(x,tag=None):
    global counter
    AF,ZC,W,DF,Q=map(float,x)
    if tag is None: tag=f"bo{counter:03d}"
    counter+=1
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini")
    ip.write_text(ini(root,AF,ZC,W,DF,Q))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=300)
    rec=dict(id=tag,AF=AF,ZC=ZC,width=W,D_floor=DF,Q=Q,
             A_s=math.exp(Q+2*TAU)/1e10,status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    objective=1.0e6
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.0
        rec["stable_subluminal"]=bool(stable)
        # Require a small positive robustness floor rather than accepting a
        # numerically marginal cs2 crossing.
        robust=stable and rec["min_cs2"]>=0.002
        rec["robust_stability"]=bool(robust)
        if robust:
            rec.update(pscore(clp)); rec["status"]="OK"
            objective=rec["chi2_planck"]
            # very small soft penalty discourages razor-thin stability floors
            if rec["min_cs2"]<0.01:
                objective += 0.05*(0.01-rec["min_cs2"])/0.008
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-700:].replace("\n"," | ")
    rec["objective"]=float(objective)
    rows.append(rec)
    print("E24BO_POINT",json.dumps(rec,sort_keys=True),flush=True)
    for p in OUT.glob(tag+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return float(objective)

# Evaluate the known successful center exactly in the same screen.
y0=evaluate(CENTER,"center")
opt=Optimizer(space,base_estimator="GP",acq_func="gp_hedge",
              n_initial_points=12,initial_point_generator="sobol",
              random_state=240924,acq_optimizer="sampling")
opt.tell(CENTER,y0)

# 38 additional BO evaluations. Batch of 2 retains some exploration while
# allowing GP updates frequently.
for it in range(19):
    xs=opt.ask(n_points=2,strategy="cl_min")
    ys=[]
    for j,x in enumerate(xs):
        ys.append(evaluate(x,f"bo{2*it+j:03d}"))
    opt.tell(xs,ys)
    ok=[r for r in rows if r.get("status")=="OK"]
    best=min(ok,key=lambda r:r["chi2_planck"]) if ok else None
    print("E24BO_PROGRESS",it,json.dumps(best,sort_keys=True) if best else "NONE",flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"e24_structural_bo.csv",index=False)
ok=df[df.status.eq("OK")].sort_values(["chi2_planck","min_cs2"],ascending=[True,False])
best=ok.head(12).to_dict("records")
summary={
 "fixed_background":{"H0":H0,"A":AL,"B":BL,"OmegaX":OX},
 "center":rows[0],
 "n_evaluated":len(rows),
 "n_stable":int(len(ok)),
 "best":best,
 "note":"Top points require exact covariant full-Plik + raw DESI promotion. SN shape is fixed to edge024-Q during this search."
}
(OUT/"e24_structural_bo_summary.json").write_text(json.dumps(summary,indent=2))
print("E24BO_BEST",json.dumps(best,sort_keys=True),flush=True)
print("E24BO_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
