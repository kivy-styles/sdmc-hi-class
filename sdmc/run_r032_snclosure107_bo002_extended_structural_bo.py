#!/usr/bin/env python3
"""
Six-dimensional Bayesian structural refinement around exact late-local006 + bo002.

Frozen low-z/SN background: late-local006
  H0=69.71275747716427
  A_late=0.0013601897284388
  B_late=0.0147392870020121

Reopened structural coordinates:
  A_F, z_c, width, D_floor, D0, z_t

Rationale:
  The preceding 4-D BO reached a stable optimum (bo002) with
  min(c_s^2) ~ 8.34e-4. Most attempts to lower Planck chi2 further crossed
  the gradient-stability wall. D0 and z_t were fixed in that search, so this
  stage tests whether those coordinates can move the stability boundary.

The exact bo002 ledger is used only for promotion ranking:
  Planck delta = -19.92492142113997
  DESI delta   = -11.902513204929505
  fair P+D     = -0.675202390114676
SN penalties of the frozen sobol107 background:
  Pantheon+ = +1.631268
  Union3    = +1.395445
  DES-Y5    = +3.275618
Thus two-SN closure needs fair P+D < -1.631268.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc,norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern,WhiteKernel,ConstantKernel
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure107_bo002_extended_structural_bo"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0=69.71275747716427
OB=0.022083219194622913; OC=0.12299536722293603
NS=0.9625227132590487; TAU=0.055202901571989066
Q=2.943463801247999; AS=math.exp(Q+2*TAU)/1e10
AL=0.00792232021316886; BL=0.012283301661722363
LAM=18.40625; DN=.5; TAUA=.25; TAUB=1.5

# exact bo002 center
X0=np.array([0.024290704212870225,3.6096407580714724,0.3631213944496205,
             0.04602317962943721,0.34231919445927034,16.173189924377947])
BO002_LITE=1021.7161397407683
BO002_EXACT_PD_FAIR=0.0
SN_PP=0.45564042031764984; SN_U3=0.4392293671844527; SN_D5=1.022321954369545

# [AF, ZC, width, D_floor, D0, z_t]
LOW=np.array([0.0185,3.20,0.285,0.0360,0.260,13.5])
HIGH=np.array([0.0315,4.35,0.470,0.0600,0.430,19.0])

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
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
    a=np.loadtxt(path); e=a[:,0].astype(int); n=int(e.max())+1; cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[e]=a[:,1]*cv; ee[e]=a[:,2]*cv; te[e]=a[:,3]*cv
    L=e.astype(float); pp[e]=a[:,5]*L*(L+1)
    return np.column_stack([e,tt[e],te[e],ee[e]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    ha,d=load_cls(path)
    def f(A):
        x=float(high.chi_squared(ha,A_planck=A))
        x+=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        x+=float(-2*lens.log_likelihood(d,**lp))
        x+=((A-1)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def ini(root,x):
    AF,ZC,W,DF,D0,ZT=map(float,x)
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
expansion_smg={ox:.17g},{LAM},{ZT:.17g},{DN},{AL},{TAUA},{BL},{TAUB}
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

def evaluate(tag,x):
    x=np.asarray(x,float); root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,x))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    r={"id":tag,"AF":x[0],"ZC":x[1],"width":x[2],"D_floor":x[3],"D0":x[4],"z_t":x[5],
       "status":"FAIL","returncode":cp.returncode}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        r.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                 max_cs2=float(bg["c_s^2"].max()),
                 max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1 and r["max_abs_noslip"]<1e-4
        r["stable_subluminal"]=bool(stable)
        if stable:
            pc,Ap=pscore(clp); r.update(status="OK",chi2_planck=pc,A_planck=Ap)
            dl=pc-BO002_LITE
            pred=BO002_EXACT_PD_FAIR+dl
            r.update(delta_lite_vs_bo002=dl,pred_screen_shift=pred,
                     pred_joint_pp=pred+SN_PP,pred_joint_u3=pred+SN_U3,pred_joint_d5=pred+SN_D5)
            js=[r["pred_joint_pp"],r["pred_joint_u3"],r["pred_joint_d5"]]
            r["pred_second_joint"]=sorted(js)[1]
            r["pred_n_sn_closed"]=sum(v<0 for v in js)
            r["goal_score"]=max(pred,r["pred_second_joint"])
    if r["status"]!="OK": r["error"]=cp.stdout[-650:].replace("\n"," | ")
    print("SN107EXT_POINT",json.dumps(r,sort_keys=True),flush=True)
    # Keep only ini and table; spectra/backgrounds are not needed after scoring.
    for p in OUT.glob(tag+"_00_*"):
        try:p.unlink()
        except:pass
    return r

rows=[evaluate("center",X0)]
# Deliberate one-dimensional wall-moving anchors around bo002.
steps=[
 [0,0,0,0, +0.020,0],[0,0,0,0,-0.020,0],
 [0,0,0,0,0,+0.8],[0,0,0,0,0,-0.8],
 [+0.0007,0,0,0,0,0],[-0.0007,0,0,0,0,0],
 [0,+0.12,0,0,0,0],[0,-0.12,0,0,0,0],
 [0,0,+0.018,0,0,0],[0,0,-0.018,0,0,0],
 [0,0,0,+0.002,0,0],[0,0,0,-0.002,0,0]
]
for i,s in enumerate(steps):
    rows.append(evaluate(f"axis{i:02d}",np.clip(X0+np.array(s),LOW,HIGH)))

# 32-point Sobol initialization.
sam=qmc.Sobol(d=6,scramble=True,seed=107602)
for i,x in enumerate(qmc.scale(sam.random_base2(m=5),LOW,HIGH)):
    rows.append(evaluate(f"seed{i:03d}",x))

def fit_gp(rows):
    ok=[r for r in rows if r.get("status")=="OK"]
    X=np.array([[r["AF"],r["ZC"],r["width"],r["D_floor"],r["D0"],r["z_t"]] for r in ok],float)
    y=np.array([r["chi2_planck"] for r in ok],float)
    Xn=(X-LOW)/(HIGH-LOW)
    ker=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=np.ones(6)*0.25,
        length_scale_bounds=(0.025,4.0),nu=2.5)+WhiteKernel(1e-6,(1e-9,1e-2))
    gp=GaussianProcessRegressor(kernel=ker,normalize_y=True,n_restarts_optimizer=4,random_state=107603)
    gp.fit(Xn,y)
    return gp,float(y.min())

rng=np.random.default_rng(107604)
for it in range(24):
    gp,ybest=fit_gp(rows)
    pool=rng.uniform(LOW,HIGH,size=(10000,6))
    ok=[r for r in rows if r.get("status")=="OK"]
    rb=min(ok,key=lambda r:r["chi2_planck"])
    xb=np.array([rb["AF"],rb["ZC"],rb["width"],rb["D_floor"],rb["D0"],rb["z_t"]])
    local=xb+rng.normal(size=(5000,6))*0.05*(HIGH-LOW)
    pool=np.vstack([pool,np.clip(local,LOW,HIGH)])
    xn=(pool-LOW)/(HIGH-LOW)
    mu,std=gp.predict(xn,return_std=True)
    imp=ybest-mu-0.001
    z=np.divide(imp,std,out=np.zeros_like(imp),where=std>1e-12)
    ei=imp*norm.cdf(z)+std*norm.pdf(z)
    Xold=np.array([[r["AF"],r["ZC"],r["width"],r["D_floor"],r["D0"],r["z_t"]] for r in rows])
    Xoldn=(Xold-LOW)/(HIGH-LOW)
    chosen=None
    for j in np.argsort(ei)[::-1]:
        if np.min(np.sqrt(np.sum((Xoldn-xn[j])**2,axis=1)))>0.010:
            chosen=pool[j]; break
    if chosen is None: chosen=pool[int(np.argmax(ei))]
    rows.append(evaluate(f"bo{it:03d}",chosen))

df=pd.DataFrame(rows)
df.to_csv(OUT/"snclosure107_bo002_extended_structural_bo.csv",index=False)
ok=df[df.status=="OK"].sort_values(["goal_score","chi2_planck"])
best=ok.head(20).to_dict("records")
summary={"background":{"H0":H0,"A":AL,"B":BL},"screen_zero_point_note":"P+D zero point intentionally unset; ranking is by Planck-lite improvement on fixed late-local006 background",
         "planck_lite_center":BO002_LITE,
         "n_total":int(len(df)),"n_ok":int(len(ok)),"best":best,
         "note":"Ranking only. Exact full-Plik + raw DESI promotion is mandatory; no raw-DESI proxy is used for the final result."}
(OUT/"snclosure107_bo002_extended_structural_bo_summary.json").write_text(json.dumps(summary,indent=2))
print("SN107EXT_BEST",json.dumps(summary,sort_keys=True),flush=True)
