#!/usr/bin/env python3
"""
Transition-width Bayesian refinement around exact sobol076 + BO002.

This stage reopens three expansion-shape coordinates that were held fixed
through the historical lambda_e-z_t scans:
    DeltaN_bg, tau_A, tau_B.

Frozen:
  ordinary cosmology = successful sobol107 solution
  low-z amplitudes A_late, B_late and H0 = sobol107
  structure = bo002
  lambda_e and z_t = current exact solution

Screen:
  * exact hi_class stability: D>0, c_s^2>0, c_s^2<=1, no-slip residual
  * Planck native-lite + lowT + lowE + lensing
  * Pantheon+, Union3, DES-Y5 full covariance shape likelihoods
  * exact BO002 fair P+D baseline is shifted only by the Planck-lite response

Raw DESI is NOT replaced by a proxy conclusion. Any finalist must be promoted
to the official raw DESI full-shape likelihood and full Plik.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap, warnings
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc, norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

warnings.filterwarnings("ignore")
OUT=Path("output/snclosure107_bo002_sobol076_transitionwidth"); OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"tmp"; TMP.mkdir(exist_ok=True)
DATA=Path("sn_data")
TCMB=2.7255; CAL_SIGMA=0.0025; OR=4.17998772e-5

# Frozen exact BO002 solution.
H0=70.00876109628007
OB=0.022083219194622913
OC=0.12299536722293603
AS=2.119721074504234e-9
NS=0.9625227132590487
TAU=0.055202901571989066
AF=0.024290704212870225
ZC=3.6096407580714724
WIDTH=0.3631213944496205
D0=0.34231919445927034
DF=0.04602317962943721
LAM=18.40625
ZT=16.173189924377947
ALATE=0.01686786537989974
BLATE=0.014543195590376853

# Historically fixed transition-shape center.
X0=np.array([0.5,0.25,1.5],float)  # DeltaN_bg, tau_A, tau_B
NAMES=["DeltaN_bg","tau_A","tau_B"]
LOW=np.array([0.32,0.14,0.90],float)
HIGH=np.array([0.72,0.42,2.20],float)

# Exact BO002 baseline.
BO002_PD_FAIR=-0.659891015199548

# local021 SN reference cosmology.
LOCAL021=dict(H0=68.56858744695782,omega_b=0.022406369378007947,
             omega_cdm=0.11824151052483357,A_s=2.05109266435559e-9,
             n_s=0.9649593164240942,tau=0.044985382026527077)

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
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*cv; ee[ell]=a[:,2]*cv; te[ell]=a[:,3]*cv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def planck_score(path):
    harr,dls=load_cls(path)
    def pieces(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

# --- SN covariance products ---
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov mismatch {path}")
    return a.reshape(n,n)

def setup_cov(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def pchi(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm)); m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=setup_cov(Cpp)

p=DATA/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float)
unp=setup_cov(read_cov(DATA/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup_cov(read_cov(DATA/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

def mu_bg(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def sn_scores(bg):
    return dict(
      pantheonplus=pchi(ppp,ppm-mu_bg(bg,ppz,ppzh)),
      union3=pchi(unp,unm-mu_bg(bg,unz,unz)),
      desy5=pchi(dep,dem-mu_bg(bg,dez,dezh))
    )

def lcdm_ini(root):
    x=LOCAL021
    return textwrap.dedent(f"""\
H0={x['H0']}
omega_b={x['omega_b']}
omega_cdm={x['omega_cdm']}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={x['A_s']}
n_s={x['n_s']}
tau_reio={x['tau']}
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

lroot=str(OUT/"local021_"); lip=OUT/"local021.ini"; lip.write_text(lcdm_ini(lroot))
cp=subprocess.run(["./class",str(lip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-3000:])
LOCAL_SN=sn_scores(table(lroot+"00_background.dat"))
print("S076TW_LOCAL021_SN",json.dumps(LOCAL_SN,sort_keys=True),flush=True)

def sdmc_ini(x,root):
    DN,TAUA,TAUB=map(float,x)
    h=H0/100.; Ox=1.-(OB+OC+OR)/(h*h)
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
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={Ox:.17g},{LAM},{ZT},{DN:.17g},{ALATE},{TAUA:.17g},{BLATE},{TAUB:.17g}
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
    x=np.asarray(x,float); root=str(TMP/(tag+"_")); ip=TMP/(tag+".ini")
    ip.write_text(sdmc_ini(x,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    DN,TAUA,TAUB=map(float,x)
    rec=dict(id=tag,DeltaN_bg=DN,tau_A=TAUA,tau_B=TAUB,status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        rec["min_D"]=float(bg["kin (D)"].min()); rec["min_cs2"]=float(bg["c_s^2"].min())
        rec["max_cs2"]=float(bg["c_s^2"].max())
        rec["max_abs_noslip"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1 and rec["max_abs_noslip"]<1e-4)
        rec["stable_subluminal"]=bool(stable)
        if stable:
            ps=planck_score(clp); sn=sn_scores(bg); rec.update(ps)
            for k,v in sn.items():
                rec[f"sn_{k}"]=v; rec[f"sn_delta_{k}"]=v-LOCAL_SN[k]
            rec["status"]="OK"
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    for q in TMP.glob(tag+"_*"):
        try:q.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return rec

center=evaluate("center",X0)
if center["status"]!="OK": raise RuntimeError(center)
CENTER_LITE=center["chi2_planck"]
print("S076TW_CENTER",json.dumps(center,sort_keys=True),flush=True)

def finish(rec):
    if rec["status"]!="OK":
        rec.update(delta_planck_lite=999.,pd_fair_proxy=999.,joint_pp_proxy=999.,
                   joint_u3_proxy=999.,joint_d5_proxy=999.,second_joint_proxy=999.,
                   n_sn_closed_proxy=0,goal_score=999.)
        return rec
    dp=rec["chi2_planck"]-CENTER_LITE
    pd=BO002_PD_FAIR+dp
    js=[pd+rec["sn_delta_pantheonplus"],pd+rec["sn_delta_union3"],pd+rec["sn_delta_desy5"]]
    rec.update(delta_planck_lite=dp,pd_fair_proxy=pd,
               joint_pp_proxy=js[0],joint_u3_proxy=js[1],joint_d5_proxy=js[2],
               second_joint_proxy=sorted(js)[1],n_sn_closed_proxy=sum(v<0 for v in js))
    rec["goal_score"]=max(pd,rec["second_joint_proxy"])
    return rec

rows=[finish(center)]
# One-axis anchors around the historical center.
anchors=[
 ("dn_lo",[0.42,0.25,1.5]),("dn_hi",[0.60,0.25,1.5]),
 ("ta_lo",[0.5,0.19,1.5]),("ta_hi",[0.5,0.33,1.5]),
 ("tb_lo",[0.5,0.25,1.20]),("tb_hi",[0.5,0.25,1.85])
]
for tag,x in anchors:
    r=finish(evaluate(tag,x)); rows.append(r); print("S076TW_POINT",json.dumps(r,sort_keys=True),flush=True)

sam=qmc.Sobol(d=3,scramble=True,seed=107803)
for i,x in enumerate(qmc.scale(sam.random_base2(m=5),LOW,HIGH)):
    r=finish(evaluate(f"seed{i:03d}",x)); rows.append(r)
    print("S076TW_POINT",json.dumps(r,sort_keys=True),flush=True)

def fit_gp(rows):
    ok=[r for r in rows if r.get("status")=="OK"]
    X=np.array([[r[n] for n in NAMES] for r in ok],float)
    y=np.array([r["goal_score"] for r in ok],float)
    Xn=(X-LOW)/(HIGH-LOW)
    ker=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=np.ones(3)*0.28,
        length_scale_bounds=(0.04,3.0),nu=2.5)+WhiteKernel(1e-6,(1e-9,1e-2))
    gp=GaussianProcessRegressor(kernel=ker,normalize_y=True,n_restarts_optimizer=3,random_state=107804)
    gp.fit(Xn,y)
    return gp,float(y.min())

rng=np.random.default_rng(107805)
for it in range(14):
    gp,ybest=fit_gp(rows)
    pool=rng.uniform(LOW,HIGH,size=(7000,3))
    ok=[r for r in rows if r.get("status")=="OK"]
    rb=min(ok,key=lambda r:r["goal_score"])
    xb=np.array([rb[n] for n in NAMES])
    local=xb+rng.normal(size=(2500,3))*0.055*(HIGH-LOW)
    pool=np.vstack([pool,np.clip(local,LOW,HIGH)])
    xn=(pool-LOW)/(HIGH-LOW)
    mu,std=gp.predict(xn,return_std=True)
    imp=ybest-mu-0.002
    z=np.divide(imp,std,out=np.zeros_like(imp),where=std>1e-12)
    ei=imp*norm.cdf(z)+std*norm.pdf(z)
    Xold=np.array([[r[n] for n in NAMES] for r in rows],float)
    Xoldn=(Xold-LOW)/(HIGH-LOW)
    pick=None
    for j in np.argsort(ei)[::-1]:
        if np.min(np.sqrt(np.sum((Xoldn-xn[j])**2,axis=1)))>0.012:
            pick=j; break
    if pick is None: pick=int(np.argmax(ei))
    r=finish(evaluate(f"bo{it:03d}",pool[pick])); rows.append(r)
    print("S076TW_BO_POINT",json.dumps(r,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"snclosure107_bo002_sobol076_transitionwidth.csv",index=False)
ok=df[df.status=="OK"].sort_values(["goal_score","second_joint_proxy","pd_fair_proxy"])
top=ok.head(15)
top.to_csv(OUT/"snclosure107_bo002_sobol076_transitionwidth_top15.csv",index=False)
summary=dict(center=X0.tolist(),bounds=dict(low=LOW.tolist(),high=HIGH.tolist()),
             bo002_exact_pd_fair=BO002_PD_FAIR,n_total=int(len(df)),n_ok=int(len(ok)),
             best=top.to_dict("records"),
             promotion_rule="exact full-Plik + official raw DESI full-shape + all three SN datasets")
(OUT/"snclosure107_bo002_sobol076_transitionwidth_summary.json").write_text(json.dumps(summary,indent=2))
print("S076TW_BEST",json.dumps(summary,sort_keys=True),flush=True)
