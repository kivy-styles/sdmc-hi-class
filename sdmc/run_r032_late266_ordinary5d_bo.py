#!/usr/bin/env python3
"""
Local 5D ordinary-sector Bayesian optimization for Planck recovery around exact late266.

Frozen:
  H0 and late background = exact late266
  SDMC structure = exactbo000

Varied:
  omega_b, omega_cdm, n_s, Q, tau
  where Q = ln(1e10 A_s) - 2 tau.

Screen objective:
  exact BO002 fair Planck+DESI baseline
+ Planck-lite change relative to BO002 center
+ candidate SN delta relative to local021.

The optimizer minimizes the second-best of the three joint SN totals, i.e. the
quantity that must cross zero to satisfy "Planck+DESI + at least two SN".

Raw DESI is NOT approximated as independent data here; top candidates must be
replayed with the official raw DESI full-shape likelihood and full Plik.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap, warnings
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler
from scipy.stats import norm

warnings.filterwarnings("ignore")
OUT=Path("output/late266_ordinary5d")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255
CAL_SIGMA=0.0025
OR=4.17998772e-5

# Frozen exact late266 geometry / exactbo000 structure.
H0=69.45160505326464
AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
DF=0.050275813910968366
LAM=18.40625
ZT=16.173189924377947
ALATE=0.013036667098291217
BLATE=0.010544837576337158
TAUA=0.25
TAUB=1.5
DNT=0.5

# Current ordinary coordinates.
OB0=0.02208511574370519
OC0=0.12298509428428875
NS0=0.962555355281866
TAU0=0.055202901571989066
AS0=2.1197359302086e-9
Q0=math.log(1e10*AS0)-2*TAU0
X0=np.array([OB0,OC0,NS0,Q0,TAU0],float)
NAMES=["omega_b","omega_cdm","n_s","Q","tau_reio"]

# Local search box: deliberately narrower than the old global 6D screen.
LOW=np.array([OB0-0.00040,OC0-0.0030,NS0-0.0100,Q0-0.0180,max(0.038,TAU0-0.0120)])
HIGH=np.array([OB0+0.00040,OC0+0.0030,NS0+0.0100,Q0+0.0180,min(0.075,TAU0+0.0120)])

# Exact baseline already established on bo002.
LOCAL021_PD=0.0
BO002_PLANCK=0.0
BO002_DESI=0.0
BO002_PD_FAIR=BO002_PLANCK+BO002_DESI-LOCAL021_PD

# local021 parameters and exact SN reference.
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
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
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
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

# ---------- SN data ----------
DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"covariance mismatch {path}")
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))
def pchi(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm)); m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=setup(Cpp)
p=DATA/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unp=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(unm)))
de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float); dep=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

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
LOCAL_SN=dict(pantheonplus=1406.2147033223882,union3=28.783140002196888,desy5=1650.6353947147727)
print("LATE266_ORD5_LOCAL021_SN",json.dumps(LOCAL_SN,sort_keys=True),flush=True)

def sdmc_ini(x,root):
    ob,oc,ns,Q,tau=map(float,x)
    As=math.exp(Q+2*tau)/1e10
    h=H0/100.; Ox=1.-(ob+oc+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={ob:.17g}
omega_cdm={oc:.17g}
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
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={Ox:.17g},{LAM},{ZT},{DNT},{ALATE},{TAUA},{BLATE},{TAUB}
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
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(sdmc_ini(x,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    ob,oc,ns,Q,tau=map(float,x); As=math.exp(Q+2*tau)/1e10
    rec=dict(id=tag,omega_b=ob,omega_cdm=oc,n_s=ns,Q=Q,tau_reio=tau,A_s=As,
             status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        rec["min_D"]=float(bg["kin (D)"].min()); rec["min_cs2"]=float(bg["c_s^2"].min())
        rec["max_cs2"]=float(bg["c_s^2"].max())
        rec["max_abs_noslip"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            ps=planck_score(clp); sn=sn_scores(bg)
            rec.update(ps); rec.update({f"sn_{k}":v for k,v in sn.items()})
            rec["sn_delta_pantheonplus"]=sn["pantheonplus"]-LOCAL_SN["pantheonplus"]
            rec["sn_delta_union3"]=sn["union3"]-LOCAL_SN["union3"]
            rec["sn_delta_desy5"]=sn["desy5"]-LOCAL_SN["desy5"]
            rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    for p in OUT.glob(tag+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return rec

center=evaluate("center",X0)
if center["status"]!="OK": raise RuntimeError(center)
CPL=center["chi2_planck"]
print("LATE266_ORD5_CENTER",json.dumps(center,sort_keys=True),flush=True)

def finish(rec):
    if rec["status"]!="OK":
        rec["delta_planck_lite"]=99999.
        rec["second_sn_delta"]=99999.
        rec["goal_score"]=99999.
        return rec
    dp=rec["chi2_planck"]-CPL
    rec["delta_planck_lite"]=float(dp)
    vals=sorted([rec["sn_delta_pantheonplus"],rec["sn_delta_union3"],rec["sn_delta_desy5"]])
    rec["second_sn_delta"]=float(vals[1])
    rec["goal_score"]=float(dp + 800.0*max(0.0,rec["second_sn_delta"]))
    return rec

center=finish(center)
rows=[center]

# Initial 64 Sobol points.
sob=qmc.Sobol(d=5,scramble=True,seed=266023)
U=sob.random_base2(m=6)
X=qmc.scale(U,LOW,HIGH)
for i,x in enumerate(X):
    r=finish(evaluate(f"seed{i:03d}",x)); rows.append(r)
    print("LATE266_ORD5_POINT",json.dumps(r,sort_keys=True),flush=True)

# Bayesian refinement: 22 sequential EI steps.
for it in range(22):
    ok=[r for r in rows if r["status"]=="OK" and np.isfinite(r["goal_score"])]
    Xa=np.array([[r[k] for k in NAMES] for r in ok],float)
    ya=np.array([r["goal_score"] for r in ok],float)
    sx=StandardScaler().fit(Xa); Xs=sx.transform(Xa)
    kernel=ConstantKernel(1.0,(1e-2,1e2))*Matern(length_scale=np.ones(5),nu=2.5)+WhiteKernel(1e-6,(1e-9,1e-2))
    gp=GaussianProcessRegressor(kernel=kernel,alpha=1e-8,normalize_y=True,n_restarts_optimizer=2,random_state=26600+it)
    gp.fit(Xs,ya)
    rng=np.random.default_rng(266200+it)
    cand=rng.uniform(LOW,HIGH,size=(5000,5))
    cs=sx.transform(cand); mu,sd=gp.predict(cs,return_std=True)
    best=float(np.min(ya)); imp=best-mu-0.002
    z=np.divide(imp,sd,out=np.zeros_like(imp),where=sd>1e-12)
    ei=imp*norm.cdf(z)+sd*norm.pdf(z)
    # discourage duplicates
    for j,c in enumerate(cand):
        if np.min(np.linalg.norm((Xa-c)/(HIGH-LOW),axis=1))<0.015: ei[j]=-1
    x=cand[int(np.argmax(ei))]
    r=finish(evaluate(f"bo{it:03d}",x)); rows.append(r)
    print("LATE266_ORD5_BO_POINT",json.dumps(r,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"late266_ordinary5d_bo.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].sort_values("goal_score")
top=ok.head(12)
top.to_csv(OUT/"late266_ordinary5d_top12.csv",index=False)
summary=dict(
  target="late266 Planck recovery with >=2 exact SN deltas <= local021",
  center_planck_lite=CPL,
  bounds={k:[float(LOW[i]),float(HIGH[i])] for i,k in enumerate(NAMES)},
  best=top.iloc[0].to_dict(),
  top12=top.to_dict("records"),
  note="SN deltas are relative to LCDM local021 exact covariance values; top points require exact raw-DESI + full-Plik replay. Exact late266 center has Planck gap +3.2273 and DESI margin -1.2418 versus local021."
)
(OUT/"late266_ordinary5d_summary.json").write_text(json.dumps(summary,indent=2))
print("LATE266_ORD5_BEST",json.dumps(top.to_dict("records"),sort_keys=True),flush=True)
