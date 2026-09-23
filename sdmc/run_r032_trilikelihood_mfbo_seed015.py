#!/usr/bin/env python3
"""
SDMC tri-likelihood multi-fidelity Bayesian search around the current
SN-friendly sobol036 + seed015 basin.

Cheap layer evaluated at every point:
  * exact Pantheon+ covariance likelihood
  * exact Union3 covariance likelihood
  * exact DES-Y5 covariance likelihood
  * DESI DR1 Gaussian BAO
  * Planck 2018 native-lite + low-l + lensing
  * exact hi_class stability diagnostics

The cheap P+D score is corrected with a Gaussian-process discrepancy model
trained on already measured full-Plik + raw-DESI points. A second GP then
optimizes the predicted exact stopping criterion:
  max(P+D_fair, P+D_fair + second-smallest SN penalty).

Top points must still be promoted through full-Plik + raw DESI full shape.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar
from scipy.stats import qmc, norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/trilikelihood_mfbo_seed015"); OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"tmp"; TMP.mkdir(exist_ok=True)
SNROOT=Path("sn_data"); BAOROOT=Path("bao_data")
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Ordinary sector: retain the small coupled sobol036 displacement.
TORD=0.014789129979908472
OB0=0.022083219194622913; OB1=0.022215084896497428
OC0=0.12299536722293603;  OC1=0.12228109714277088
NS0=0.9625227132590487;  NS1=0.9647922897702084
AS0=2.119721074504234e-9; AS1=2.120753981008663e-9
OB=OB0+TORD*(OB1-OB0); OC=OC0+TORD*(OC1-OC0)
NS=NS0+TORD*(NS1-NS0); AS=AS0+TORD*(AS1-AS0)
TAU=0.055202901571989066

# Seed015 structural center. AF and width are frozen in this first joint BO;
# the two strongest edge directions ZC and D_floor remain open.
AF=0.0235834146537818
WIDTH=0.3689727572351694
D0=0.34231919445927034
LAM=18.40625; ZT=16.173189924377947; DNT=.5; TAUA=.25; TAUB=1.5

# Reference exact/cheap values at sobol036 + seed015.
A_REF=0.011249459300190211
B_REF=0.01280436261370778
H_REF=69.80878993044607
ZC_REF=3.4274108000472188
DF_REF=0.050407338682562114
PLITE_REF=1021.2202619980889
BAO_REF=13.204567845148436
EXACT_PD_REF=0.05337920225869297

LOCAL_BAO_CHI=13.400290718544086
LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}
RD_FIXED=145.8653484988842
ZSTAR=1089.9

# Five-dimensional local search box: late A,B,H0 plus the two strongest
# structural edge directions from seed015.
LOW=np.array([0.0065,0.0115,69.68,3.38,0.045])
HIGH=np.array([0.0120,0.0145,69.98,3.52,0.056])

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov mismatch {path}: {len(a)} != {n*n}")
    return a.reshape(n,n)

def setup_cov(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def chi_profile(pack,r):
    cf,one,den=pack
    y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

# SN products
pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm))
m=ppz>0.01; ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]
pppack=setup_cov(Cpp)

p=SNROOT/"Union3/lcparam_full.txt"
cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float)
unpack=setup_cov(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv")
cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
depack=setup_cov(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

# DESI BAO
bao_rows=[]
for ln in (BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
    if not ln.strip() or ln.startswith("#"): continue
    z,v,q=ln.split(); bao_rows.append((float(z),float(v),q))
bao_obs=np.array([x[1] for x in bao_rows])
bao_inv=np.linalg.inv(np.loadtxt(BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_cov.txt"))

# Planck native-lite
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def interp(bg,col,z):
    zz=bg.z.to_numpy(); yy=bg[col].to_numpy(); o=np.argsort(zz)
    return float(np.interp(z,zz[o],yy[o]))

def mu_bg(bg,z,zh):
    zz=bg.z.to_numpy(); dm=bg["comov. dist."].to_numpy(); o=np.argsort(zz)
    DM=np.interp(z,zz[o],dm[o]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def sn_scores(bg):
    return {
      "pantheonplus":chi_profile(pppack,ppm-mu_bg(bg,ppz,ppzh)),
      "union3":chi_profile(unpack,unm-mu_bg(bg,unz,unz)),
      "desy5":chi_profile(depack,dem-mu_bg(bg,dez,dezh)),
    }

def bao_score(bg):
    pred=[]
    for z,_,q in bao_rows:
        DM=interp(bg,"comov. dist.",z); H=interp(bg,"H [1/Mpc]",z)
        DH=1/H; DV=(z*DM*DM*DH)**(1/3)
        pred.append({"DM_over_rs":DM/RD_FIXED,"DH_over_rs":DH/RD_FIXED,"DV_over_rs":DV/RD_FIXED}[q])
    d=np.array(pred)-bao_obs
    return float(d@bao_inv@d)

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

def ini_text(root,A,B,H0,ZC,DF):
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
parameters_smg={AF:.17g},{ZC:.17g},{WIDTH:.17g},{D0:.17g},1.0,{DF:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},{DNT},{A:.17g},{TAUA},{B:.17g},{TAUB}
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

# Exact discrepancy anchors. residual = exact fair P+D - local cheap proxy.
ANCHORS=[
 # A,B,H0,ZC,DF, Planck-lite, BAO, exact fair P+D
 (A_REF,B_REF,H_REF,ZC_REF,DF_REF,PLITE_REF,BAO_REF,0.05337920225869297),
 (A_REF,B_REF,H_REF,3.6096407580714724,0.04602317962943721,1021.6070095999479,BAO_REF,0.37387272101582525),
 (A_REF,B_REF,H_REF,3.4096407580714723,0.04602317962943721,1021.4456072757796,BAO_REF,0.22398869626735163),
 (0.00792232021316886,0.012283301661722363,69.71275747716427,3.6096407580714724,0.04602317962943721,1021.7161397407683,13.251270395581628,0.8931202224049706),
 (0.0013601897284388,0.0147392870020121,70.09917331933976,3.6096407580714724,0.04602317962943721,1021.210755132605,13.433004816867827,-0.675200032792783),
]
def local_proxy(plite,bao):
    return EXACT_PD_REF+(plite-PLITE_REF)+0.25*(bao-BAO_REF)

Xa=[]; ya=[]
for A,B,H,ZC,DF,pl,ba,ex in ANCHORS:
    Xa.append([A,B,H,ZC,DF]); ya.append(ex-local_proxy(pl,ba))
Xa=np.asarray(Xa,float); ya=np.asarray(ya,float)
Xan=(Xa-LOW)/(HIGH-LOW)
corr_gp=GaussianProcessRegressor(
 kernel=ConstantKernel(1.0,(1e-3,10))*Matern(length_scale=np.ones(5)*0.7,length_scale_bounds=(0.08,5.0),nu=2.5)+WhiteKernel(0.02,(1e-5,0.5)),
 normalize_y=True,n_restarts_optimizer=2,random_state=1501)
corr_gp.fit(Xan,ya)

def evaluate(tag,x):
    A,B,H0,ZC,DF=map(float,x)
    root=str(TMP/(tag+"_")); ip=TMP/(tag+".ini")
    ip.write_text(ini_text(root,A,B,H0,ZC,DF))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    r={"id":tag,"A":A,"B":B,"H0":H0,"ZC":ZC,"D_floor":DF,"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        r.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                 max_cs2=float(bg["c_s^2"].max()),
                 max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1
        r["stable_subluminal"]=bool(stable)
        if stable:
            sn=sn_scores(bg)
            for k,v in sn.items():
                r["sn_"+k]=v; r["sn_delta_"+k]=v-LOCAL_SN[k]
            r["bao_chi2"]=bao_score(bg)
            pc,Ap=pscore(clp); r["chi2_planck"]=pc; r["A_planck"]=Ap
            r["pd_proxy_local"]=local_proxy(pc,r["bao_chi2"])
            xn=((np.array([[A,B,H0,ZC,DF]])-LOW)/(HIGH-LOW))
            cm,cs=corr_gp.predict(xn,return_std=True)
            r["exact_corr_mu"]=float(cm[0]); r["exact_corr_std"]=float(cs[0])
            r["pd_exact_pred"]=r["pd_proxy_local"]+r["exact_corr_mu"]
            ss=sorted([r["sn_delta_pantheonplus"],r["sn_delta_union3"],r["sn_delta_desy5"]])
            r["second_sn_penalty"]=float(ss[1])
            r["joint_pp_pred"]=r["pd_exact_pred"]+r["sn_delta_pantheonplus"]
            r["joint_u3_pred"]=r["pd_exact_pred"]+r["sn_delta_union3"]
            r["joint_d5_pred"]=r["pd_exact_pred"]+r["sn_delta_desy5"]
            r["second_joint_pred"]=sorted([r["joint_pp_pred"],r["joint_u3_pred"],r["joint_d5_pred"]])[1]
            r["goal_pred"]=max(r["pd_exact_pred"],r["second_joint_pred"])
            r["n_sn_closed_pred"]=sum(v<0 for v in [r["joint_pp_pred"],r["joint_u3_pred"],r["joint_d5_pred"]])
            r["status"]="OK"
    if r["status"]!="OK": r["error"]=cp.stdout[-600:].replace("\n"," | ")
    print("TRILIK_POINT",json.dumps(r,sort_keys=True),flush=True)
    return r

# Initial design includes the exact center and key anchors plus Sobol seeds.
pts=[
 ("seed015_center",np.array([A_REF,B_REF,H_REF,ZC_REF,DF_REF])),
 ("sobol036_center",np.array([A_REF,B_REF,H_REF,3.6096407580714724,0.04602317962943721])),
 ("sobol036_zcm",np.array([A_REF,B_REF,H_REF,3.4096407580714723,0.04602317962943721])),
 ("late006",np.array([0.00792232021316886,0.012283301661722363,69.71275747716427,3.6096407580714724,0.04602317962943721])),
]
sob=qmc.Sobol(d=5,scramble=True,seed=1502)
for i,x in enumerate(qmc.scale(sob.random_base2(m=4),LOW,HIGH)):
    pts.append((f"seed{i:03d}",np.asarray(x,float)))
rows=[evaluate(tag,x) for tag,x in pts]

def fit_goal(rows):
    ok=[r for r in rows if r.get("status")=="OK"]
    X=np.array([[r["A"],r["B"],r["H0"],r["ZC"],r["D_floor"]] for r in ok],float)
    y=np.array([r["goal_pred"] for r in ok],float)
    Xn=(X-LOW)/(HIGH-LOW)
    gp=GaussianProcessRegressor(
      kernel=ConstantKernel(1.0,(1e-3,10))*Matern(length_scale=np.ones(5)*0.35,length_scale_bounds=(0.04,3.0),nu=2.5)+WhiteKernel(1e-5,(1e-8,0.1)),
      normalize_y=True,n_restarts_optimizer=2,random_state=1503)
    gp.fit(Xn,y)
    return gp,float(y.min())

rng=np.random.default_rng(1504)
for it in range(20):
    gp,ybest=fit_goal(rows)
    pool=rng.uniform(LOW,HIGH,size=(9000,5))
    ok=[r for r in rows if r.get("status")=="OK"]
    rb=min(ok,key=lambda r:r["goal_pred"])
    xb=np.array([rb["A"],rb["B"],rb["H0"],rb["ZC"],rb["D_floor"]])
    local=xb+rng.normal(size=(3000,5))*0.055*(HIGH-LOW)
    pool=np.vstack([pool,np.clip(local,LOW,HIGH)])
    xn=(pool-LOW)/(HIGH-LOW)
    mu,std=gp.predict(xn,return_std=True)
    imp=ybest-mu-0.002
    z=np.divide(imp,std,out=np.zeros_like(imp),where=std>1e-12)
    ei=imp*norm.cdf(z)+std*norm.pdf(z)
    # Encourage candidates in regions where exact P+D discrepancy remains uncertain.
    _,corrstd=corr_gp.predict(xn,return_std=True)
    acq=ei+0.04*corrstd
    Xold=np.array([[r["A"],r["B"],r["H0"],r["ZC"],r["D_floor"]] for r in rows],float)
    Xoldn=(Xold-LOW)/(HIGH-LOW)
    chosen=None
    for j in np.argsort(acq)[::-1]:
        if np.min(np.sqrt(np.sum((Xoldn-xn[j])**2,axis=1)))>0.010:
            chosen=pool[j]; break
    rows.append(evaluate(f"bo{it:03d}",chosen))

df=pd.DataFrame(rows)
df.to_csv(OUT/"trilikelihood_mfbo.csv",index=False)
ok=df[df.status.eq("OK")].sort_values(["goal_pred","second_joint_pred","pd_exact_pred"])
best=ok.head(12).to_dict("records")
summary={
 "search":"multi-fidelity Bayesian Planck+DESI+SN",
 "cheap_layer":["Planck native-lite","DESI DR1 Gaussian BAO","Pantheon+","Union3","DES-Y5","stability"],
 "exact_correction_anchors":len(ANCHORS),
 "n_ok":int(len(ok)),
 "best":best,
 "promotion_rule":"promote leading stable points to exact full-Plik + raw DESI; stop only if P+D fair<0 and at least two SN joint gaps<0"
}
(OUT/"trilikelihood_mfbo_summary.json").write_text(json.dumps(summary,indent=2))
print("TRILIK_BEST",json.dumps(summary,sort_keys=True),flush=True)
