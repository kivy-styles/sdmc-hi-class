#!/usr/bin/env python3
"""
Joint local Bayesian optimization around the current SN107+bo002 basin.

Reopens simultaneously:
  late background: (A_late, B_late, H0)
  structure:       (A_F, z_c, width, D_floor)

Every point is required to pass exact hi_class background/stability checks.
Pantheon+, Union3 and DES-Y5 shape likelihoods are evaluated with their full
covariance products. Planck is evaluated with the native-lite likelihood.

Raw DESI is NOT replaced by the proxy. For ranking only, the local raw-DESI
movement is estimated from the measured exact relation across three already
promoted backgrounds:
  sobol223, sobol107 and edge024-Q.
Across those exact points the raw DESI delta is very nearly linear in H0,
with slope d(Delta chi2_DESI)/dH0 = -2.13245454 per km/s/Mpc.
Finalists from this search must be promoted to the real raw-DESI likelihood.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc, norm
from scipy.linalg import cho_factor, cho_solve
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure107_bo002_joint_local"); OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"tmp"; TMP.mkdir(exist_ok=True)
SNROOT=Path("sn_data"); BAOROOT=Path("bao_data")
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Ordinary cosmological sector kept at the successful SN107 solution.
OB=0.022083219194622913
OC=0.12299536722293603
NS=0.9625227132590487
TAU=0.055202901571989066
Q=2.943463801247999
AS=math.exp(Q+2*TAU)/1e10

# Fixed SDMC constants.
D0=0.34231919445927034
LAM=18.40625; ZT=16.173189924377947; DN=.5; TAUA=.25; TAUB=1.5

# Current exact best center: sobol107 late background + bo002 structure.
A0=0.0013601897284388; B0=0.0147392870020121; H00=70.09917331933976
AF0=0.024290704212870225; ZC0=3.6096407580714724
W0=0.3631213944496205; DF0=0.04602317962943721

# Exact anchors already measured for the center.
CENTER_PLANCK_LITE=1021.2107551326056
CENTER_PD_FAIR=-0.6752023860694756
DESI_H0_SLOPE=-2.13245454

LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}
RD_FIXED=145.8653484988842
LOCAL_BAO_CHI=13.400290718544086

# Broad enough to interpolate toward the high-H0/Planck-friendly sobol072
# region while remaining local enough for the empirically calibrated DESI proxy.
LOW=np.array([0.0000,0.0110,69.88, 0.0235,3.50,0.315,0.0435])
HIGH=np.array([0.0160,0.0215,70.66, 0.0300,4.18,0.410,0.0550])

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

# SN products.
pp=pd.read_csv(SNROOT/'PantheonPlus/Pantheon+SH0ES.dat',sep=r'\s+',header=0,engine='python')
ppm=pp['m_b_corr'].to_numpy(float); ppz=pp['zHD'].to_numpy(float); ppzh=pp['zHEL'].to_numpy(float)
Cpp=read_cov(SNROOT/'PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov',len(ppm))
m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]
pppack=setup_cov(Cpp)

p=SNROOT/'Union3/lcparam_full.txt'
cols=p.read_text().splitlines()[0].removeprefix('#').split()
un=pd.read_csv(p,sep=r'\s+',comment='#',names=cols,engine='python')
cm={c.lower():c for c in un.columns}
unm=un[cm['mb']].to_numpy(float); unz=un[cm['zcmb']].to_numpy(float)
unpack=setup_cov(read_cov(SNROOT/'Union3/mag_covmat.txt',len(unm)))

de=pd.read_csv(SNROOT/'DESY5/DES-SN5YR_HD.csv')
cm={c.lower():c for c in de.columns}
dem=de[cm['mu']].to_numpy(float); dez=de[cm['zhd']].to_numpy(float); dezh=de[cm['zhel']].to_numpy(float)
derr=de[cm['muerr_final']].to_numpy(float)
depack=setup_cov(read_cov(SNROOT/'DESY5/covsys_000.txt',len(dem))+np.diag(derr*derr))

# Compressed BAO is retained only as a sanity diagnostic, not as the DESI objective.
bao_rows=[]
for ln in (BAOROOT/'desi_2024_gaussian_bao_ALL_GCcomb_mean.txt').read_text().splitlines():
    if not ln.strip() or ln.startswith('#'): continue
    z,v,q=ln.split(); bao_rows.append((float(z),float(v),q))
bao_obs=np.array([x[1] for x in bao_rows])
bao_cov=np.loadtxt(BAOROOT/'desi_2024_gaussian_bao_ALL_GCcomb_cov.txt')
bao_inv=np.linalg.inv(bao_cov)

def interp(bg,col,z):
    zz=bg.z.to_numpy(); yy=bg[col].to_numpy(); o=np.argsort(zz)
    return float(np.interp(z,zz[o],yy[o]))

def mu_bg(bg,z,zh):
    zz=bg.z.to_numpy(); dm=bg['comov. dist.'].to_numpy(); o=np.argsort(zz)
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
        DM=interp(bg,'comov. dist.',z); H=interp(bg,'H [1/Mpc]',z)
        DH=1/H; DV=(z*DM*DM*DH)**(1/3)
        pred.append({"DM_over_rs":DM/RD_FIXED,"DH_over_rs":DH/RD_FIXED,"DV_over_rs":DV/RD_FIXED}[q])
    d=np.array(pred)-bao_obs
    return float(d@bao_inv@d)

def load_cls(path):
    a=np.loadtxt(path); e=a[:,0].astype(int); n=int(e.max())+1; cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); ppv=np.zeros(n)
    tt[e]=a[:,1]*cv; ee[e]=a[:,2]*cv; te[e]=a[:,3]*cv
    L=e.astype(float); ppv[e]=a[:,5]*L*(L+1)
    return np.column_stack([e,tt[e],te[e],ee[e]]),{"tt":tt,"te":te,"ee":ee,"pp":ppv}

def pscore(path):
    ha,d=load_cls(path)
    def f(A):
        x=float(high.chi_squared(ha,A_planck=A))
        x+=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        x+=float(-2*lens.log_likelihood(d,**lp))
        x+=((A-1)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def ini(root,A,B,H0,AF,ZC,W,DF):
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
expansion_smg={ox:.17g},{LAM},{ZT},{DN},{A:.17g},{TAUA},{B:.17g},{TAUB}
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

def run(tag,x):
    A,B,H0,AF,ZC,W,DF=map(float,x)
    root=str(TMP/(tag+"_")); ip=TMP/(tag+".ini"); ip.write_text(ini(root,A,B,H0,AF,ZC,W,DF))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    r={"id":tag,"A":A,"B":B,"H0":H0,"AF":AF,"ZC":ZC,"width":W,"D_floor":DF,"status":"FAIL"}
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
                r[f"sn_{k}"]=v; r[f"sn_delta_{k}"]=v-LOCAL_SN[k]
            r["bao_chi2"]=bao_score(bg)
            pc,Ap=pscore(clp); r["chi2_planck"]=pc; r["A_planck"]=Ap
            dP=pc-CENTER_PLANCK_LITE
            dD=DESI_H0_SLOPE*(H0-H00)
            r["desi_local_proxy_shift"]=dD
            r["pd_fair_proxy"]=CENTER_PD_FAIR+dP+dD
            js=[r["pd_fair_proxy"]+r["sn_delta_pantheonplus"],
                r["pd_fair_proxy"]+r["sn_delta_union3"],
                r["pd_fair_proxy"]+r["sn_delta_desy5"]]
            r["joint_pp_proxy"],r["joint_u3_proxy"],r["joint_d5_proxy"]=js
            r["second_joint_proxy"]=sorted(js)[1]
            r["n_sn_closed_proxy"]=sum(v<0 for v in js)
            # Preserve P+D closure while requiring at least two SN closures.
            r["goal_score"]=max(r["pd_fair_proxy"],r["second_joint_proxy"])
            # Soft sanity penalty only for grossly bad compressed BAO.
            if r["bao_chi2"]>18:
                r["goal_score"]+=0.15*(r["bao_chi2"]-18)
            r["status"]="OK"
    if r["status"]!="OK": r["error"]=cp.stdout[-650:].replace("\n"," | ")
    print("SN107JOINT_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in TMP.glob(tag+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

# Deterministic seeds include current best and useful neighboring late backgrounds.
seeds=[
 ("center",[A0,B0,H00,AF0,ZC0,W0,DF0]),
 ("107_s009",[A0,B0,H00,0.024221895329654217,4.145109392143786,0.3370092310383916,0.04453643597662449]),
 ("107_s013",[A0,B0,H00,0.02582530555780977,3.8840197445824742,0.3657897564768791,0.043821363696828486]),
 ("107_s006",[A0,B0,H00,0.02559563393890858,3.8428389857523144,0.388143997490406,0.04659694366157055]),
 ("072_edge",[0.0134815868735313,0.0181835683900862,70.60292498096824,0.03200255395658314,4.007449422683567,0.34799921004101636,0.04607192634791136]),
 ("223_s110",[0.018589897081255913,0.013815921754576268,69.85635133907199,0.026252536563202738,3.79610941009596,0.3229390993900597,0.046828035046346486])
]
rows=[run(tag,np.array(x,float)) for tag,x in seeds]

# 64-point joint Sobol design.
sam=qmc.Sobol(d=7,scramble=True,seed=107702)
for i,x in enumerate(qmc.scale(sam.random_base2(m=6),LOW,HIGH)):
    rows.append(run(f"sobol{i:03d}",x))

def fit_gp(rows):
    ok=[r for r in rows if r.get("status")=="OK"]
    X=np.array([[r["A"],r["B"],r["H0"],r["AF"],r["ZC"],r["width"],r["D_floor"]] for r in ok])
    y=np.array([r["goal_score"] for r in ok])
    Xn=(X-LOW)/(HIGH-LOW)
    ker=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=np.ones(7)*0.3,
        length_scale_bounds=(0.04,3.0),nu=2.5)+WhiteKernel(1e-6,(1e-9,1e-2))
    gp=GaussianProcessRegressor(kernel=ker,normalize_y=True,n_restarts_optimizer=2,random_state=107703)
    gp.fit(Xn,y)
    return gp,float(y.min())

rng=np.random.default_rng(107704)
for it in range(18):
    gp,ybest=fit_gp(rows)
    pool=rng.uniform(LOW,HIGH,size=(9000,7))
    ok=[r for r in rows if r.get("status")=="OK"]
    rb=min(ok,key=lambda r:r["goal_score"])
    xb=np.array([rb["A"],rb["B"],rb["H0"],rb["AF"],rb["ZC"],rb["width"],rb["D_floor"]])
    local=xb+rng.normal(size=(3000,7))*0.055*(HIGH-LOW)
    local=np.clip(local,LOW,HIGH)
    pool=np.vstack([pool,local])
    xn=(pool-LOW)/(HIGH-LOW)
    mu,std=gp.predict(xn,return_std=True)
    imp=ybest-mu-0.003
    z=np.divide(imp,std,out=np.zeros_like(imp),where=std>1e-12)
    ei=imp*norm.cdf(z)+std*norm.pdf(z)
    Xold=np.array([[r["A"],r["B"],r["H0"],r["AF"],r["ZC"],r["width"],r["D_floor"]] for r in rows])
    Xoldn=(Xold-LOW)/(HIGH-LOW)
    pick=None
    for j in np.argsort(ei)[::-1]:
        if np.min(np.sqrt(np.sum((Xoldn-xn[j])**2,axis=1)))>0.014:
            pick=j; break
    rows.append(run(f"bo{it:03d}",pool[pick]))

df=pd.DataFrame(rows)
df.to_csv(OUT/"snclosure107_bo002_joint_local.csv",index=False)
ok=df[df.status=="OK"].sort_values(["goal_score","second_joint_proxy","pd_fair_proxy"])
best=ok.head(15).to_dict("records")
summary={
 "center":{"A":A0,"B":B0,"H0":H00,"AF":AF0,"ZC":ZC0,"width":W0,"D_floor":DF0,
           "exact_pd_fair":CENTER_PD_FAIR},
 "desi_proxy":{"slope_per_H0":DESI_H0_SLOPE,
               "basis":"linear fit to exact raw-DESI results at sobol223, sobol107 and edge024-Q; promotion still mandatory"},
 "n_ok":int(len(ok)),"best":best,
 "success_proxy":"pd_fair_proxy < 0 and at least two of the three joint SN proxy gaps < 0",
 "promotion_rule":"top stable points must be rerun with exact full-Plik + raw DESI + all three SN datasets"
}
(OUT/"snclosure107_bo002_joint_local_summary.json").write_text(json.dumps(summary,indent=2))
print("SN107JOINT_BEST",json.dumps(summary,sort_keys=True),flush=True)
