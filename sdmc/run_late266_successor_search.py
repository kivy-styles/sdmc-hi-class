#!/usr/bin/env python3
"""
Late266-successor search.

Searches the manuscript-A late266 basin for a stable candidate that improves the
remaining CMB cost without giving up the exact Pantheon+, Union3 and DES-Y5
advantages. The late background is frozen to late266; the search uses the two
directions that survived the Part-IX derivative audit:
  1) decreasing A_F in the independent-kinetic sector;
  2) a narrow primordial A_s-tau-n_s degeneracy direction.

This is a screening stage. No point is called a new winner until promoted to
full Plik + official DESI DR1 full shape and then replayed through the complete
Part-VIII/IX battery.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/late266_successor"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
TCMB=2.7255; CAL_SIGMA=0.0025; OR=4.17998772e-5

# Accepted-action late266 center from Manuscript A.
H0=69.45160505326464
OB=0.02208511574370519
OC=0.12298509428428875
NS0=0.962555355281866
TAU0=0.055202901571989066
AS0=2.1197359302086e-9
AEFF0=AS0*math.exp(-2*TAU0)

AF0=0.023604633340554453
ZC0=3.4328776987879466
W0=0.36879721635160295
D0=0.34231919445927034
DF0=0.050275813910968366
POWER=1.0

LAM=18.40625
ZT=16.173189924377947
DNT=0.5
ALATE=0.013036667098291217
BLATE=0.010544837576337158
TAUA=0.25; TAUB=1.5

LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}

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
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def planck(path):
    harr,dls=load_cls(path)
    def pieces(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pieces(A)[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-9})
    p=pieces(op.x)
    return dict(A_planck=float(op.x),chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck_lite=p[5])

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov mismatch {path}")
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False); return cf,one,float(one@u)
def pchi(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r@y-(one@y)**2/den)

pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm))
m=ppz>0.01; ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=setup(Cpp)

p=SNROOT/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unp=setup(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

def mu_bg(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def sn(bg):
    return {
      "pantheonplus":pchi(ppp,ppm-mu_bg(bg,ppz,ppzh)),
      "union3":pchi(unp,unm-mu_bg(bg,unz,unz)),
      "desy5":pchi(dep,dem-mu_bg(bg,dez,dezh))
    }

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; P=a[:,1]; x=8*k
    W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
    return float(np.sqrt(np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))))

def ini(root,AF,zc,w,df,ns,tau):
    As=AEFF0*math.exp(2*tau)
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
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
parameters_smg={AF:.17g},{zc:.17g},{w:.17g},{D0:.17g},{POWER:.17g},{df:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM:.17g},{ZT:.17g},{DNT:.17g},{ALATE:.17g},{TAUA:.17g},{BLATE:.17g},{TAUB:.17g}
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

def evaluate(label,AF,zc,w,df,ns,tau):
    As=AEFF0*math.exp(2*tau)
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,AF,zc,w,df,ns,tau))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    r=dict(id=label,A_F=float(AF),z_c=float(zc),width=float(w),D_floor=float(df),
           n_s=float(ns),tau=float(tau),A_s=float(As),status="FAIL")
    try:
        bg=table(root+"00_background.dat")
        r["D_min"]=float(bg["kin (D)"].min()); r["cs2_min"]=float(bg["c_s^2"].min()); r["cs2_max"]=float(bg["c_s^2"].max())
        r["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1
        r["stable_subluminal"]=bool(stable)
        if stable:
            r.update(planck(root+"00_cl_lensed.dat"))
            s=sn(bg)
            for k,v in s.items():
                r[f"sn_{k}"]=float(v); r[f"sn_delta_{k}"]=float(v-LOCAL_SN[k])
            pks=sorted(OUT.glob(label+"_*pk.dat"))
            r["sigma8"]=sigma8(pks[0])
            i0=int(np.argmin(np.abs(bg.z.to_numpy()))); r["Omega_m0"]=float(bg.iloc[i0]["Omega_m(z)"])
            r["S8"]=r["sigma8"]*math.sqrt(r["Omega_m0"]/0.3)
            r["status"]="OK"
    except Exception as e:r["error"]=repr(e)
    if r["status"]!="OK":r["error"]=r.get("error",cp.stdout[-700:].replace("\n"," | "))
    print("LATE266_SUCCESSOR_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

anchors=[
 ("late266",AF0,ZC0,W0,DF0,NS0,TAU0),
 ("afm10",AF0*.90,ZC0,W0,DF0,NS0,TAU0),
 ("afm15",AF0*.85,ZC0,W0,DF0,NS0,TAU0),
 ("afm18",AF0*.82,ZC0,W0,DF0,NS0,TAU0),
]
rows=[evaluate(*a) for a in anchors]
center=[r for r in rows if r["id"]=="late266"][0]
if center["status"]!="OK": raise RuntimeError(center)
P0=float(center["chi2_planck_lite"]); S80=float(center["S8"])

# Box concentrates on the manuscript's Planck-improving AF edge but allows
# zc/width/Dfloor and primordial coordinates to recover stability and growth.
lo=np.array([AF0*.80,ZC0-0.20,W0-0.055,DF0-0.004,NS0-0.0020,0.044])
hi=np.array([AF0*.94,ZC0+0.22,W0+0.055,DF0+0.006,NS0+0.0040,0.057])
sob=qmc.Sobol(d=6,scramble=True,seed=2662026)
for i,x in enumerate(qmc.scale(sob.random_base2(m=6),lo,hi)):
    rows.append(evaluate(f"s{i:03d}",*map(float,x)))

df=pd.DataFrame(rows)
df["d_planck_lite"]=df["chi2_planck_lite"]-P0
df["d_S8"]=df["S8"]-S80
df["sn_worst"]=df[["sn_delta_pantheonplus","sn_delta_union3","sn_delta_desy5"]].max(axis=1)
# Promotion target: improve Planck, keep all 3 SN below local021, and do not
# sit directly on the scalar stability edge.
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
gate=ok[(ok.cs2_min>=0.006)&(ok.sn_worst<=0.0)&(ok.d_planck_lite<0.0)].copy()
if len(gate)<6:
    gate=ok[(ok.cs2_min>=0.004)&(ok.sn_worst<=0.03)&(ok.d_planck_lite<0.10)].copy()
gate["score"]=gate["d_planck_lite"] + 40*np.maximum(gate["sn_worst"],0) + 10*np.maximum(gate["d_S8"],0) + 0.02/np.maximum(gate["cs2_min"],1e-4)
prom=gate.sort_values(["score","d_planck_lite","S8"]).head(10)
df.to_csv(OUT/"late266_successor_scan.csv",index=False)
prom.to_csv(OUT/"late266_successor_promote.csv",index=False)
summary={
 "late266_center":{"planck_lite":P0,"S8":S80,
   "sn_deltas":{k:center[f"sn_delta_{k}"] for k in ["pantheonplus","union3","desy5"]}},
 "n_ok":int(len(ok)),"n_promote":int(len(prom)),
 "promotion_rule":"Full Plik + official DESI DR1 full shape + exact 3-SN replay required before any successor claim.",
 "promotion":prom.to_dict("records")
}
(OUT/"late266_successor_summary.json").write_text(json.dumps(summary,indent=2))
print("LATE266_SUCCESSOR_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
