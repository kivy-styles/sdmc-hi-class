#!/usr/bin/env python3
"""
FJ014 background bridge.

Hold the successful fj014 perturbative/structural sector fixed and interpolate
only the ordinary + late distance sector from fj014 toward late266. This asks
whether the very small fj014 SN penalties can be removed before the exact
Planck+DESI gain is lost.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/fj014_background_bridge"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data"); T=2.7255; SIG=.0025; OR=4.17998772e-5; LAM=18.40625
LOCAL_SN={"pantheonplus":1406.2147033223882,"union3":28.783140002196888,"desy5":1650.6353947147727}

# Freeze fj014 perturbative/structural sector.
P=dict(As=2.1096255943551658e-9,ns=0.9639264293648303,tau=0.05461631424725056,
       AF=0.01948851037286222,zc=3.644913007412106,width=0.35893273874651643,
       Dfloor=0.05569390742853284,D0=0.3865188211016357,power=1.0877414136193693)
FJ=dict(H0=69.51135345487855,ob=0.0221336738740094,oc=0.12258112205294892,
        zt=16.06170386472717,Alate=0.011254649444762618,Blate=0.00994506031619385)
L266=dict(H0=69.45160505326464,ob=0.02208511574370519,oc=0.12298509428428875,
          zt=16.173189924377947,Alate=0.013036667098291217,Blate=0.010544837576337158)

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")
def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(T*1e6)**2
    tt=np.zeros(n);ee=np.zeros(n);te=np.zeros(n);pp=np.zeros(n)
    tt[ell]=a[:,1]*conv;ee[ell]=a[:,2]*conv;te[ell]=a[:,3]*conv
    L=ell.astype(float);pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}
def pscore(path):
    h,s=load_cls(path)
    def f(A):
        v=float(high.chi_squared(h,A_planck=A))
        v+=float(-2*lowT.log_likelihood(s["tt"],calib=A))
        v+=float(-2*lowE.log_likelihood(s["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        v+=float(-2*lens.log_likelihood(s,**lp))+((A-1)/SIG)**2
        return v
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:a=a[1:]
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False);one=np.ones(C.shape[0]);u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one@u)
def chi(pack,r):
    cf,one,den=pack;y=cho_solve(cf,r,check_finite=False)
    return float(r@y-(one@y)**2/den)
pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",engine="python")
ppm=pp["m_b_corr"].to_numpy(float);ppz=pp["zHD"].to_numpy(float);ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm));m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=setup(Cpp)
p=SNROOT/"Union3/lcparam_full.txt";cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python");cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float);unz=un[cm["zcmb"]].to_numpy(float);unp=setup(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))
de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv");cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float);dez=de[cm["zhd"]].to_numpy(float);dezh=de[cm["zhel"]].to_numpy(float);derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))
def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]);DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
    return dict(pantheonplus=chi(ppp,ppm-mu(bg,ppz,ppzh)),
                union3=chi(unp,unm-mu(bg,unz,unz)),
                desy5=chi(dep,dem-mu(bg,dez,dezh)))
def sig8(path):
    a=np.loadtxt(path);k=a[:,0];P=a[:,1];x=8*k
    W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
    return float(np.sqrt(np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))))
def ini(root,c):
    h=c["H0"]/100.;ox=1.-(c["ob"]+c["oc"]+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={c['H0']:.17g}
omega_b={c['ob']:.17g}
omega_cdm={c['oc']:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={P['As']:.17e}
n_s={P['ns']:.17g}
tau_reio={P['tau']:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={P['AF']},{P['zc']},{P['width']},{P['D0']},{P['power']},{P['Dfloor']}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{c['zt']:.17g},0.5,{c['Alate']:.17g},0.25,{c['Blate']:.17g},1.5
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
def run(t):
    c={k:FJ[k]+t*(L266[k]-FJ[k]) for k in FJ}
    label=f"b{int(round(t*1000)):04d}";root=str(OUT/(label+"_"));ip=OUT/(label+".ini");ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    r=dict(id=label,t=float(t),**{k:float(v) for k,v in c.items()},status="FAIL")
    try:
        bg=table(root+"00_background.dat"); r["D_min"]=float(bg["kin (D)"].min());r["cs2_min"]=float(bg["c_s^2"].min());r["cs2_max"]=float(bg["c_s^2"].max())
        stable=cp.returncode==0 and r["D_min"]>0 and r["cs2_min"]>0 and r["cs2_max"]<=1
        if stable:
            pc,Ap=pscore(root+"00_cl_lensed.dat");s=sns(bg);pk=sorted(OUT.glob(label+"_*pk.dat"));sg=sig8(pk[0])
            i0=int(np.argmin(np.abs(bg.z.to_numpy())));Om=float(bg.iloc[i0]["Omega_m(z)"])
            r.update(status="OK",planck_lite=pc,A_planck=Ap,S8=sg*math.sqrt(Om/0.3),sigma8=sg,
                     sn_pp_delta=s["pantheonplus"]-LOCAL_SN["pantheonplus"],
                     sn_u3_delta=s["union3"]-LOCAL_SN["union3"],
                     sn_d5_delta=s["desy5"]-LOCAL_SN["desy5"])
    except Exception as e:r["error"]=repr(e)
    print("FJ014_BG_POINT",json.dumps(r,sort_keys=True),flush=True)
    for q in OUT.glob(label+"_*"):
        try:q.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r
rows=[run(float(t)) for t in np.linspace(0,1,41)]
df=pd.DataFrame(rows);df.to_csv(OUT/"fj014_bg_bridge.csv",index=False)
# Favor near-zero SN with best Planck screen; exact DESI decides final closure.
df["sn_worst"]=df[["sn_pp_delta","sn_u3_delta","sn_d5_delta"]].max(axis=1)
ok=df[(df.status=="OK")&(df.cs2_min>=0.006)]
prom=pd.concat([
    ok.sort_values(["sn_worst","planck_lite"]).head(8),
    ok[ok.sn_worst<=0.05].sort_values("planck_lite").head(8),
    ok[ok.sn_worst<=0.0].sort_values("planck_lite").head(8)
]).drop_duplicates("id").head(16)
prom.to_csv(OUT/"fj014_bg_promote.csv",index=False)
summary={"promotion":prom.to_dict("records"),
         "rule":"Exact full Plik + official DESI must be evaluated; objective is P+D+each SN against local021, not SN-only."}
(OUT/"fj014_bg_summary.json").write_text(json.dumps(summary,indent=2))
print("FJ014_BG_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
