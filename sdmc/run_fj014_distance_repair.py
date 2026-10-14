#!/usr/bin/env python3
"""
fj014 distance-sector repair.
Freeze the strong fj014 structural/kinetic/primordial coordinates and move only
the late/background distance sector from fj014 back toward late266, with a small
controlled extrapolation beyond late266.  This directly targets the remaining
Pantheon+/Union3/DES-Y5 deficits while measuring the Planck-lite/BAO/age cost.
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

OUT=Path("output/fj014_distance_repair"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data"); BAOROOT=Path("bao_data")
T=2.7255; SIG=.0025; OR=4.17998772e-5; LAM=18.40625
LOCAL_SN={"pantheonplus":1406.2147033223882,"union3":28.783140002196888,"desy5":1650.6353947147727}
LOCAL_BAO=13.400290718544086; LOCAL_AGE=13.75343

# Freeze fj014 structural / kinetic / primordial sector.
FIX=dict(As=2.1096255943551658e-9,ns=0.9639264293648303,tau=0.05461631424725056,
         AF=0.01948851037286222,zc=3.644913007412106,width=0.35893273874651643,
         Dfloor=0.05569390742853284,D0=0.3865188211016357,power=1.0877414136193693)

# fj014 distance sector -> late266 distance sector.
A=dict(H0=69.51135345487855,ob=0.0221336738740094,oc=0.12258112205294892,
       zt=16.06170386472717,Alate=0.011254649444762618,Blate=0.00994506031619385)
B=dict(H0=69.45160505326464,ob=0.02208511574370519,oc=0.12298509428428875,
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
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:a=a[1:]
    if len(a)!=n*n: raise RuntimeError(path)
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0]); u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one@u)
def pchi(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False); return float(r@y-(one@y)**2/den)

pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
C=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm)); m=ppz>0.01
ppm,ppz,ppzh,C=ppm[m],ppz[m],ppzh[m],C[np.ix_(m,m)]; ppp=setup(C)
p=SNROOT/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unp=setup(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))
de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float); derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

bao_rows=[]
for ln in (BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
    if ln.strip() and not ln.startswith("#"):
        z,v,q=ln.split(); bao_rows.append((float(z),float(v),q))
bao_cov=np.loadtxt(BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_cov.txt"); bao_inv=np.linalg.inv(bao_cov); bao_obs=np.array([x[1] for x in bao_rows])

def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sn(bg):
    return {"pantheonplus":pchi(ppp,ppm-mu(bg,ppz,ppzh)),
            "union3":pchi(unp,unm-mu(bg,unz,unz)),
            "desy5":pchi(dep,dem-mu(bg,dez,dezh))}
def rd(bg,th):
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy()
    for i in range(len(zz)-1):
        if (kb[i]-1)*(kb[i+1]-1)<=0 and kb[i]!=kb[i+1]:
            q=(1-kb[i])/(kb[i+1]-kb[i]); zd=zz[i]+q*(zz[i+1]-zz[i]); return float(np.interp(zd,bg.z,bg["comov.snd.hrz."])),float(zd)
    raise RuntimeError("drag")
def bao(bg,r):
    pred=[]
    for z,_,q in bao_rows:
        DM=float(np.interp(z,bg.z,bg["comov. dist."])); H=float(np.interp(z,bg.z,bg["H [1/Mpc]"])); DH=1/H; DV=(z*DM*DM*DH)**(1/3)
        pred.append({"DM_over_rs":DM/r,"DH_over_rs":DH/r,"DV_over_rs":DV/r}[q])
    d=np.array(pred)-bao_obs; return float(d@bao_inv@d)

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(T*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv; L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}
def planck(path):
    h,s=load_cls(path)
    def f(Ap):
        v=float(high.chi_squared(h,A_planck=Ap)); v+=float(-2*lowT.log_likelihood(s["tt"],calib=Ap)); v+=float(-2*lowE.log_likelihood(s["ee"],calib=Ap))
        lp={lens.calibration_param:Ap} if getattr(lens,"calibration_param",None) else {}; v+=float(-2*lens.log_likelihood(s,**lp)); v+=((Ap-1)/SIG)**2
        return v
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9}); return float(op.fun),float(op.x)

def ini(root,c):
    h=c["H0"]/100.; ox=1.-(c["ob"]+c["oc"]+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={c['H0']:.17g}
omega_b={c['ob']:.17g}
omega_cdm={c['oc']:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={FIX['As']:.17e}
n_s={FIX['ns']:.17g}
tau_reio={FIX['tau']:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={FIX['AF']:.17g},{FIX['zc']:.17g},{FIX['width']:.17g},{FIX['D0']:.17g},{FIX['power']:.17g},{FIX['Dfloor']:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM:.17g},{c['zt']:.17g},0.5,{c['Alate']:.17g},0.25,{c['Blate']:.17g},1.5
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
write thermodynamics=yes
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
    c={k:A[k]+t*(B[k]-A[k]) for k in A}; label=f"r{int(round(t*100)):03d}"; root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
    out=dict(id=label,t=float(t),**{k:float(v) for k,v in c.items()},status="FAIL")
    try:
        bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat"); D=bg["kin (D)"].to_numpy(); cs=bg["c_s^2"].to_numpy()
        out["D_min"]=float(D.min()); out["cs2_min"]=float(cs.min()); out["cs2_max"]=float(cs.max()); out["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        if cp.returncode==0 and out["D_min"]>0 and out["cs2_min"]>0 and out["cs2_max"]<=1:
            pc,Ap=planck(root+"00_cl_lensed.dat"); ss=sn(bg); rr,zd=rd(bg,th); bc=bao(bg,rr); age=float(np.interp(0,bg.z,bg["proper time [Gyr]"]))
            out.update(status="OK",planck_lite=pc,A_planck=Ap,sn_pp_delta=ss["pantheonplus"]-LOCAL_SN["pantheonplus"],sn_u3_delta=ss["union3"]-LOCAL_SN["union3"],sn_d5_delta=ss["desy5"]-LOCAL_SN["desy5"],bao_delta=bc-LOCAL_BAO,rd=rr,z_drag=zd,age=age)
    except Exception as e: out["error"]=repr(e)
    print("FJ014_REPAIR_POINT",json.dumps(out,sort_keys=True),flush=True)
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return out

rows=[run(float(t)) for t in np.linspace(-0.25,1.25,31)]
df=pd.DataFrame(rows)
base=float(df[np.isclose(df.t,0)].iloc[0].planck_lite)
df["d_planck_vs_fj014"]=df["planck_lite"]-base
df["sn_worst"]=df[["sn_pp_delta","sn_u3_delta","sn_d5_delta"]].max(axis=1)
df["joint_proxy"]=df["d_planck_vs_fj014"]+df["sn_worst"]+0.1*np.maximum(df["bao_delta"],0)
df.to_csv(OUT/"fj014_distance_repair.csv",index=False)
prom=df[(df.status=="OK")&(df.cs2_min>=0.006)].sort_values(["joint_proxy","planck_lite"]).head(10)
prom.to_csv(OUT/"fj014_distance_repair_promote.csv",index=False)
summary={"promotion":prom.to_dict("records"),"rule":"Promote best points to exact full Plik + DESI + separate SN."}
(OUT/"fj014_distance_repair_summary.json").write_text(json.dumps(summary,indent=2))
print("FJ014_REPAIR_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
