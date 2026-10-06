#!/usr/bin/env python3
"""
A034 <-> late266 background corridor.

Fix the successful g022-a034 perturbative/structural sector and continuously
interpolate only the ordinary + late-background geometry from late266 (t=0)
to the exact a034/g022 background (t=1).  This searches for a point that
retains late266's stronger SN distances while inheriting enough of a034's
very strong official DESI full-shape advantage to cross the joint gates.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/a034_late266_corridor"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data"); BAOROOT=Path("bao_data")
T=2.7255; SIG=.0025; OR=4.17998772e-5

# late266 endpoint
L=dict(H0=69.45160505326464,ob=0.02208511574370519,oc=0.12298509428428875,
       zt=16.173189924377947,Alate=0.013036667098291217,Blate=0.010544837576337158)
# g022/a034 endpoint
G=dict(H0=69.71482083084993,ob=0.022083219194622913,oc=0.12299536722293603,
       zt=16.19317622240633,Alate=0.024822032167576252,Blate=0.01264704344701022)

# a034 perturbative/structural sector
AS=2.1047019076951746e-09; NS=0.9628416801318527; TAU=0.053167105823755265
AF=0.02048146490100771; ZC=3.927876388467848; WIDTH=0.33114133956842123
D0=0.34231919445927034; POWER=1.0; DF=0.05181542627513409
LAM=18.40625

LOCAL_SN={"pantheonplus":1406.2147033223882,"union3":28.783140002196888,"desy5":1650.6353947147727}
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
  return a.reshape(n,n)
def setup(C):
  cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0]); u=cho_solve(cf,one,check_finite=False)
  return cf,one,float(one@u)
def prof(pack,r):
  cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
  return float(r@y-(one@y)**2/den)

pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm)); m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=setup(Cpp)
p=SNROOT/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unp=setup(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))
de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float); derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))
def mu(bg,z,zh):
  DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z); return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
  return {"pantheonplus":prof(ppp,ppm-mu(bg,ppz,ppzh)),
          "union3":prof(unp,unm-mu(bg,unz,unz)),
          "desy5":prof(dep,dem-mu(bg,dez,dezh))}

def load_cls(path):
  a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(T*1e6)**2
  tt=np.zeros(n);ee=np.zeros(n);te=np.zeros(n);ppv=np.zeros(n)
  tt[ell]=a[:,1]*conv;ee[ell]=a[:,2]*conv;te[ell]=a[:,3]*conv;L=ell.astype(float);ppv[ell]=a[:,5]*L*(L+1)
  return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":ppv}
def pscore(path):
  h,s=load_cls(path)
  def f(A):
    v=float(high.chi_squared(h,A_planck=A))+float(-2*lowT.log_likelihood(s["tt"],calib=A))+float(-2*lowE.log_likelihood(s["ee"],calib=A))
    lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}; v+=float(-2*lens.log_likelihood(s,**lp)); v+=((A-1)/SIG)**2; return v
  op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9}); return float(op.fun),float(op.x)
def sig8(path):
  a=np.loadtxt(path);k=a[:,0];P=a[:,1];x=8*k;W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
  return float(np.sqrt(np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))))

def ini(root,c):
  h=c["H0"]/100.; ox=1-(c["ob"]+c["oc"]+OR)/(h*h)
  return textwrap.dedent(f"""\
H0={c['H0']:.17g}
omega_b={c['ob']:.17g}
omega_cdm={c['oc']:.17g}
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
parameters_smg={AF},{ZC},{WIDTH},{D0},{POWER},{DF}
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
  c={k:L[k]+t*(G[k]-L[k]) for k in L}; label=f"c{int(round(100*t)):03d}"; root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,c))
  cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=360)
  r=dict(id=label,t=float(t),**{k:float(v) for k,v in c.items()},status="FAIL")
  try:
    bg=table(root+"00_background.dat"); r["Dmin"]=float(bg["kin (D)"].min());r["cs2min"]=float(bg["c_s^2"].min());r["cs2max"]=float(bg["c_s^2"].max())
    if cp.returncode==0 and r["Dmin"]>0 and r["cs2min"]>0 and r["cs2max"]<=1:
      pc,Ap=pscore(root+"00_cl_lensed.dat"); ss=sns(bg); pk=sorted(OUT.glob(label+"_*pk.dat"))[0]; s=sig8(pk); i0=int(np.argmin(abs(bg.z.to_numpy())));Om=float(bg.iloc[i0]["Omega_m(z)"])
      r.update(status="OK",planck_lite=pc,A_planck=Ap,S8=s*math.sqrt(Om/0.3),
               sn_pp_delta=ss["pantheonplus"]-LOCAL_SN["pantheonplus"],
               sn_u3_delta=ss["union3"]-LOCAL_SN["union3"],sn_d5_delta=ss["desy5"]-LOCAL_SN["desy5"])
  except Exception as e:r["error"]=repr(e)
  print("A034_LATE266_CORRIDOR",json.dumps(r,sort_keys=True),flush=True)
  for p in OUT.glob(label+"_*"):
    try:p.unlink()
    except:pass
  try:ip.unlink()
  except:pass
  return r
rows=[run(float(t)) for t in np.linspace(0,1,21)]
df=pd.DataFrame(rows);df.to_csv(OUT/"corridor.csv",index=False)
# Keep points spanning SN-negative through a034-like geometry; exact DESI decides.
ok=df[df.status=="OK"].copy()
ok["sn_worst"]=ok[["sn_pp_delta","sn_u3_delta","sn_d5_delta"]].max(axis=1)
prom=ok.sort_values(["sn_worst","planck_lite"]).head(12)
# Ensure mid/high-t points are represented because DESI is expected to improve toward a034.
extra=ok.iloc[[8,10,12,14,16,18,20]] if len(ok)>=21 else ok.tail(7)
prom=pd.concat([prom,extra]).drop_duplicates("id").sort_values("t")
prom.to_csv(OUT/"promotion.csv",index=False)
summary={"promotion":prom.to_dict("records"),"rule":"Exact DESI + full Plik + separate SN determine the best joint point; SN-only negativity is not required."}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("A034_LATE266_CORRIDOR_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
