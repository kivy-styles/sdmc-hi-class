#!/usr/bin/env python3
"""
Ordinary-parameter reoptimization on the SN-closed sobol029 late geometry.
Vary (omega_b, omega_cdm, n_s, Q) while holding H0 and late A/B fixed.
Exact SN covariances and Planck-lite are evaluated for every stable point.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor,cho_solve
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure029_o4d_refine");OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
TCMB=2.7255;CAL=.0025;OR=4.17998772e-5
H0=69.71482083084993;AL=0.024822032167576252;BL=0.01264704344701022
OB0=0.022083219194622913;OC0=0.12299536722293603;NS0=0.9625227132590487
TAU=0.055202901571989066;Q0=2.942463801247999
AF=0.03200255395658314;ZC=4.007449422683567;W=0.34799921004101636;D0=0.34231919445927034;DF=0.04607192634791136
LAM=18.40625;ZT=16.173189924377947
LOCAL_SN={"pantheonplus":1406.2147033223882,"union3":28.783140002196888,"desy5":1650.6353947147727}
EDGE024_PD=-1.016724593277559
BASE_PLANCK_PENALTY=2.2215495955675805 # sobol029 relative to edge024 center in the late scan
BASE_BAO_DELTA=-2.375826707090633 # bao_chi2 sobol029 - edge024 bao chi2

def table(path):
 l=Path(path).read_text().splitlines();h=[x for x in l if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
 ms=list(re.finditer(r"(\d+)\s*:\s*",h));n=[]
 for i,m in enumerate(ms):n.append(h[m.end():(ms[i+1].start() if i+1<len(ms) else len(h))].strip())
 return pd.DataFrame(np.loadtxt(path),columns=n)
def rcov(path,n):
 a=np.asarray(np.loadtxt(path),float).reshape(-1)
 if len(a)==n*n+1 and int(round(a[0]))==n:a=a[1:]
 return a.reshape(n,n)
def setup(C):
 cf=cho_factor(C,lower=True,check_finite=False);one=np.ones(C.shape[0]);u=cho_solve(cf,one,check_finite=False);return cf,one,float(one@u)
def prof(pack,r):
 cf,one,den=pack;y=cho_solve(cf,r,check_finite=False);return float(r@y-(one@y)**2/den)

pp=pd.read_csv(SNROOT/'PantheonPlus/Pantheon+SH0ES.dat',sep=r'\s+',engine='python')
ppm=pp.m_b_corr.to_numpy(float);ppz=pp.zHD.to_numpy(float);ppzh=pp.zHEL.to_numpy(float);Cpp=rcov(SNROOT/'PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov',len(ppm))
m=ppz>0.01;ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)];pppack=setup(Cpp)
p=SNROOT/'Union3/lcparam_full.txt';cols=p.read_text().splitlines()[0].removeprefix('#').split()
un=pd.read_csv(p,sep=r'\s+',comment='#',names=cols,engine='python');cm={c.lower():c for c in un.columns};unm=un[cm['mb']].to_numpy(float);unz=un[cm['zcmb']].to_numpy(float);unpack=setup(rcov(SNROOT/'Union3/mag_covmat.txt',len(unm)))
de=pd.read_csv(SNROOT/'DESY5/DES-SN5YR_HD.csv');cm={c.lower():c for c in de.columns};dem=de[cm['mu']].to_numpy(float);dez=de[cm['zhd']].to_numpy(float);dezh=de[cm['zhel']].to_numpy(float);derr=de[cm['muerr_final']].to_numpy(float);depack=setup(rcov(SNROOT/'DESY5/covsys_000.txt',len(dem))+np.diag(derr*derr))

high=TTTEEE_lite_native(packages_path="planck_packages");lowT=TT(packages_path="planck_packages");lowE=EE(packages_path="planck_packages");lens=LensingNative(packages_path="planck_packages")
def mu(bg,z,zh):
 zz=bg.z.to_numpy();dm=bg['comov. dist.'].to_numpy();o=np.argsort(zz);DM=np.interp(z,zz[o],dm[o]);DA=DM/(1+z);return 5*np.log10((1+zh)*(1+z)*DA)
def sn(bg):
 return {"pantheonplus":prof(pppack,ppm-mu(bg,ppz,ppzh)),"union3":prof(unpack,unm-mu(bg,unz,unz)),"desy5":prof(depack,dem-mu(bg,dez,dezh))}
def pscore(path):
 a=np.loadtxt(path);ell=a[:,0].astype(int);n=int(ell.max())+1;conv=(TCMB*1e6)**2;L=ell.astype(float)
 tt=np.zeros(n);ee=np.zeros(n);te=np.zeros(n);pp=np.zeros(n);tt[ell]=a[:,1]*conv;ee[ell]=a[:,2]*conv;te[ell]=a[:,3]*conv;pp[ell]=a[:,5]*L*(L+1)
 harr=np.column_stack([ell,tt[ell],te[ell],ee[ell]]);d={"tt":tt,"te":te,"ee":ee,"pp":pp}
 def f(A):
  lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
  return float(high.chi_squared(harr,A_planck=A)-2*lowT.log_likelihood(tt,calib=A)-2*lowE.log_likelihood(ee,calib=A)-2*lens.log_likelihood(d,**lp)+((A-1)/CAL)**2)
 op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9});return float(op.fun),float(op.x)
def ini(root,ob,oc,ns,Q):
 As=math.exp(Q+2*TAU)/1e10;h=H0/100;ox=1-(ob+oc+OR)/(h*h)
 return textwrap.dedent(f"""H0={H0}
omega_b={ob}
omega_cdm={oc}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={As}
n_s={ns}
tau_reio={TAU}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={AF},{ZC},{W},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={ox},{LAM},{ZT},0.5,{AL},0.25,{BL},1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
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
def run(label,ob,oc,ns,Q):
 root=str(OUT/(label+"_"));ip=OUT/(label+".ini");ip.write_text(ini(root,ob,oc,ns,Q));cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
 r={"id":label,"omega_b":float(ob),"omega_cdm":float(oc),"n_s":float(ns),"Q":float(Q),"A_s":float(math.exp(Q+2*TAU)/1e10),"status":"FAIL"}
 bgp=Path(root+"00_background.dat");clp=Path(root+"00_cl_lensed.dat")
 if cp.returncode==0 and bgp.exists() and clp.exists():
  bg=table(bgp);r.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
  if r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1:
   s=sn(bg);r.update({f"sn_{k}":v for k,v in s.items()});r.update({f"sn_delta_{k}":s[k]-LOCAL_SN[k] for k in s});r["chi2_planck"],r["A_planck"]=pscore(clp);r["status"]="OK"
 print("SN029O4D_POINT",json.dumps(r,sort_keys=True),flush=True);return r
center=run("center",OB0,OC0,NS0,Q0);C0=center["chi2_planck"]
sob=qmc.Sobol(d=4,scramble=True,seed=94230);u=sob.random_base2(m=6);rows=[]
for i,x in enumerate(u):
 ob=OB0+(2*x[0]-1)*0.00045;oc=OC0+(2*x[1]-1)*0.0025;ns=NS0+(2*x[2]-1)*0.007;Q=Q0+(2*x[3]-1)*0.0035
 r=run(f"s{i:03d}",ob,oc,ns,Q)
 if r["status"]=="OK":
  r["delta_planck"]=r["chi2_planck"]-C0
  # conservative P+D proxy: baseline sobol029 P+D proxy + Planck movement only
  r["pd_proxy"]=0.6108683255173633+r["delta_planck"]
  r["joint_pp"]=r["pd_proxy"]+r["sn_delta_pantheonplus"];r["joint_u3"]=r["pd_proxy"]+r["sn_delta_union3"];r["joint_d5"]=r["pd_proxy"]+r["sn_delta_desy5"]
  r["second_sn"]=sorted([r["joint_pp"],r["joint_u3"],r["joint_d5"]])[1];r["goal"]=max(r["pd_proxy"],r["second_sn"])
 rows.append(r)
df=pd.DataFrame([center]+rows);df.to_csv(OUT/"snclosure029_o4d_refine.csv",index=False)
ok=pd.DataFrame(rows);ok=ok[ok.status=="OK"].sort_values(["goal","second_sn","pd_proxy"]);best=ok.head(12).to_dict("records")
print("SN029O4D_BEST",json.dumps(best,sort_keys=True),flush=True)
(OUT/"snclosure029_o4d_summary.json").write_text(json.dumps({"center":center,"best":best},indent=2))
