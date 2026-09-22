#!/usr/bin/env python3
from pathlib import Path
import json,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative
OUT=Path("output/snclosure029_afthreshold"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL=.0025; OR=4.17998772e-5
H0=69.71482083084993; AL=0.024822032167576252; BL=0.01264704344701022
OB=0.022083219194622913; OC=0.12299536722293603
AS=2.117602412937069e-9; NS=0.9625227132590487; TAU=0.055202901571989066
AF0=0.03200255395658314; ZC=4.007449422683567; W=0.34799921004101636
D0=0.34231919445927034; DF=0.04607192634791136; LAM=18.40625; ZT=16.173189924377947
high=TTTEEE_lite_native(packages_path="planck_packages"); lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages"); lens=LensingNative(packages_path="planck_packages")
def table(p):
 l=Path(p).read_text().splitlines(); h=[x for x in l if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
 ms=list(re.finditer(r"(\d+)\s*:\s*",h)); n=[]
 for i,m in enumerate(ms): n.append(h[m.end():(ms[i+1].start() if i+1<len(ms) else len(h))].strip())
 return pd.DataFrame(np.loadtxt(p),columns=n)
def score(p):
 a=np.loadtxt(p); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
 tt=np.zeros(n);ee=np.zeros(n);te=np.zeros(n);pp=np.zeros(n); L=ell.astype(float)
 tt[ell]=a[:,1]*conv;ee[ell]=a[:,2]*conv;te[ell]=a[:,3]*conv;pp[ell]=a[:,5]*L*(L+1)
 harr=np.column_stack([ell,tt[ell],te[ell],ee[ell]]); d={"tt":tt,"te":te,"ee":ee,"pp":pp}
 def f(A):
  lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
  return float(high.chi_squared(harr,A_planck=A)-2*lowT.log_likelihood(tt,calib=A)-2*lowE.log_likelihood(ee,calib=A)-2*lens.log_likelihood(d,**lp)+((A-1)/CAL)**2)
 op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
 return float(op.fun),float(op.x)
def ini(root,AF):
 h=H0/100; ox=1-(OB+OC+OR)/(h*h)
 return textwrap.dedent(f"""H0={H0}
omega_b={OB}
omega_cdm={OC}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={AS}
n_s={NS}
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
rows=[]
for i,AF in enumerate(np.linspace(0.0205,0.0210,11)):
 root=str(OUT/f"a{i:02d}_"); ip=OUT/f"a{i:02d}.ini"; ip.write_text(ini(root,AF))
 cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
 r={"id":f"a{i:02d}","A_F":float(AF),"status":"FAIL"}
 bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
 if cp.returncode==0 and bgp.exists() and clp.exists():
  bg=table(bgp); r.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
  if r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1:
   r["chi2_planck"],r["A_planck"]=score(clp);r["status"]="OK"
 rows.append(r); print("SN029AFTH_POINT",json.dumps(r,sort_keys=True),flush=True)
df=pd.DataFrame(rows); center=float(df.iloc[np.argmin(abs(df.A_F-AF0))].chi2_planck)
df["delta_vs_nearest_center"]=df.chi2_planck-center; df.to_csv(OUT/"snclosure029_afthreshold.csv",index=False)
ok=df[df.status=="OK"].sort_values("chi2_planck"); best=ok.head(8).to_dict("records")
print("SN029AFTH_BEST",json.dumps(best,sort_keys=True),flush=True)
