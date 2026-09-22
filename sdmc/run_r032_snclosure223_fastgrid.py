#!/usr/bin/env python3
from pathlib import Path
import json, math, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure223_fastgrid"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5
H0=69.85635133907199; A_LATE=0.018589897081255913; B_LATE=0.013815921754576268
OB=0.022083219194622913; OC=0.12299536722293603
AF=0.03200255395658314; ZC=4.007449422683567; WIDTH=0.34799921004101636
D0=0.34231919445927034; DF=0.04607192634791136; LAM=18.40625; ZT=16.173189924377947
Q0=2.942463801247999; NS0=0.9625227132590487; TAU0=0.055202901571989066

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def load_cls(path):
 a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
 tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
 tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
 L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
 return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
 harr,dls=load_cls(path)
 def pc(A):
  ch=float(high.chi_squared(harr,A_planck=A))
  ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
  ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
  lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
  cl=float(-2*lens.log_likelihood(dls,**lp))
  cp=((A-1.)/CAL_SIGMA)**2
  return ch+ct+ce+cl+cp
 op=minimize_scalar(lambda A:pc(float(A)),bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
 return float(op.x),float(op.fun)

def ini(root,Q,ns,tau):
 As=math.exp(Q+2*tau)/1e10; h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
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
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},0.5,{A_LATE:.17g},0.25,{B_LATE:.17g},1.5
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

def run(tag,Q,ns,tau):
 root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,Q,ns,tau))
 cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
 rec={"id":tag,"Q":Q,"n_s":ns,"tau":tau,"A_s":math.exp(Q+2*tau)/1e10,"status":"FAIL"}
 cl=Path(root+"00_cl_lensed.dat")
 if cp.returncode==0 and cl.exists():
  A,chi=pscore(cl); rec.update(A_planck=A,chi2_planck=chi,status="OK")
 else: rec["error"]=cp.stdout[-600:].replace("\n"," | ")
 for p in OUT.glob(tag+"_*"):
  try:p.unlink()
  except:pass
 try:ip.unlink()
 except:pass
 print("FASTGRID_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec

rows=[]
for iq,dq in enumerate([-0.0015,0.0,0.0015]):
 for ins,dns in enumerate([-0.0030,0.0,0.0030]):
  for it,dt in enumerate([-0.0040,0.0,0.0040]):
   rows.append(run(f"g{iq}{ins}{it}",Q0+dq,NS0+dns,TAU0+dt))
center=run("center",Q0,NS0,TAU0); rows.append(center)
df=pd.DataFrame(rows); C0=float(center["chi2_planck"]); df["delta_vs_center"]=df.chi2_planck-C0
df.to_csv(OUT/"snclosure223_fastgrid.csv",index=False)
ok=df[df.status=="OK"].sort_values("chi2_planck")
summary={"center_chi2":C0,"best":ok.iloc[0].to_dict(),"best_delta_vs_center":float(ok.iloc[0].chi2_planck-C0),"top10":ok.head(10).to_dict("records")}
(OUT/"snclosure223_fastgrid_summary.json").write_text(json.dumps(summary,indent=2))
print("FASTGRID_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
