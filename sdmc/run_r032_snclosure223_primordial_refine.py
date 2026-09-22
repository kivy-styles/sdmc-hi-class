#!/usr/bin/env python3
"""
Primordial re-closure around SN-aware sobol223 late-background candidate.

Late geometry and SDMC structure are frozen at the SN-improved point:
  A=0.018589897081255913
  B=0.013815921754576268
  H0=69.85635133907199
with the successful edge024 structural action.

We vary only (Q, n_s, tau), where Q=ln(1e10 A_s)-2 tau.
These directions leave the SN luminosity-distance curve unchanged, making
this the cleanest attempt to recover the Planck penalty without sacrificing
the SN gains.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure223_primordial_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0=69.85635133907199
A_LATE=0.018589897081255913
B_LATE=0.013815921754576268
OB=0.022083219194622913
OC=0.12299536722293603

AF=0.03200255395658314
ZC=4.007449422683567
WIDTH=0.34799921004101636
D0=0.34231919445927034
DF=0.04607192634791136
LAM=18.40625
ZT=16.173189924377947
DNT=0.5
TAUA=0.25
TAUB=1.5

Q0=2.942463801247999
NS0=0.9625227132590487
TAU0=0.055202901571989066

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
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
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,Q,ns,tau):
    As=math.exp(Q+2*tau)/1e10
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
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},{DNT},{A_LATE:.17g},{TAUA},{B_LATE:.17g},{TAUB}
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

def run(label,Q,ns,tau):
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,Q,ns,tau))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=label,Q=float(Q),n_s=float(ns),tau=float(tau),
             A_s=float(math.exp(Q+2*tau)/1e10),status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        rec.update(pscore(clp)); rec["status"]="OK"
    else:
        rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    print("SN223PRIM_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

center=run("center",Q0,NS0,TAU0)
if center["status"]!="OK": raise RuntimeError(center)
C0=center["chi2_planck"]

# Stage A: broad 64-point Sobol box.
sob=qmc.Sobol(d=3,scramble=True,seed=94223)
u=sob.random_base2(m=6)
broad=[]
for i,x in enumerate(u):
    Q=Q0 + (2*x[0]-1)*0.0045
    ns=NS0 + (2*x[1]-1)*0.0100
    tau=TAU0 + (2*x[2]-1)*0.0120
    broad.append(run(f"b{i:03d}",Q,ns,tau))
bdf=pd.DataFrame(broad)
bok=bdf[bdf.status=="OK"].sort_values("chi2_planck")
best0=bok.iloc[0]

# Stage B: local 64-point refinement around broad winner.
sob2=qmc.Sobol(d=3,scramble=True,seed=94224)
u2=sob2.random_base2(m=6)
fine=[]
for i,x in enumerate(u2):
    Q=float(best0.Q)+(2*x[0]-1)*0.0015
    ns=float(best0.n_s)+(2*x[1]-1)*0.0035
    tau=float(best0.tau)+(2*x[2]-1)*0.0040
    # physical/simple bounds
    Q=min(max(Q,Q0-0.007),Q0+0.007)
    ns=min(max(ns,0.945),0.980)
    tau=min(max(tau,0.035),0.080)
    fine.append(run(f"f{i:03d}",Q,ns,tau))

allrows=[center]+broad+fine
df=pd.DataFrame(allrows)
df["delta_vs_center"]=df.chi2_planck-C0
df.to_csv(OUT/"snclosure223_primordial_refine.csv",index=False)
ok=df[df.status=="OK"].sort_values("chi2_planck")
best=ok.head(12).to_dict("records")
summary={
 "late_candidate":{"A":A_LATE,"B":B_LATE,"H0":H0},
 "center":{"Q":Q0,"n_s":NS0,"tau":TAU0,"chi2_planck":C0},
 "best":best[0],
 "best_delta_vs_center":float(best[0]["chi2_planck"]-C0),
 "top12":best,
 "required_recovery_note":"sobol223 needed about 0.63 chi2 additional P+D improvement to make Pantheon+ and Union3 jointly negative under the conservative first-stage proxy"
}
(OUT/"snclosure223_primordial_refine_summary.json").write_text(json.dumps(summary,indent=2))
print("SN223PRIM_BEST",json.dumps(best,sort_keys=True),flush=True)
print("SN223PRIM_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
