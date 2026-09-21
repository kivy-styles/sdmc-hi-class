#!/usr/bin/env python3
from pathlib import Path
import os, re, json, math, subprocess
import numpy as np

H0=float(os.environ["H0_SCAN"])
ZT=float(os.environ["ZT_SCAN"])
TAG=os.environ.get("TAG_SCAN",str(ZT).replace(".","p"))

OB=0.022083219194622913
OC=0.12299536722293603
OR=4.17998772e-5
AF=0.03773318954370916
ZC=3.6879637566395105
WIDTH=0.4136390263400972
D0=0.34231919445927034
DF=0.04570108858030289
LAM=18.40625
NS=0.9625227132590487
TAU=0.055202901571989066
Q=2.944463801247999
AS=math.exp(Q+2*TAU)/1e10
TCMB=2.7255
CAL_SIGMA=.0025

OUT=Path("output")/f"edge028_exact_zt_{TAG}"
OUT.mkdir(parents=True,exist_ok=True)
h=H0/100.
OX=1.-(OB+OC+OR)/(h*h)

def run(cmd,timeout=600):
    cp=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=timeout)
    print(cp.stdout,flush=True)
    if cp.returncode:
        raise RuntimeError(f"command failed {cmd}: rc={cp.returncode}")
    return cp

target=OUT/"target.ini"
target.write_text(f"""H0 = {H0:.17g}
omega_b = {OB:.17g}
omega_cdm = {OC:.17g}
N_ncdm = 0
N_ur = 3.046
T_cmb = 2.7255
YHe = 0.2453
A_s = {AS:.17e}
n_s = {NS:.17g}
tau_reio = {TAU:.17g}
Omega_Lambda = 0
Omega_fld = 3.1443554e-8
fluid_equation_of_state = SDMC_TRACKER
cs2_fld = 0.003
use_ppf = no
Omega_smg = -1
gravity_model = sdmc_v3_independent_kinetic
parameters_smg = {AF:.17g}, {ZC:.17g}, {WIDTH:.17g}, {D0:.17g}, 1.0, {DF:.17g}
expansion_model = sdmc_full
expansion_smg = {OX:.17g}, {LAM:.17g}, {ZT:.17g}, 0.5, 0.01105624999, 0.25, 0.01951933685, 1.5
pert_initial_conditions_smg = zero
method_qs_smg = fully_dynamic
a_ini_over_a_today_default = 1.e-8
a_ini_test_qs_smg = 1.e-8
pert_ic_ini_z_ref_smg = 1.e7
a_min_stability_test_smg = 1.e-8
output_background_smg = 3
write background = yes
write thermodynamics = no
root = {OUT}/cov_exact_target_
format = class
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
spectra_verbose = 0
output_verbose = 0
""")
run(["./class",str(target)],300)

# The reconstruction script reads output/cov_exact_target_00_background.dat.
src=OUT/"cov_exact_target_00_background.dat"
dst=Path("output/cov_exact_target_00_background.dat")
dst.write_bytes(src.read_bytes())

rp=Path("sdmc/run_r032_partial_zeq_ztbest_covariant.py")
s=rp.read_text()
s,n=re.subn(r"H0=[0-9.eE+-]+/299792\.458",f"H0={H0:.17g}/299792.458",s,count=1)
if n!=1: raise RuntimeError("failed to patch reconstruction H0")
rp.write_text(s)
run(["python3",str(rp)],300)

# The reconstructed model uses the target expansion history explicitly.
p=Path("source/background.c"); s=p.read_text()
pat=re.compile(
    r'double Ox\s*=\s*pba->parameters_smg\[0\];\s*'
    r'double lambda_e\s*=\s*pba->parameters_smg\[1\];\s*'
    r'double zt\s*=\s*pba->parameters_smg\[2\];\s*'
    r'double dNt\s*=\s*pba->parameters_smg\[3\];\s*'
    r'double A\s*=\s*pba->parameters_smg\[4\];\s*'
    r'double tauA\s*=\s*pba->parameters_smg\[5\];\s*'
    r'double B\s*=\s*pba->parameters_smg\[6\];\s*'
    r'double tauB\s*=\s*pba->parameters_smg\[7\];'
)
repl=f"""double Ox = {OX:.17g};
    double lambda_e = {LAM:.17g};
    double zt = {ZT:.17g};
    double dNt = 0.5;
    double A = 0.01105624999;
    double tauA = 0.25;
    double B = 0.01951933685;
    double tauB = 1.5;"""
s2,n=pat.subn(repl,s,count=1)
if n!=1: raise RuntimeError(f"tracker helper replacement count={n}")
p.write_text(s2)

run(["make","clean"],300)
run(["make","-j2","class"],600)

cand=OUT/"candidate.ini"
cand.write_text(f"""H0 = {H0:.17g}
omega_b = {OB:.17g}
omega_cdm = {OC:.17g}
N_ncdm = 0
N_ur = 3.046
T_cmb = 2.7255
YHe = 0.2453
A_s = {AS:.17e}
n_s = {NS:.17g}
tau_reio = {TAU:.17g}
Omega_Lambda = 0
Omega_fld = 3.1443554e-8
fluid_equation_of_state = SDMC_TRACKER
cs2_fld = 0.003
use_ppf = no
Omega_smg = -1
gravity_model = sdmc_v3_covariant_exact
parameters_smg = 0.0
pert_initial_conditions_smg = zero
method_qs_smg = fully_dynamic
a_ini_over_a_today_default = 1.e-8
a_ini_test_qs_smg = 1.e-8
pert_ic_ini_z_ref_smg = 1.e7
a_min_stability_test_smg = 1.e-8
output_background_smg = 3
modes = s
output = tCl,pCl,lCl
lensing = yes
l_max_scalars = 3000
write background = yes
root = {OUT}/covariant_
format = class
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
spectra_verbose = 0
lensing_verbose = 0
output_verbose = 0
""")
run(["./class",str(cand)],300)

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

bg=table(OUT/"covariant_00_background.dat")
mnD=float(np.min(bg["kin (D)"]))
mncs=float(np.min(bg["c_s^2"]))
mxcs=float(np.max(bg["c_s^2"]))
slip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
stable=(mnD>0 and mncs>0 and mxcs<=1)
if not stable:
    print("EXACT_ZT_RESULT",json.dumps(dict(tag=TAG,H0=H0,z_t=ZT,status="UNSTABLE",
          min_D=mnD,min_cs2=mncs,max_cs2=mxcs,max_abs_noslip=slip),sort_keys=True),flush=True)
    raise SystemExit(0)

from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

a=np.loadtxt(OUT/"covariant_00_cl_lensed.dat")
ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
harr=np.column_stack([ell,tt[ell],te[ell],ee[ell]])
dls={"tt":tt,"te":te,"ee":ee,"pp":pp}

def pieces(A):
    ch=float(high.chi_squared(harr,A_planck=A))
    ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
    ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
    lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
    cl=float(-2*lens.log_likelihood(dls,**lp))
    cp=((A-1.)/CAL_SIGMA)**2
    return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp

op=minimize_scalar(lambda A:pieces(float(A))[-1],bounds=(.97,1.03),method="bounded",
                   options={"xatol":1e-10})
A=float(op.x); z=pieces(A)
rec=dict(tag=TAG,H0=H0,z_t=ZT,status="OK",A_planck=A,
         chi2_high=z[0],chi2_lowT=z[1],chi2_lowE=z[2],chi2_lensing=z[3],
         chi2_cal=z[4],chi2_planck=z[5],min_D=mnD,min_cs2=mncs,max_cs2=mxcs,
         max_abs_noslip=slip,Omega_x=OX)
(OUT/"result.json").write_text(json.dumps(rec,indent=2))
print("EXACT_ZT_RESULT",json.dumps(rec,sort_keys=True),flush=True)
