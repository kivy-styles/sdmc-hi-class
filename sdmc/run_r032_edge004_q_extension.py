#!/usr/bin/env python3
"""
Fine Q re-closure around the edge004 structural winner.

All geometry, matter, the edge004 SDMC structure, n_s and tau are fixed.
Only Q = ln(10^10 A_s) - 2 tau is varied. The purpose is to re-optimize
the primordial-amplitude combination after the structural minimum moved.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/edge004_q_extension"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0=70.78157258819029
OB=0.022083219194622913
OC=0.12299536722293603
NS=0.9625227132590487
TAU=0.055202901571989066
Q0=2.944463801247999

AF=0.034172075465787204
ZC=3.943745281267911
WIDTH=0.36045085035264496
D0=0.34231919445927034
DF=0.0479059932352975
LAM=18.40625
ZT=15.418161678034812

OFFSETS=np.array([-0.0080,-0.0075,-0.0070,-0.0065,-0.0060,-0.0055,
                   -0.0050,-0.0045,-0.0040,-0.0035,-0.0030,-0.0025,
                   -0.0020,0.0])

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
    return pd.DataFrame(np.loadtxt(path),columns=names)

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
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,Q):
    lnAs=Q+2.*TAU
    As=math.exp(lnAs)/1e10
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
n_s={NS:.17g}
tau_reio={TAU:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
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

def run(i,dq):
    Q=Q0+dq; root=str(OUT/(f"q{i:02d}_")); ip=OUT/(f"q{i:02d}.ini")
    ip.write_text(ini(root,Q))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=f"q{i:02d}",dq=float(dq),Q=float(Q),ln10As=float(Q+2*TAU),
             A_s=float(math.exp(Q+2*TAU)/1e10),status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1)
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(pscore(clp)); rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-800:].replace("\n"," | ")
    print("E04QX_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[run(i,dq) for i,dq in enumerate(OFFSETS)]
df=pd.DataFrame(rows)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
center=float(ok.loc[np.isclose(ok.dq,0.0),"chi2_planck"].iloc[0])
ok["delta_vs_edge004"]=ok.chi2_planck-center
df=pd.DataFrame(rows).merge(ok[["id","delta_vs_edge004"]],on="id",how="left")
df.to_csv(OUT/"edge004_q_extension.csv",index=False)
best=ok.nsmallest(8,"chi2_planck").to_dict("records")
summary={"Q0":Q0,"center_chi2":center,"best":best[0],
         "best_delta_vs_edge004":float(ok.chi2_planck.min()-center),
         "left_edge_delta":float(ok.sort_values("Q").iloc[0].chi2_planck-center),
         "right_edge_delta":float(ok.sort_values("Q").iloc[-1].chi2_planck-center)}
(OUT/"edge004_q_extension_summary.json").write_text(json.dumps(summary,indent=2))
print("E04QX_BEST",json.dumps(best,sort_keys=True),flush=True)
print("E04QX_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
