#!/usr/bin/env python3
"""
Structural No-Slip refinement around exact SN-close010.

Late SN-friendly background and ordinary parameters are frozen:
  A_late=0.062007905058562754, B_late=0.018307528030127286,
  H0=70.31907490585.

Reopen only (A_F,z_c,width,D_floor), which changes the covariant/perturbation
sector but not the prescribed late distance history.  The goal is to recover
>=~0.46 chi2 in Planck+DESI while preserving SN-close010's SN shape.
This is a Planck-lite structural screen; winners require exact full-Plik
and raw DESI promotion.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclose010_structure_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0=70.31907490585
OB=0.022083219194622913
OC=0.12299536722293603
NS=0.9625227132590487
TAU=0.055202901571989066
Q=2.942463801247999
AS=math.exp(Q+2*TAU)/1e10

AF0=0.03200255395658314
ZC0=4.007449422683567
W0=0.34799921004101636
D0=0.34231919445927034
DF0=0.04607192634791136

LAM=18.40625
ZT=16.173189924377947
A_LATE=0.062007905058562754
B_LATE=0.018307528030127286
OX=1.-(OB+OC+OR)/(H0/100.)**2

LOW=np.array([0.0240,3.72,0.300,0.0430],float)
HIGH=np.array([0.0400,4.24,0.415,0.0550],float)

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
                       options={"xatol":1e-9})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,af,zc,w,df):
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
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
parameters_smg={af:.17g},{zc:.17g},{w:.17g},{D0:.17g},1.0,{df:.17g}
expansion_model=sdmc_full
expansion_smg={OX:.17g},{LAM:.17g},{ZT:.17g},0.5,{A_LATE:.17g},0.25,{B_LATE:.17g},1.5
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

def run(tag,af,zc,w,df):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,af,zc,w,df))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=tag,A_F=float(af),z_c=float(zc),width=float(w),D_floor=float(df),
             status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(pscore(clp)); rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-700:].replace("\n"," | ")
    # keep only summaries
    for p in OUT.glob(tag+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    print("SN010STRUCT_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[]
base=run("center",AF0,ZC0,W0,DF0)
if base.get("status")!="OK": raise RuntimeError(f"center failed {base}")
BASE=float(base["chi2_planck"]); rows.append(base)

anchors=[
 ("afm",AF0-0.002,ZC0,W0,DF0),("afp",AF0+0.002,ZC0,W0,DF0),
 ("zcm",AF0,ZC0-0.12,W0,DF0),("zcp",AF0,ZC0+0.12,W0,DF0),
 ("wm",AF0,ZC0,W0-0.025,DF0),("wp",AF0,ZC0,W0+0.025,DF0),
 ("dfm",AF0,ZC0,W0,DF0-0.002),("dfp",AF0,ZC0,W0,DF0+0.003)]
for a in anchors: rows.append(run(*a))

sam=qmc.Sobol(d=4,scramble=True,seed=10066)
pts=qmc.scale(sam.random_base2(m=7),LOW,HIGH)
for i,p in enumerate(pts): rows.append(run(f"sobol{i:03d}",*map(float,p)))

for r in rows:
    if r.get("status")=="OK":
        r["delta_planck_vs_center"]=r["chi2_planck"]-BASE
        r["required_exact_pd_recovery_margin"]=-(r["delta_planck_vs_center"])-0.1150406719710304

df=pd.DataFrame(rows); df.to_csv(OUT/"snclose010_structure_refine.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
best=ok.nsmallest(16,"chi2_planck").to_dict("records")
summary={"center":base,"n_total":int(len(df)),"n_stable":int(len(ok)),
         "best":best[0] if best else None,
         "best_delta_planck_vs_center":float(ok.chi2_planck.min()-BASE) if len(ok) else None,
         "target_note":"need roughly >0.46 combined P+D recovery over SN-close010 to close P+D and >=2 SN exactly"}
(OUT/"snclose010_structure_refine_summary.json").write_text(json.dumps(summary,indent=2))
print("SN010STRUCT_BEST",json.dumps(best,sort_keys=True),flush=True)
print("SN010STRUCT_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
