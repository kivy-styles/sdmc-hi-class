#!/usr/bin/env python3
"""
Structural Planck-lite refinement at the SN-friendly sobol223 background.

The late expansion geometry is frozen at the SN-improved sobol223 point:
  H0 = 69.85635133907199
  A_late = 0.018589897081255913
  B_late = 0.013815921754576268

Only the No-Slip structural coordinates (A_F,z_c,width,D_floor,z_t) vary.
Thus the low-z SN distance curve is effectively fixed; the purpose is to
recover CMB/growth likelihood without sacrificing the SN gain.
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

OUT=Path("output/snclosure223_structural_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0=69.85635133907199
OB=0.022083219194622913
OC=0.12299536722293603
NS=0.9625227132590487
TAU=0.055202901571989066
Q=2.942463801247999
AS=math.exp(Q+2*TAU)/1e10

A_LATE=0.018589897081255913
B_LATE=0.013815921754576268
LAM=18.40625

AF0=0.03200255395658314
ZC0=4.007449422683567
W0=0.34799921004101636
D0=0.34231919445927034
DF0=0.04607192634791136
ZT0=16.173189924377947

LOW=np.array([0.0250,3.70,0.305,0.0420,15.90],float)
HIGH=np.array([0.0400,4.35,0.405,0.0560,16.45],float)

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
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,af,zc,w,df,zt):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
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
expansion_smg={ox:.17g},{LAM:.17g},{zt:.17g},0.5,{A_LATE:.17g},0.25,{B_LATE:.17g},1.5
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

def run(label,af,zc,w,df,zt):
    root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,af,zc,w,df,zt))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=label,A_F=float(af),z_c=float(zc),width=float(w),D_floor=float(df),z_t=float(zt),
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
    if rec["status"]!="OK": rec["error"]=cp.stdout[-800:].replace("\n"," | ")
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    print("SN223STRUCT_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[run("center",AF0,ZC0,W0,DF0,ZT0)]
if rows[0]["status"]!="OK": raise RuntimeError(rows[0])
C0=rows[0]["chi2_planck"]

# axis probes
for name,v in [
 ("afm",(AF0-0.003,ZC0,W0,DF0,ZT0)),("afp",(AF0+0.003,ZC0,W0,DF0,ZT0)),
 ("zcm",(AF0,ZC0-0.18,W0,DF0,ZT0)),("zcp",(AF0,ZC0+0.18,W0,DF0,ZT0)),
 ("wm",(AF0,ZC0,W0-0.025,DF0,ZT0)),("wp",(AF0,ZC0,W0+0.025,DF0,ZT0)),
 ("dfm",(AF0,ZC0,W0,DF0-0.003,ZT0)),("dfp",(AF0,ZC0,W0,DF0+0.004,ZT0)),
 ("ztm",(AF0,ZC0,W0,DF0,ZT0-0.18)),("ztp",(AF0,ZC0,W0,DF0,ZT0+0.18))]:
    rows.append(run(name,*v))

sob=qmc.Sobol(d=5,scramble=True,seed=223501)
pts=qmc.scale(sob.random_base2(m=6),LOW,HIGH)
for i,p in enumerate(pts):
    rows.append(run(f"s{i:03d}",*map(float,p)))

df=pd.DataFrame(rows)
df["delta_vs_center"]=df.chi2_planck-C0
df.to_csv(OUT/"snclosure223_structural_refine.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].sort_values("chi2_planck")
best=ok.head(12).to_dict("records")
summary={"center_chi2":C0,"n_total":int(len(df)),"n_stable":int(len(ok)),
         "best":best[0],"best_delta_vs_center":float(best[0]["chi2_planck"]-C0),
         "top12":best}
(OUT/"snclosure223_structural_refine_summary.json").write_text(json.dumps(summary,indent=2))
print("SN223STRUCT_BEST",json.dumps(best,sort_keys=True),flush=True)
print("SN223STRUCT_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
