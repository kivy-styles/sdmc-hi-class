#!/usr/bin/env python3
"""
Retired Kp terminal physical-operator refinement.

The physical reciprocal terminal operator has four runtime coordinates:
  Cmin, x1, x2, width
where x is the normalized star-to-drag residual acoustic coordinate.

Stage T1:
  - keep the historical frozen Kp tracker/core cosmology fixed;
  - compare the exact null (operator disabled) with 64 deterministic Sobol
    terminal profiles plus manuscript-motivated seeds;
  - score Planck TTTEEE-lite + lowT + lowE + lensing + calibration prior;
  - record the full spectrum-shape response.

A successful terminal profile is only a seed. It must subsequently be combined
with the stable Kp core / ordinary cosmology and promoted through full Plik,
DESI and exact-covariant gates.
"""
from pathlib import Path
import json, math, re, subprocess
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/kp_terminal_physical_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025
CONFIG=Path("sdmc/config/kp_tracker_clustered_cs0003.ini")
PARTIX_SCREEN=1028.6314402250505

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp},a

def pscore(path):
    harr,dls,a=load_cls(path)
    def pc(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.96,1.04),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5]),a

def make_ini(tag,enabled,Cmin,x1,x2,width):
    txt=CONFIG.read_text()
    txt=re.sub(r"(?m)^\s*root\s*=.*$",f"root = {OUT/(tag+'_')}",txt)
    txt += f"""
kp_terminal_enabled = {1 if enabled else 0}
kp_terminal_Cmin = {Cmin:.17g}
kp_terminal_x1 = {x1:.17g}
kp_terminal_x2 = {x2:.17g}
kp_terminal_width = {width:.17g}
"""
    ip=OUT/(tag+".ini"); ip.write_text(txt); return ip

def run(tag,enabled,Cmin,x1,x2,width,base_spec=None):
    ip=make_ini(tag,enabled,Cmin,x1,x2,width)
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=tag,enabled=bool(enabled),Cmin=Cmin,x1=x1,x2=x2,width=width,
             status="FAIL",returncode=cp.returncode)
    clp=OUT/(tag+"_00_cl_lensed.dat")
    if cp.returncode==0 and clp.exists():
        try:
            sc,a=pscore(clp); rec.update(sc); rec["status"]="OK"
            rec["delta_vs_partix_screen"]=rec["chi2_planck"]-PARTIX_SCREEN
            if base_spec is not None:
                # raw D_ell ratios at representative multipoles
                ell=a[:,0]
                for L in (200,500,1000,1500,2000,2500):
                    i=int(np.argmin(np.abs(ell-L)))
                    ib=int(np.argmin(np.abs(base_spec[:,0]-L)))
                    for j,nm in ((1,"TT"),(2,"TE"),(3,"EE")):
                        den=float(base_spec[ib,j])
                        rec[f"{nm}_ratio_{L}"]=float(a[i,j]/den) if den!=0 else float("nan")
        except Exception as e:
            rec["error"]=str(e)
    if rec["status"]!="OK" and "error" not in rec:
        rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    print("KPTERM_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

# Exact null.
base=run("terminal_null",False,1.0,0.0207,0.9793,0.012,None)
if base["status"]!="OK": raise RuntimeError("Kp terminal null failed")
_,base_spec=pscore(OUT/"terminal_null_00_cl_lensed.dat")
BASE=base["chi2_planck"]
rows=[base]

# Manuscript / prior-inspired explicit seeds.
seeds=[
 ("seed_c054",0.539812,0.0207,0.9793,0.012),
 ("seed_wide",0.539812,0.0207,0.9793,0.025),
 ("seed_soft",0.70,0.0207,0.9793,0.012),
 ("seed_late",0.539812,0.05,0.96,0.012),
]
for name,c,x1,x2,w in seeds:
    rows.append(run(name,True,c,x1,x2,w,base_spec))

# [Cmin,x1,x2,width]
LOW=np.array([0.35,0.002,0.82,0.003],float)
HIGH=np.array([0.98,0.18,0.998,0.060],float)
sam=qmc.Sobol(d=4,scramble=False)
pts=qmc.scale(sam.random_base2(m=6),LOW,HIGH)
for i,p in enumerate(pts):
    c,x1,x2,w=map(float,p)
    if x1>=x2:
        continue
    rows.append(run(f"sobol{i:03d}",True,c,x1,x2,w,base_spec))

for rec in rows:
    if rec.get("status")=="OK":
        rec["delta_vs_terminal_null"]=rec["chi2_planck"]-BASE

df=pd.DataFrame(rows)
df.to_csv(OUT/"kp_terminal_physical_refine.csv",index=False)
ok=df[df.status=="OK"].copy()
best=ok.nsmallest(12,"chi2_planck").to_dict("records")
summary={
 "terminal_null_chi2":float(BASE),
 "partix_reference_screen":PARTIX_SCREEN,
 "n_ok":int(len(ok)),
 "best":best[0] if best else None,
 "best_delta_vs_terminal_null":float(ok.chi2_planck.min()-BASE) if len(ok) else None,
 "best_delta_vs_partix_screen":float(ok.chi2_planck.min()-PARTIX_SCREEN) if len(ok) else None
}
(OUT/"kp_terminal_physical_summary.json").write_text(json.dumps(summary,indent=2))
print("KPTERM_BEST_PLANCK",json.dumps(best,sort_keys=True),flush=True)
print("KPTERM_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
