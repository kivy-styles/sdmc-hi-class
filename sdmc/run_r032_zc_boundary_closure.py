#!/usr/bin/env python3
"""
Boundary continuation in z_c around the final r032 zt-best + Qfine-best point.

Purpose:
  The final joint z_t closure winner landed at z_c=2.98625 while that
  search box ended at z_c=3.00.  This test extends that direction rather
  than assuming it is enclosed.

At each z_c:
  * keep the final matter densities and SDMC structural parameters fixed,
  * keep (n_s,tau,Q) at the Qfine-best solution,
  * re-solve H0 to the established acoustic target ell_A,
  * enforce D>0, c_s^2>0, c_s^2<=1,
  * score Planck TTTEEE-lite + lowT + lowE + lensing + calibration prior.

Any improved point is only a screen winner and must later be promoted to
full Plik + raw DESI + exact-covariant replay.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/r032_zc_boundary_closure"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5
TARGET=301.6798271349195

H00=70.79608815124146
OB=0.022083219194622913
OC=0.12299536722293603
NS=0.9625227132590487
TAU=0.055202901571989066
Q=2.944463801247999
LNAS=Q+2.*TAU
AS=math.exp(LNAS)/1e10

AF=0.04984375
ZC0=2.98625
WIDTH=0.515
D0=0.34231919445927034
DF=0.041875
LAM=18.40625
ZT=16.496875

ZCS=[2.86,2.90,2.94,2.96,2.98,2.98625,3.00,3.02,3.04,3.06,3.08,3.10,3.125,3.15,3.20,3.25,3.30]

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

def derived(bg,th):
    b=bg.sort_values("z")
    imax=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,b["z"],b["comov. dist."]))
    rs=float(np.interp(zstar,b["z"],b["comov.snd.hrz."]))
    return dict(z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=float(np.pi*DM/rs))

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

def ini(root,H0,zc,full):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    obs="""modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
""" if full else ""
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
parameters_smg={AF:.17g},{zc:.17g},{WIDTH:.17g},{D0:.17g},1.0,{DF:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM:.17g},{ZT:.17g},0.5,0.01105624999,0.25,0.01951933685,1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
output_background_smg=3
{obs}
write background=yes
write thermodynamics=yes
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

def run_class(tag,H0,zc,full):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini")
    ip.write_text(ini(root,H0,zc,full))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=300)
    if cp.returncode:
        raise RuntimeError(cp.stdout[-1200:])
    bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat")
    der=derived(bg,th)
    stab=dict(min_D=float(bg["kin (D)"].min()),
              min_cs2=float(bg["c_s^2"].min()),
              max_cs2=float(bg["c_s^2"].max()),
              max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    return root,bg,der,stab

def candidate(tag,zc):
    vals=[]
    for H in np.linspace(H00-0.35,H00+0.35,8):
        try:
            _,_,d,st=run_class(f"{tag}_root_{H:.5f}".replace(".","p"),float(H),zc,False)
            stable=st["min_D"]>0 and st["min_cs2"]>0 and st["max_cs2"]<=1.
            if stable: vals.append((float(H),float(d["ell_A"]-TARGET)))
        except Exception:
            pass
    vals=sorted(vals)
    pair=None
    for a,b in zip(vals[:-1],vals[1:]):
        if a[1]==0 or a[1]*b[1] <= 0:
            pair=(a,b); break
    if pair is None:
        if not vals:
            rec=dict(id=tag,z_c=zc,status="NO_STABLE_ACOUSTIC_POINTS")
            print("ZCBC_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
        Hbest,gbest=min(vals,key=lambda q:abs(q[1]))
        if abs(gbest)>0.04:
            rec=dict(id=tag,z_c=zc,status="NO_ACOUSTIC_ROOT",
                     closest_H0=Hbest,closest_delta_ellA=gbest)
            print("ZCBC_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
        H0=Hbest
    else:
        (h1,g1),(h2,g2)=pair
        H0=float(h1+(0.-g1)*(h2-h1)/(g2-g1))
    try:
        root,bg,der,st=run_class(tag,H0,zc,True)
    except Exception as e:
        rec=dict(id=tag,z_c=zc,H0=H0,status="FINAL_CLASS_FAIL",error=str(e)[-600:])
        print("ZCBC_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
    stable=st["min_D"]>0 and st["min_cs2"]>0 and st["max_cs2"]<=1.
    rec=dict(id=tag,z_c=float(zc),H0=H0,delta_ellA=der["ell_A"]-TARGET,
             stable_subluminal=bool(stable),**der,**st,status="OK" if stable else "UNSTABLE")
    if stable: rec.update(pscore(root+"00_cl_lensed.dat"))
    print("ZCBC_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

base=candidate("current_zc",ZC0)
if base.get("status")!="OK": raise RuntimeError("current z_c baseline failed")
BASE=float(base["chi2_planck"])
rows=[base]
for i,zc in enumerate(ZCS):
    if abs(zc-ZC0)<1e-10: continue
    rows.append(candidate(f"zc_{i:02d}",float(zc)))

for r in rows:
    if r.get("status")=="OK":
        r["delta_planck_vs_current"]=r["chi2_planck"]-BASE
        r["delta_high_vs_current"]=r["chi2_high"]-base["chi2_high"]
        r["delta_lowT_vs_current"]=r["chi2_lowT"]-base["chi2_lowT"]
        r["delta_lowE_vs_current"]=r["chi2_lowE"]-base["chi2_lowE"]
        r["delta_lensing_vs_current"]=r["chi2_lensing"]-base["chi2_lensing"]
        r["delta_cal_vs_current"]=r["chi2_cal"]-base["chi2_cal"]

df=pd.DataFrame(rows)
df.to_csv(OUT/"r032_zc_boundary_closure.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
best=ok.nsmallest(10,"chi2_planck").to_dict("records")
summary={
  "target_ellA":TARGET,
  "baseline":base,
  "best":best[0] if best else None,
  "best_delta_vs_current":float(ok.chi2_planck.min()-BASE) if len(ok) else None,
  "max_tested_zc":max(ZCS),
  "best_at_upper_edge":bool(len(ok) and abs(float(ok.nsmallest(1,"chi2_planck").iloc[0].z_c)-max(ZCS))<1e-9)
}
(OUT/"r032_zc_boundary_closure_summary.json").write_text(json.dumps(summary,indent=2))
print("ZCBC_BEST",json.dumps(best,sort_keys=True),flush=True)
print("ZCBC_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
