#!/usr/bin/env python3
from pathlib import Path
import json,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/bridge6796_earlytracker_grid"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5
H0=67.96
OB=0.022083219194622913; OC=0.12299536722293603
AS=2.119721074504234e-9; NS=0.9625227132590487; TAU=0.055202901571989066
AF=0.024290704212870225; ZC=3.6096407580714724; WIDTH=0.3631213944496205
D0=0.34231919445927034; DF=0.04602317962943721
AL=0.011061627057101114; BL=0.010557478912554592
TAUA=.25; TAUB=1.5; DN=.5
h=H0/100.; OX=1.-(OB+OC+OR)/(h*h)

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
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
    a=np.loadtxt(path); e=a[:,0].astype(int); n=int(e.max())+1; cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[e]=a[:,1]*cv; ee[e]=a[:,2]*cv; te[e]=a[:,3]*cv
    L=e.astype(float); pp[e]=a[:,5]*L*(L+1)
    return np.column_stack([e,tt[e],te[e],ee[e]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    ha,d=load_cls(path)
    def pieces(A):
        x=float(high.chi_squared(ha,A_planck=A))
        x+=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        x+=float(-2*lens.log_likelihood(d,**lp))
        x+=((A-1)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(pieces,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def ini(root,lam,zt):
    return textwrap.dedent(f"""\
H0={H0}
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
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={OX:.17g},{lam:.17g},{zt:.17g},{DN},{AL},{TAUA},{BL},{TAUB}
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

rows=[]
lams=[7.0,8.0,9.0,10.0,11.0,12.0,13.5,15.0,16.5,18.40625,20.5,23.0,26.0]
zts=[8.0,12.0,16.173189924377947,22.0,30.0,45.0]
for lam in lams:
  for zt in zts:
    tag=f"l{lam:g}_z{zt:g}".replace(".","p")
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,lam,zt))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    r={"lambda_e":lam,"z_t":zt,"status":"FAIL","returncode":cp.returncode}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
      try:
        bg=table(bgp)
        r["min_D"]=float(bg["kin (D)"].min()); r["min_cs2"]=float(bg["c_s^2"].min()); r["max_cs2"]=float(bg["c_s^2"].max())
        r["max_abs_noslip"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1 and r["max_abs_noslip"]<1e-4
        r["stable_subluminal"]=bool(stable)
        if stable:
          chi,Apl=pscore(clp); r.update(status="OK",chi2_planck_lite=chi,A_planck=Apl)
      except Exception as e: r["error"]=str(e)
    if r["status"]!="OK" and "error" not in r: r["error"]=cp.stdout[-900:].replace("\n"," | ")
    rows.append(r); print("EARLYTRACKER_GRID_POINT",json.dumps(r,sort_keys=True),flush=True)
    for p in OUT.glob(tag+"_*"):
      try:p.unlink()
      except:pass
df=pd.DataFrame(rows)
df.to_csv(OUT/"bridge6796_earlytracker_grid.csv",index=False)
ok=df[df.status.eq("OK")].sort_values("chi2_planck_lite")
best=ok.head(20)
best.to_csv(OUT/"bridge6796_earlytracker_top20.csv",index=False)
summary={"n_total":len(df),"n_ok":len(ok),"best":best.to_dict("records")}
(OUT/"bridge6796_earlytracker_summary.json").write_text(json.dumps(summary,indent=2))
print("EARLYTRACKER_GRID_BEST",json.dumps(summary,sort_keys=True),flush=True)
