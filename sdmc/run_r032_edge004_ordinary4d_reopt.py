#!/usr/bin/env python3
"""
Acoustically locked ordinary-parameter re-optimization around edge004.

The structural sector is fixed at edge004.  Re-open:
  omega_b, omega_cdm, n_s, Q = ln(10^10 A_s)-2 tau
with tau fixed.  H0 is re-solved at every point to the established acoustic
scale ell_A.  This tests whether the large structural shift moved the optimum
of the ordinary cosmological parameters.

This is a candidate screen using Planck TTTEEE-lite + lowT + lowE + lensing
and the calibration prior.  Any useful winner must be promoted to exact
covariant full Plik + raw DESI.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/r032_edge004_ordinary4d_reopt"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5
TARGET=301.6798271349195

H00=70.78157258819029
OB0=0.022083219194622913
OC0=0.12299536722293603
NS0=0.9625227132590487
TAU=0.055202901571989066
Q0=2.941463801247999

AF=0.034172075465787204
ZC=3.943745281267911
WIDTH=0.36045085035264496
D0=0.34231919445927034
DF=0.0479059932352975
LAM=18.40625
ZT=15.418161678034812

# omega_b, omega_cdm, n_s, Q
LOW=np.array([0.02190,0.1216,0.9600,2.9360],float)
HIGH=np.array([0.02228,0.1244,0.9660,2.9445],float)

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

def ini(root,H0,ob,oc,ns,Q,full):
    As=math.exp(Q+2.*TAU)/1e10
    h=H0/100.; ox=1.-(ob+oc+OR)/(h*h)
    obs="""modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
""" if full else ""
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={ob:.17g}
omega_cdm={oc:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={As:.17e}
n_s={ns:.17g}
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

def run_class(tag,H0,ob,oc,ns,Q,full):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini")
    ip.write_text(ini(root,H0,ob,oc,ns,Q,full))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=300)
    if cp.returncode: raise RuntimeError(cp.stdout[-1200:])
    bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat")
    der=derived(bg,th)
    st=dict(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
            max_cs2=float(bg["c_s^2"].max()),
            max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    return root,der,st

def candidate(tag,ob,oc,ns,Q):
    vals=[]
    for H in np.linspace(H00-1.2,H00+1.2,13):
        try:
            _,d,st=run_class(f"{tag}_root_{H:.5f}".replace(".","p"),float(H),ob,oc,ns,Q,False)
            if st["min_D"]>0 and st["min_cs2"]>0 and st["max_cs2"]<=1.:
                vals.append((float(H),float(d["ell_A"]-TARGET)))
        except Exception: pass
    vals=sorted(vals); pair=None
    for a,b in zip(vals[:-1],vals[1:]):
        if a[1]==0 or a[1]*b[1]<=0: pair=(a,b); break
    if pair is None:
        rec=dict(id=tag,omega_b=ob,omega_cdm=oc,n_s=ns,Q=Q,status="NO_ACOUSTIC_ROOT")
        print("O4D_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
    (h1,g1),(h2,g2)=pair
    H0=float(h1+(0.-g1)*(h2-h1)/(g2-g1))
    try:
        root,der,st=run_class(tag,H0,ob,oc,ns,Q,True)
    except Exception as e:
        rec=dict(id=tag,omega_b=ob,omega_cdm=oc,n_s=ns,Q=Q,H0=H0,status="FINAL_FAIL",error=str(e)[-500:])
        print("O4D_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
    stable=st["min_D"]>0 and st["min_cs2"]>0 and st["max_cs2"]<=1.
    rec=dict(id=tag,omega_b=float(ob),omega_cdm=float(oc),n_s=float(ns),Q=float(Q),
             A_s=float(math.exp(Q+2*TAU)/1e10),H0=H0,stable_subluminal=bool(stable),
             delta_ellA=der["ell_A"]-TARGET,**der,**st,status="OK" if stable else "UNSTABLE")
    if stable: rec.update(pscore(root+"00_cl_lensed.dat"))
    print("O4D_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec

rows=[]
base=candidate("center",OB0,OC0,NS0,Q0)
if base.get("status")!="OK": raise RuntimeError("center failed")
BASE=float(base["chi2_planck"]); rows.append(base)

anchors=[
 ("obm",OB0-0.00010,OC0,NS0,Q0),("obp",OB0+0.00010,OC0,NS0,Q0),
 ("ocm",OB0,OC0-0.0007,NS0,Q0),("ocp",OB0,OC0+0.0007,NS0,Q0),
 ("nsm",OB0,OC0,NS0-0.0012,Q0),("nsp",OB0,OC0,NS0+0.0012,Q0),
 ("qm",OB0,OC0,NS0,Q0-0.0010),("qp",OB0,OC0,NS0,Q0+0.0010),
]
for a in anchors: rows.append(candidate(*a))

sam=qmc.Sobol(d=4,scramble=True,seed=20260927)
pts=qmc.scale(sam.random_base2(m=5),LOW,HIGH)
for i,p in enumerate(pts):
    rows.append(candidate(f"sobol{i:03d}",*map(float,p)))

for r in rows:
    if r.get("status")=="OK":
        r["delta_planck_vs_center"]=r["chi2_planck"]-BASE
        r["delta_high_vs_center"]=r["chi2_high"]-base["chi2_high"]
        r["delta_lowT_vs_center"]=r["chi2_lowT"]-base["chi2_lowT"]
        r["delta_lowE_vs_center"]=r["chi2_lowE"]-base["chi2_lowE"]
        r["delta_lensing_vs_center"]=r["chi2_lensing"]-base["chi2_lensing"]
        r["delta_cal_vs_center"]=r["chi2_cal"]-base["chi2_cal"]

df=pd.DataFrame(rows); df.to_csv(OUT/"r032_edge004_ordinary4d_reopt.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
best=ok.nsmallest(12,"chi2_planck").to_dict("records")
summary={"baseline":base,"n_total":int(len(df)),"n_stable":int(len(ok)),
         "best":best[0] if best else None,
         "best_delta_vs_center":float(ok.chi2_planck.min()-BASE) if len(ok) else None}
(OUT/"r032_edge004_ordinary4d_reopt_summary.json").write_text(json.dumps(summary,indent=2))
print("O4D_BEST",json.dumps(best,sort_keys=True),flush=True)
print("O4D_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
