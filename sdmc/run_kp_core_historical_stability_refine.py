#!/usr/bin/env python3
"""
Historical frozen-Kp core stabilization/refinement.

Stage 1 of the requested Kp campaign:
- Keep the historical Kp ordinary cosmology fixed.
- Reopen the full active Kp core structure, including tau_A which the earlier
  18D broad screen kept fixed.
- Find a stable/sub-luminal Kp-like seed before reopening ordinary cosmology.

Variables:
 A_F, z_c, DeltaN, D0, p_D, D_floor,
 lambda_e, z_t, DeltaN_t, A_late, tau_A, B_late, tau_B.

Reference:
 current exact-covariant Part-IX zt-best point is scored in the same screen.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/kp_core_historical_stability"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Historical frozen-Kp ordinary cosmology.
ORD=dict(H0=70.8514,ob=0.02239952,oc=0.12000048,ns=0.965,
         tau=0.0544,lnAs=math.log(2.10e-9*1e10))

# Historical core coordinates.
KP=dict(af=0.02586,zc=4.38,width=0.4467,D0=0.40,pD=1.0,df=1e-4,
        lam=20.0,zt=19.5,dNt=0.5,A=0.021223342,tauA=0.25,
        B=0.0116038358,tauB=1.5)

# Current exact-covariant Part-IX leader (zt-best), for same-pipeline reference.
LEADER=dict(H0=70.79608815124146,ob=0.022083219194622913,
            oc=0.12299536722293603,ns=0.9625227132590487,
            tau=0.055202901571989066,lnAs=3.0598696043919773,
            af=0.04984375,zc=2.98625,width=0.515,
            D0=0.34231919445927034,pD=1.0,df=0.041875,
            lam=18.40625,zt=16.496875,dNt=0.5,
            A=0.01105624999,tauA=0.25,B=0.01951933685,tauB=1.5)

NAMES=["af","zc","width","D0","pD","df","lam","zt","dNt","A","tauA","B","tauB"]
LOW=np.array([0.012,3.10,0.28,0.16,0.45,0.012,18.3,17.0,0.28,0.010,0.08,0.005,0.65],float)
HIGH=np.array([0.045,5.55,0.78,0.68,1.65,0.105,21.7,22.1,0.82,0.038,0.65,0.023,2.45],float)

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
    b=bg.sort_values("z"); z=b["z"].to_numpy()
    matter=b["(.)rho_b"].to_numpy()+b["(.)rho_cdm"].to_numpy()
    rad=b["(.)rho_g"].to_numpy()+b["(.)rho_ur"].to_numpy()
    ratio=matter/rad; q=ratio-1.; ii=np.where(q[:-1]*q[1:]<=0)[0]
    zeq=float("nan")
    if len(ii):
        i=ii[-1]; zeq=float(z[i]+(1-ratio[i])*(z[i+1]-z[i])/(ratio[i+1]-ratio[i]))
    imax=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,z,b["comov. dist."]))
    rs=float(np.interp(zstar,z,b["comov.snd.hrz."]))
    return dict(z_eq_exact=zeq,z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=float(np.pi*DM/rs))

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    harr,dls=load_cls(path)
    def pc(Ap):
        ch=float(high.chi_squared(harr,A_planck=Ap))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=Ap))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=Ap))
        lp={lens.calibration_param:Ap} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((Ap-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda x:pc(float(x))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,c):
    H0=float(c["H0"]); ob=float(c["ob"]); oc=float(c["oc"]); h=H0/100.
    ox=1.-(ob+oc+OR)/(h*h); As=math.exp(float(c["lnAs"]))/1e10
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={ob:.17g}
omega_cdm={oc:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={As:.17e}
n_s={float(c["ns"]):.17g}
tau_reio={float(c["tau"]):.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={float(c["af"]):.17g},{float(c["zc"]):.17g},{float(c["width"]):.17g},{float(c["D0"]):.17g},{float(c["pD"]):.17g},{float(c["df"]):.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{float(c["lam"]):.17g},{float(c["zt"]):.17g},{float(c["dNt"]):.17g},{float(c["A"]):.17g},{float(c["tauA"]):.17g},{float(c["B"]):.17g},{float(c["tauB"]):.17g}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
output_background_smg=3
modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
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

def run(tag,c):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=tag,**{k:float(c[k]) for k in ["H0","ob","oc","ns","tau","lnAs"]+NAMES},
             status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and thp.exists() and clp.exists():
        bg=table(bgp); th=table(thp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1 and rec["max_abs_noslip"]<1e-5)
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(derived(bg,th)); rec.update(pscore(clp)); rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    print("KPHIST_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[]
lead=run("partix_ztbest_reference",LEADER); rows.append(lead)
if lead["status"]!="OK": raise RuntimeError("Part IX reference failed")
BASE=float(lead["chi2_planck"])

# Historical exact mapped point.
hist=dict(ORD,**KP)
rows.append(run("historical_kp",hist))

# Explicit stability ladder around the historical structural point.
for j,df in enumerate([0.015,0.025,0.035,0.045,0.060,0.080,0.100]):
    c=dict(hist); c["df"]=df
    rows.append(run(f"historical_floor_{j:02d}",c))

# 128 deterministic Sobol structural points with historical ordinary cosmology.
sam=qmc.Sobol(d=len(NAMES),scramble=False)
pts=qmc.scale(sam.random_base2(m=7),LOW,HIGH)
for i,p in enumerate(pts):
    c=dict(ORD); c.update({k:float(v) for k,v in zip(NAMES,p)})
    rows.append(run(f"sobol{i:03d}",c))

for rec in rows:
    if rec.get("status")=="OK":
        rec["delta_planck_vs_partix_ztbest"]=rec["chi2_planck"]-BASE

df=pd.DataFrame(rows); df.to_csv(OUT/"kp_core_historical_stability_refine.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
best=ok.nsmallest(15,"chi2_planck").to_dict("records")
print("KPHIST_BEST_PLANCK",json.dumps(best,sort_keys=True),flush=True)
summary={"partix_reference_chi2":BASE,"n_total":int(len(df)),"n_stable":int(len(ok)),
         "best":best[0] if best else None,
         "best_delta_vs_partix_ztbest":float(ok.chi2_planck.min()-BASE) if len(ok) else None}
(OUT/"kp_core_historical_stability_summary.json").write_text(json.dumps(summary,indent=2))
print("KPHIST_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
