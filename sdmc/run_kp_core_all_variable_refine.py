#!/usr/bin/env python3
"""
Frozen-Kp-core all-variable refinement (terminal operator OFF).

Purpose
-------
Re-open the historical frozen Kp background/No-Slip core inside the present
hi_class likelihood pipeline, but do NOT activate the retired terminal C(eta)
operator. The historical values are reference coordinates, not hard priors.

Search variables (18D):
 H0, omega_b, omega_cdm, n_s, tau, ln10As,
 A_F, z_c, DeltaN, D0, p_D, D_floor,
 lambda_e, z_t, DeltaN_t, A_late, B_late, tau_B.

tau_A remains at the historical short-lobe scale 0.25 in this first broad
screen; all other active core coordinates are reopened. A follow-up local
screen may release tau_A too if a viable basin is found.

Hard gates
----------
  D>0, 0<c_s^2<=1 over the returned background grid, exact No-Slip closure.

Likelihood screen
-----------------
 Planck TTTEEE-lite + lowT + lowE + lensing + calibration prior.
The best stable candidates are promoted later to full Plik and raw DESI.
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

OUT=Path("output/kp_core_allvar"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Historical Kp reference quantities from Parts VI-VIII.
H0_KP=70.8514
OB_KP=0.02239952
OM_KP=0.1424000
OC_KP=OM_KP-OB_KP
AF_KP=0.02586
ZC_KP=4.38
W_KP=0.4467
LAM_KP=20.0
A_KP=0.024283
B_KP=0.0132767
SLATE_KP=0.874

# Present exact-covariant no-terminal leader, scored in the same screen.
LEADER=dict(H0=70.69065645344024,ob=0.022063253393818805,
            oc=0.1232890837695007,ns=0.9625227132590487,
            tau=0.055202901571989066,lnAs=3.0598696043919773,
            af=0.052875,zc=2.878125,width=0.54625,
            D0=0.34231919445927034,pD=1.0,df=0.03825,
            lam=18.278125,zt=16.9625,dNt=0.5,
            A=0.01105624999,B=0.01951933685,tauB=1.5)

# 18D broad all-variable box.
LOW=np.array([
  68.6, 0.02175, 0.1160, 0.953, 0.042, 3.015,
  0.015, 2.2, 0.42, 0.18, 0.25, 0.020,
  17.0, 15.5, 0.30, 0.004, 0.004, 0.90
],float)
HIGH=np.array([
  72.4, 0.02265, 0.1265, 0.976, 0.068, 3.085,
  0.080, 5.6, 1.55, 0.52, 1.50, 0.080,
  21.2, 22.0, 0.85, 0.042, 0.027, 2.20
],float)

NAMES=["H0","ob","oc","ns","tau","lnAs","af","zc","width","D0","pD","df",
       "lam","zt","dNt","A","B","tauB"]

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
    z=b["z"].to_numpy()
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
    return dict(z_eq_exact=zeq,z_star=zstar,D_M_star=DM,r_s_star=rs,
                ell_A=float(np.pi*DM/rs))

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
    Ap=float(op.x); p=pc(Ap)
    return dict(A_planck=Ap,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def ini(root,c):
    H0=float(c["H0"]); ob=float(c["ob"]); oc=float(c["oc"])
    h=H0/100.; ox=1.-(ob+oc+OR)/(h*h)
    As=math.exp(float(c["lnAs"]))/1e10
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
expansion_smg={ox:.17g},{float(c["lam"]):.17g},{float(c["zt"]):.17g},{float(c["dNt"]):.17g},{float(c["A"]):.17g},0.25,{float(c["B"]):.17g},{float(c["tauB"]):.17g}
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
    rec=dict(id=tag,**{k:float(c[k]) for k in NAMES},status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and thp.exists() and clp.exists():
        bg=table(bgp); th=table(thp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
                and rec["max_abs_noslip"]<1e-5)
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(derived(bg,th)); rec.update(pscore(clp)); rec["status"]="OK"
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    print("KPCORE_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[]

# Present no-terminal leader, exact same coordinates as the current Part IX leader.
leader=run("current_leader",LEADER); rows.append(leader)
if leader["status"]!="OK": raise RuntimeError("current leader failed same-pipeline screen")
BASE=float(leader["chi2_planck"])
print("KPCORE_CURRENT_LEADER",json.dumps(leader,sort_keys=True),flush=True)

# Historical Kp mapped into the present sdmc_full coordinates. zt=19.5 and
# dNt=0.5 are embedding coordinates, not claimed historical observables.
# Use the final Part-VIII late amplitude scaling s_late=0.874.
kp=dict(H0=H0_KP,ob=OB_KP,oc=OC_KP,ns=0.965,tau=0.0544,lnAs=math.log(2.10e-9*1e10),
        af=AF_KP,zc=ZC_KP,width=W_KP,D0=0.40,pD=1.0,df=1e-4,
        lam=LAM_KP,zt=19.5,dNt=0.5,A=A_KP*SLATE_KP,B=B_KP*SLATE_KP,tauB=1.5)
hist=run("historical_kp_mapped",kp); rows.append(hist)

# Also include the unscaled historical two-lobe target.
kp2=dict(kp); kp2["A"]=A_KP; kp2["B"]=B_KP
rows.append(run("historical_kp_unscaled",kp2))

# 64 deterministic Sobol points spanning every reopened core variable.
sam=qmc.Sobol(d=len(NAMES),scramble=False)
pts=qmc.scale(sam.random_base2(m=6),LOW,HIGH)
for i,p in enumerate(pts):
    c={k:float(v) for k,v in zip(NAMES,p)}
    rec=run(f"sobol{i:03d}",c)
    if rec.get("status")=="OK":
        rec["delta_planck_vs_current_leader"]=rec["chi2_planck"]-BASE
    rows.append(rec)

# add deltas to explicit references
for rec in rows[:3]:
    if rec.get("status")=="OK":
        rec["delta_planck_vs_current_leader"]=rec["chi2_planck"]-BASE

df=pd.DataFrame(rows)
df.to_csv(OUT/"kp_core_all_variable_refine.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
best=ok.nsmallest(12,"chi2_planck").to_dict("records")
print("KPCORE_BEST_PLANCK",json.dumps(best,sort_keys=True),flush=True)
summary={
  "current_leader_chi2_screen":BASE,
  "historical_kp_mapped":hist,
  "n_total":int(len(df)),
  "n_stable":int(len(ok)),
  "best":best[0] if best else None,
  "best_delta_vs_current_leader":float(ok.chi2_planck.min()-BASE) if len(ok) else None
}
(OUT/"kp_core_all_variable_summary.json").write_text(json.dumps(summary,indent=2))
print("KPCORE_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
