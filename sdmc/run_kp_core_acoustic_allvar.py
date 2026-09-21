#!/usr/bin/env python3
"""
Second-generation frozen-Kp-core all-variable refinement.

This closes the gap left by the two earlier Kp scans:
  * the broad 18-D scan kept tau_A fixed;
  * the historical-stability scan varied tau_A but froze ordinary cosmology.

Here all active Kp-core coordinates are reopened together with
(omega_b, omega_cdm, n_s, tau, Q), where
    Q = ln(10^10 A_s) - 2 tau.
H0 is not frozen: it is solved candidate-by-candidate to preserve the
current Part-IX acoustic scale ell_A. Thus the search does not throw away
most of its points merely by shifting the CMB peak positions.

The search is deliberately Kp-like: z_c and the late-lobe coordinates remain
in the historical/stabilized Kp basin rather than collapsing immediately back
onto the present geometry-refined leader.
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

OUT=Path("output/kp_core_acoustic_allvar"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Current fine-Q geometry-refined point: screen reference only.
REF=dict(
 H0=70.79608815124146,ob=0.022083219194622913,oc=0.12299536722293603,
 ns=0.9625227132590487,tau=0.055202901571989066,Q=2.944463801247999,
 af=0.04984375,zc=2.98625,width=0.515,D0=0.34231919445927034,pD=1.0,df=0.041875,
 lam=18.40625,zt=16.496875,dNt=0.5,A=0.01105624999,tauA=0.25,
 B=0.01951933685,tauB=1.5
)

# Best non-reference stable Kp-like seed from the historical-stability scan.
KPSEED=dict(
 H0=70.8514,ob=0.02239952,oc=0.12000048,ns=0.965,tau=0.0544,
 Q=3.044522437723423-2*0.0544,
 af=0.02540625,zc=5.3203125,width=0.420625,D0=0.20875,pD=1.1625,df=0.09628125,
 lam=19.04375,zt=19.071875,dNt=0.701875,A=0.016125,tauA=0.5965625,
 B=0.0089375,tauB=1.71875
)

# All free coordinates except H0, which is acoustically solved.
NAMES=["ob","oc","ns","tau","Q","af","zc","width","D0","pD","df",
       "lam","zt","dNt","A","tauA","B","tauB"]
LOW=np.array([
 0.02190,0.1185,0.9580,0.0480,2.925,
 0.0180,4.20,0.30,0.14,0.75,0.045,
 18.0,17.2,0.40,0.008,0.18,0.0045,0.90
],float)
HIGH=np.array([
 0.02250,0.1245,0.9700,0.0620,2.965,
 0.0400,5.80,0.68,0.38,1.45,0.130,
 20.6,20.8,0.88,0.027,0.78,0.0170,2.35
],float)

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
    imax=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,z,b["comov. dist."]))
    rs=float(np.interp(zstar,z,b["comov.snd.hrz."]))
    matter=b["(.)rho_b"].to_numpy()+b["(.)rho_cdm"].to_numpy()
    rad=b["(.)rho_g"].to_numpy()+b["(.)rho_ur"].to_numpy()
    ratio=matter/rad; q=ratio-1.; ii=np.where(q[:-1]*q[1:]<=0)[0]
    zeq=float("nan")
    if len(ii):
        i=ii[-1]; zeq=float(z[i]+(1-ratio[i])*(z[i+1]-z[i])/(ratio[i+1]-ratio[i]))
    return dict(z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=float(np.pi*DM/rs),z_eq_exact=zeq)

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

def ini(root,c,full):
    H0=float(c["H0"]); ob=float(c["ob"]); oc=float(c["oc"]); h=H0/100.
    ox=1.-(ob+oc+OR)/(h*h)
    lnAs=float(c["Q"])+2.*float(c["tau"]); As=math.exp(lnAs)/1e10
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

def one(tag,c,full=False):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,c,full))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat")
    if cp.returncode==0 and bgp.exists() and thp.exists():
        bg=table(bgp); th=table(thp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1 and rec["max_abs_noslip"]<1e-5
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(derived(bg,th)); rec["status"]="OK"
            if full:
                rec.update(pscore(root+"00_cl_lensed.dat"))
    if rec["status"]!="OK": rec["error"]=cp.stdout[-700:].replace("\n"," | ")
    return rec

# Current reference and target acoustic scale.
refc=dict(REF); ref=one("reference",refc,True)
if ref["status"]!="OK": raise RuntimeError("reference failed")
TARGET=float(ref["ell_A"]); BASE=float(ref["chi2_planck"])
print("KPALL_REFERENCE",json.dumps({**ref,**REF},sort_keys=True),flush=True)

def acoustic_candidate(tag,c0):
    c=dict(c0)
    vals=[]
    for j,H in enumerate([68.8,70.0,71.2,72.4,73.6]):
        cc=dict(c); cc["H0"]=H
        r=one(f"{tag}_root{j}",cc,False)
        if r.get("status")=="OK":
            vals.append((H,float(r["ell_A"]-TARGET)))
    vals=sorted(vals)
    pair=None
    for a,b in zip(vals[:-1],vals[1:]):
        if a[1]==0 or a[1]*b[1] <=0:
            pair=(a,b); break
    if pair is None:
        rec=dict(id=tag,status="NO_STABLE_ACOUSTIC_ROOT",n_root=len(vals))
        print("KPALL_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
    (h1,g1),(h2,g2)=pair
    H=float(h1+(0.-g1)*(h2-h1)/(g2-g1))
    c["H0"]=H
    r=one(tag,c,True)
    rec=dict(id=tag,**{k:float(c[k]) for k in ["H0"]+NAMES},**r)
    rec["lnAs"]=float(c["Q"])+2.*float(c["tau"])
    if rec.get("status")=="OK":
        rec["delta_planck_vs_qfine"]=rec["chi2_planck"]-BASE
        rec["delta_ellA"]=rec["ell_A"]-TARGET
    print("KPALL_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[]
rows.append(acoustic_candidate("kp_seed",KPSEED))

sam=qmc.Sobol(d=len(NAMES),scramble=False)
pts=qmc.scale(sam.random_base2(m=5),LOW,HIGH)  # 32 all-variable Kp-like points
for i,p in enumerate(pts):
    c={k:float(v) for k,v in zip(NAMES,p)}
    rows.append(acoustic_candidate(f"sobol{i:03d}",c))

df=pd.DataFrame(rows); df.to_csv(OUT/"kp_core_acoustic_allvar.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
best=ok.nsmallest(12,"chi2_planck").to_dict("records") if len(ok) else []
summary={
 "reference_chi2":BASE,"target_ellA":TARGET,"n_total":int(len(df)),"n_stable":int(len(ok)),
 "best":best[0] if best else None,
 "best_delta_vs_qfine":float(ok.chi2_planck.min()-BASE) if len(ok) else None
}
(OUT/"kp_core_acoustic_allvar_summary.json").write_text(json.dumps(summary,indent=2))
print("KPALL_BEST",json.dumps(best,sort_keys=True),flush=True)
print("KPALL_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
