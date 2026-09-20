#!/usr/bin/env python3
"""
Focused four-dimensional Planck screen:
    log10(c_phi^2), ln(1e10 A_s), tau_reio, n_s
at the frozen exact-leader SDMC acoustic geometry.

Motivation:
  * exact-leader full Plik pulls A_planck to ~1.00825;
  * lowering the previously frozen tracker cs2_fld from 0.003 to 1e-4
    improves Planck-lite+lensing by ~4.74;
  * H0/omega_b/omega_c moves away from the exact leader are very costly.

This is a screening calculation, not a final model-comparison result.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/coldtracker_primordial")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255
SIG=0.0025

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
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path)
    ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def pieces(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/SIG)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x)
    return (A,)+pieces(A)

def sigma8(path):
    a=np.loadtxt(path)
    k=a[:,0]; p=a[:,1]
    q=(k>1e-4)&(p>0)&np.isfinite(p)
    k=k[q]; p=p[q]
    x=8*k; w=np.ones_like(x); m=x!=0
    w[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    y=k**3*p*w*w/(2*np.pi**2)
    return float(np.sqrt(np.trapezoid(y,x=np.log(k))))

def ini(cs2,lnAs,tau,ns,root):
    As=math.exp(lnAs)/1e10
    return textwrap.dedent(f"""\
    H0 = 70.8514
    omega_b = 0.02239952
    omega_cdm = 0.12444227328918850
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {As:.17e}
    n_s = {ns:.15g}
    tau_reio = {tau:.15g}

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = {cs2:.17g}
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = 0.0715,5.00,1.426,0.34231919445927034,1.0,0.045
    expansion_model = sdmc_full
    expansion_smg = 0.7073985893,17.925,17.775,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic

    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
    write background = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    lensing_verbose = 0
    output_verbose = 0
    """)

# Anchors make the screen interpretable and guarantee exact inclusion of the
# published current closure and the best point from the first cs2 response.
anchors=[
    ("current",math.log10(0.003),3.076,0.0544,0.964),
    ("cold1e4",math.log10(1e-4),3.076,0.0544,0.964),
    ("cold_planckAs",math.log10(1e-4),math.log(1e10*2.10e-9),0.0544,0.9649),
    ("cold_lowAs",math.log10(1e-4),3.020,0.0544,0.964),
    ("cold_lowTau",math.log10(1e-4),3.076,0.045,0.964),
    ("cold_lowAsTau",math.log10(1e-4),3.035,0.045,0.964),
]

# 64 deterministic Sobol points in the potentially important residual sector.
# c_phi^2 is sampled logarithmically toward the cold-clustering limit.
lo=np.array([-8.0,3.000,0.035,0.950],float)
hi=np.array([-3.0,3.085,0.070,0.978],float)
sob=qmc.Sobol(d=4,scramble=True,seed=20260920)
x=qmc.scale(sob.random_base2(m=6),lo,hi)
design=anchors+[(f"sobol{i+1:03d}",*row) for i,row in enumerate(x)]
(OUT/"design.json").write_text(json.dumps({
  "variables":["log10_cs2_fld","ln10As","tau_reio","n_s"],
  "low":lo.tolist(),"high":hi.tolist(),"n_design":len(design),
  "geometry":{"H0":70.8514,"omega_b":0.02239952,"omega_cdm":0.12444227328918850},
},indent=2))

rows=[]
for i,item in enumerate(design,1):
    tag,lgc,lnAs,tau,ns=item
    cs2=10.**float(lgc)
    root=str(OUT/(f"{i:03d}_{tag}_"))
    ip=OUT/(f"{i:03d}_{tag}.ini")
    ip.write_text(ini(cs2,float(lnAs),float(tau),float(ns),root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec={
      "id":tag,"index":i,"log10_cs2_fld":float(lgc),"cs2_fld":cs2,
      "ln10As":float(lnAs),"A_s":math.exp(float(lnAs))/1e10,
      "tau_reio":float(tau),"n_s":float(ns),
      "primary_log_amp":float(lnAs)-2.*float(tau),
      "returncode":cp.returncode,"status":"FAIL",
    }
    bgp=Path(root+"00_background.dat")
    clp=Path(root+"00_cl_lensed.dat")
    pkp=Path(root+"00_pk.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update({
          "min_D_all":float(bg["kin (D)"].min()),
          "min_cs2_smg_all":float(bg["c_s^2"].min()),
          "max_cs2_smg_all":float(bg["c_s^2"].max()),
          "max_abs_noslip_all":float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
        })
    if cp.returncode==0 and clp.exists() and pkp.exists():
        s=score(clp)
        rec.update({
          "status":"OK","A_planck":s[0],
          "chi2_high":s[1],"chi2_lowT":s[2],"chi2_lowE":s[3],
          "chi2_lensing":s[4],"chi2_cal":s[5],"chi2_planck":s[6],
          "sigma8":sigma8(pkp),
        })
    else:
        rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    print("COLDPRIM_POINT",json.dumps(rec,sort_keys=True),flush=True)
    rows.append(rec)

df=pd.DataFrame(rows)
df.to_csv(OUT/"coldtracker_primordial_screen.csv",index=False)
ok=df[df.status.eq("OK")].copy()
if ok.empty:
    raise SystemExit("No cold-tracker primordial points completed")

base=ok[ok.id.eq("current")].iloc[0]
for c in ["chi2_planck","chi2_high","chi2_lensing","chi2_lowT","chi2_lowE","chi2_cal"]:
    ok["delta_"+c+"_from_current"]=ok[c]-float(base[c])
ok["Aplanck_pull_sigma"]=(ok["A_planck"]-1.)/SIG
ok=ok.sort_values("chi2_planck")
ok.to_csv(OUT/"coldtracker_primordial_ranked.csv",index=False)

cols=["id","cs2_fld","ln10As","tau_reio","n_s","primary_log_amp",
      "A_planck","Aplanck_pull_sigma","chi2_planck","delta_chi2_planck_from_current",
      "chi2_high","chi2_lensing","chi2_cal","sigma8",
      "min_cs2_smg_all","max_cs2_smg_all"]
print("COLDPRIM_CURRENT",base.to_dict(),flush=True)
print("COLDPRIM_TOP20",ok.nsmallest(20,"chi2_planck")[cols].to_dict("records"),flush=True)

# Preserve a compact promotion set: top 12 + exact anchors.
keep=set(ok.nsmallest(12,"chi2_planck").id.astype(str))
keep.update(a[0] for a in anchors)
promo=ok[ok.id.astype(str).isin(keep)].sort_values("chi2_planck")
promo.to_csv(OUT/"coldtracker_primordial_shortlist.csv",index=False)
