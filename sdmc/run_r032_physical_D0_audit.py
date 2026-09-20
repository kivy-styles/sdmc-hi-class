#!/usr/bin/env python3
"""
Audit the hidden D0 coupling at the best final-r032 physical-lensing transition.

The transition point is held fixed at the derived Sobol optimum:
  A_F=0.06616764205507934
  z_c=3.0664408076554537
  width=0.5638187950477004

We scan D0 at fixed Dfloor=0.045, including:
  * the inherited old-r032 D0=0.34231919445927034,
  * the D(z=0)-matched value D0=0.28860640116889824,
  * nearby values.

This cleanly separates the effect of F(a) timing from the effect of changing the
kinetic normalization. Planck likelihood is the same Plik-lite screening stack
used in the other diagnostics.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/physical_D0_audit"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.8514; OB=0.02239952; OC=0.12444227328918850
AS=math.exp(3.076)/1e10; NS=0.964; TAU=0.0544
OR=4.17998772e-5; OX=1.-(OB+OC+OR)/(H0/100.)**2
LAM=17.925; ZT=17.775; DFLOOR=0.045
AF=0.06616764205507934
ZC=3.0664408076554537
WIDTH=0.5638187950477004
D_TODAY_BASE=0.31146874444858513
Nc=-math.log1p(ZC)
S0=0.5*(1+math.tanh((0-Nc)/(2*WIDTH)))
D0_MATCH=(D_TODAY_BASE-DFLOOR)/S0
D0_OLD=0.34231919445927034
D0_VALUES=[0.24,0.265,D0_MATCH,0.305,0.325,D0_OLD,0.36]

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
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def pcs(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda x:pcs(x)[-1],bounds=(0.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]
    q=(k>1e-4)&(p>0)&np.isfinite(p); k=k[q]; p=p[q]
    x=8*k; W=np.ones_like(x); m=x!=0
    W[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    return float(np.sqrt(np.trapezoid(k**3*p*W**2/(2*np.pi**2),x=np.log(k))))

def lcdm_ini(root):
    return textwrap.dedent(f"""\
    H0 = 67.36
    omega_b = 0.02237
    omega_cdm = 0.1200
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = 2.10e-9
    n_s = 0.9649
    tau_reio = 0.0544
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
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

def ini(root,D0):
    return textwrap.dedent(f"""\
    H0 = {H0}
    omega_b = {OB}
    omega_cdm = {OC}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {AS:.17e}
    n_s = {NS}
    tau_reio = {TAU}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF:.17g}, {ZC:.17g}, {WIDTH:.17g}, {D0:.17g}, 1.0, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g}, {LAM}, {ZT}, 0.5, 0.01105624999, 0.25, 0.01951933685, 1.5
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

lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,
                  stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
base=score(Path(lr+"00_cl_lensed.dat"))
print("D0_AUDIT_LCDM",base,flush=True)
print("D0_AUDIT_INVARIANT",
      json.dumps(dict(S0=S0,D_today_target=D_TODAY_BASE,D0_matched=D0_MATCH,
                      D0_inherited=D0_OLD,
                      D_today_inherited=DFLOOR+D0_OLD*S0),sort_keys=True),flush=True)

rows=[]
for i,D0 in enumerate(D0_VALUES,1):
    tag=("matched" if abs(D0-D0_MATCH)<1e-12 else
         "inherited" if abs(D0-D0_OLD)<1e-12 else f"D0_{D0:.6f}")
    root=str(OUT/(f"{i:02d}_{tag}_")); ip=OUT/(f"{i:02d}_{tag}.ini")
    ip.write_text(ini(root,D0))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec=dict(id=tag,D0=float(D0),S0=S0,D_today=float(DFLOOR+D0*S0),
             D_today_target=D_TODAY_BASE,returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat"); pkp=Path(root+"00_pk.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()))
    if cp.returncode==0 and clp.exists() and pkp.exists():
        stable=rec.get("min_D",0)>0 and rec.get("min_cs2",0)>0 and rec.get("max_cs2",2)<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            sc=score(clp); rec.update(sc); rec["status"]="OK"; rec["sigma8"]=sigma8(pkp)
            for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
                rec["delta_"+k]=rec[k]-base[k]
    if rec["status"]!="OK": rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    rows.append(rec); print("D0_AUDIT_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"physical_D0_audit.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
if not ok.empty:
    ok.sort_values("chi2_total").to_csv(OUT/"physical_D0_audit_shortlist.csv",index=False)
    print("D0_AUDIT_BEST",ok.sort_values("chi2_total").to_dict("records"),flush=True)
