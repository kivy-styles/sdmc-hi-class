#!/usr/bin/env python3
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/s49_structural_local"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025
H0=70.6661848725751
OB=0.02213832986768335
OC=0.12419066331610083
NS=0.9614044621046632
TAU=0.05464880801877007
LNAS=3.051477308589965
AS=math.exp(LNAS)/1e10
OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2
D0=0.34231919445927034; DFLOOR=.045
LAM=17.925; ZT=17.775

AF0=0.05709640212077647
ZC0=2.8532181778922676
W0=0.575869657304138

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
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
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
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    o=minimize_scalar(lambda A:pieces(A)[-1],bounds=(.97,1.03),method="bounded",
                      options={"xatol":1e-10})
    A=float(o.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def ini(root,AF,zc,width):
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
    parameters_smg = {AF:.17g},{zc:.17g},{width:.17g},{D0:.17g},1.0,{DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    a_ini_over_a_today_default = 1.e-8
    a_ini_test_qs_smg = 1.e-8
    pert_ic_ini_z_ref_smg = 1.e7
    a_min_stability_test_smg = 1.e-8
    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
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

AFS=[AF0-0.002,AF0,AF0+0.002]
ZCS=[ZC0-0.15,ZC0,ZC0+0.15]
WS=[W0-0.05,W0,W0+0.05]
rows=[]
for AF in AFS:
  for zc in ZCS:
    for width in WS:
      tag=f"af{AF:.6f}_z{zc:.4f}_w{width:.4f}".replace(".","p")
      root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,AF,zc,width))
      cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
      rec=dict(id=tag,A_F=AF,z_c=zc,width=width,returncode=cp.returncode,status="FAIL")
      bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
      if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
      if cp.returncode==0 and bgp.exists() and clp.exists():
        rec["stable_subluminal"]=bool(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.)
        if rec["stable_subluminal"]:
          rec.update(score(clp)); rec["status"]="OK"
      if rec["status"]!="OK": rec["error"]=cp.stdout[-900:].replace("\n"," | ")
      rows.append(rec); print("S49_STRUCT_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
base=df[(np.isclose(df.A_F,AF0))&(np.isclose(df.z_c,ZC0))&(np.isclose(df.width,W0))].iloc[0]
for c in ["chi2_total","chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal"]:
    df["delta_"+c+"_vs_base"]=df[c]-base[c]
df.to_csv(OUT/"sobol049_structural_local.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
best=ok.nsmallest(12,"chi2_total")
best.to_csv(OUT/"sobol049_structural_shortlist.csv",index=False)
print("S49_STRUCT_BASE",json.dumps(base.to_dict(),default=float,sort_keys=True),flush=True)
print("S49_STRUCT_BEST",json.dumps(best.to_dict("records"),default=float,sort_keys=True),flush=True)
