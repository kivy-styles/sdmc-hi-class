#!/usr/bin/env python3
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/noslip_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.66173072; OB=0.02224; OC=0.12444227328918850
AS=math.exp(3.074)/1e10; NS=0.964; TAU=0.061
OX=0.7061451696512845
D0=0.34231919445927034; DFLOOR=0.045
LAMBDA=17.925; ZT=17.775

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
    def pieces(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def ini(AF,zc,width,root):
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
    parameters_smg = {AF:.12g},{zc:.12g},{width:.12g},{D0:.17g},1.0,{DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX},{LAMBDA},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic

    modes = s
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
    output_background_smg = 3
    write background = yes
    write thermodynamics = yes
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

# LCDM reference
lr=str(OUT/"lcdm_"); ip=OUT/"lcdm.ini"; ip.write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2500:])
lcdm=score(Path(lr+"00_cl_lensed.dat"))
print("NOSLIP_REFINE_LCDM",lcdm,flush=True)

# Joint physical refinement around the stable low-AF solution, now using
# the acoustic/amplitude ordinary-cosmology candidate.
AFS=[0.0400,0.0450,0.0475,0.0500,0.0525,0.0550]
ZCS=[3.0,3.4,3.8,4.2,4.6]
design=[("ordinary_only",0.0715,5.0,1.426),
        ("structural_seed",0.0600,5.0,0.9573146853146852)]
for af in AFS:
    wc=af/(0.0715/1.426)
    for fac in [0.50,0.60,0.65,0.70]:
        w=wc*fac
        for zc in ZCS:
            design.append((f"joint_af{af:.4f}_zc{zc:.2f}_w{w:.3f}",af,zc,w))

rows=[]
for i,(tag,af,zc,w) in enumerate(design,1):
    root=str(OUT/(f"{i:03d}_{tag}_"))
    ip=OUT/(f"{i:03d}_{tag}.ini"); ip.write_text(ini(af,zc,w,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec={"id":tag,"A_F":af,"z_c":zc,"width":w,"returncode":cp.returncode}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   F0=float(bg.iloc[np.argmin(np.abs(bg.z.to_numpy()))]["M*^2_smg"]),
                   alphaM0=float(bg.iloc[np.argmin(np.abs(bg.z.to_numpy()))]["M2_running_smg"]),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    if cp.returncode==0 and clp.exists():
        stable=rec.get("min_D",-1)>0 and rec.get("min_cs2",-1)>0 and rec.get("max_cs2",2)<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            sc=score(clp); rec.update(sc)
            for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
                rec["delta_"+k]=sc[k]-lcdm[k]
    else:
        rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    rows.append(rec)
    print("NOSLIP_REFINE_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"noslip_profile_refine.csv",index=False)
ok=df[(df.returncode==0)&(df.stable_subluminal==True)&np.isfinite(df.chi2_total)].copy()
print("NOSLIP_REFINE_BEST_TOTAL",ok.nsmallest(12,"chi2_total").to_dict("records"),flush=True)
print("NOSLIP_REFINE_BEST_HIGH",ok.nsmallest(12,"chi2_high").to_dict("records"),flush=True)
print("NOSLIP_REFINE_BEST_LENSING",ok.nsmallest(12,"chi2_lensing").to_dict("records"),flush=True)
