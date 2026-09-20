#!/usr/bin/env python3
"""
Combined acoustic-degeneracy + localized-braiding refinement.

Ordinary cosmology is fixed to the best acoustic-screen point (sobol002):
H0=69.7095092787552, omega_b=0.02211778594478965,
omega_cdm=0.1258596659367904, n_s=0.9637720508351921,
tau=0.06456347454525531, ln(1e10 As)=3.0773358772153037.

We scan a positive localized departure from exact No-Slip because the first
sign/localization diagnostic found all tested negative amplitudes unstable,
while A_lens=+0.02, z_lens=1.5, sigma_lens=0.8 was stable and lowered the
Planck-lite score.

This remains an EFT-level diagnostic. Finalists require:
  (1) covariant G3 reconstruction,
  (2) full Plik nuisance profiling,
  (3) raw DESI DR1 full-shape,
  (4) separate SN compilation checks.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/acoustic_braid"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025

H0=69.7095092787552
OB=0.02211778594478965
OC=0.1258596659367904
NS=0.9637720508351921
TAU=0.06456347454525531
LNAS=3.0773358772153037
AS=math.exp(LNAS)/1e10
OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2

AF=0.0715; ZC=5.; WIDTH=1.426
D0=0.34231919445927034; DFLOOR=.045
LAM=17.925; ZT=17.775

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=ell.max()+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); te=np.zeros(n); ee=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp},a

def score(path):
    harr,dls,a=load_cls(path)
    def pieces(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pieces(A)[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pieces(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5]),a

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]; q=(k>1e-4)&(p>0)&np.isfinite(p)
    k=k[q]; p=p[q]; x=8*k
    W=np.ones_like(x); m=x!=0
    W[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    return float(np.sqrt(np.trapezoid(k**3*p*W**2/(2*np.pi**2),x=np.log(k))))

def acoustic(bg,th):
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy()
    cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i]); cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross: return {}
    zd=float(cross[0])
    rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    rd=float(np.interp(zd,bg.z,bg["comov.snd.hrz."]))
    dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return dict(z_star=zstar,z_drag=zd,rs_star=rs,rd=rd,ell_A=math.pi*dm/rs)

def ini(root,A,zl,sl):
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {OB:.15g}
    omega_cdm = {OC:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {AS:.17e}
    n_s = {NS:.15g}
    tau_reio = {TAU:.15g}

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_lens_braiding
    parameters_smg = {AF},{ZC},{WIDTH},{D0:.17g},1.0,{DFLOOR},{A:.17g},{zl:.17g},{sl:.17g}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
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

amps=[0.0,0.004,0.008,0.012,0.016,0.020,0.024]
centers=[1.0,1.25,1.5,1.75,2.0]
widths=[0.6,0.8,1.0,1.2]
design=[(0.0,1.5,0.8)]
for A in amps[1:]:
    for z in centers:
        for w in widths:
            design.append((A,z,w))

rows=[]; control=None; control_a=None
for i,(A,zl,sl) in enumerate(design,1):
    tag=f"a{A:.3f}_z{zl:.2f}_s{sl:.2f}".replace(".","p")
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,A,zl,sl))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,index=i,A_lens=A,z_lens=zl,sigma_lens=sl,returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat")
    clp=Path(root+"00_cl_lensed.dat"); pkp=Path(root+"00_pk.dat")
    if bgp.exists():
        bg=tab(bgp); zz=bg.z.to_numpy()
        rec.update(
            min_D=float(bg["kin (D)"].min()),
            min_cs2=float(bg["c_s^2"].min()),
            max_cs2=float(bg["c_s^2"].max()),
            max_abs_slip_driver=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
            alphaM_zlens=float(np.interp(zl,zz,bg["M2_running_smg"])),
            alphaB_zlens=float(np.interp(zl,zz,bg["braiding_smg"])),
        )
    if cp.returncode==0 and bgp.exists() and thp.exists() and clp.exists() and pkp.exists():
        stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.0)
        rec["stable_subluminal"]=bool(stable)
        if stable:
            bg=tab(bgp); th=tab(thp); sc,a=score(clp)
            rec.update(sc); rec.update(acoustic(bg,th)); rec["sigma8"]=sigma8(pkp); rec["status"]="OK"
            if A==0:
                control=sc.copy(); control_a=a.copy()
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    rows.append(rec)
    print("ACOUSTIC_BRAID_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
zero=df[df.A_lens==0].iloc[0]
for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
    df["delta_"+k+"_vs_zero"]=df[k]-zero[k]
df.to_csv(OUT/"acoustic_braiding_refine.csv",index=False)

ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
print("ACOUSTIC_BRAID_ZERO",zero.to_dict(),flush=True)
print("ACOUSTIC_BRAID_BEST_TOTAL",ok.nsmallest(20,"chi2_total").to_dict("records"),flush=True)
print("ACOUSTIC_BRAID_BEST_LENS",ok.nsmallest(20,"chi2_lensing").to_dict("records"),flush=True)
print("ACOUSTIC_BRAID_STABLE_COUNT",len(ok),"OF",len(df),flush=True)
