#!/usr/bin/env python3
from pathlib import Path
import csv, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/tracker_cs2")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255
CAL_SIGMA=0.0025
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

CS2_VALUES=[1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1.0]

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
    A=float(opt.x)
    return (A,)+pieces(A)

def sigma8(path):
    a=np.loadtxt(path)
    k=a[:,0]; p=a[:,1]
    q=(k>1e-4)&(p>0)&np.isfinite(p)
    k=k[q]; p=p[q]
    x=8*k
    w=np.ones_like(x)
    m=x!=0
    w[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    y=k**3*p*w*w/(2*np.pi**2)
    return float(np.sqrt(np.trapezoid(y,x=np.log(k))))

def ini(cs2,root):
    return textwrap.dedent(f"""\
    H0 = 70.8514
    omega_b = 0.02239952
    omega_cdm = 0.12444227328918850
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = 2.16715426219737e-9
    n_s = 0.964
    tau_reio = 0.0544

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = {cs2:.17g}
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = 0.0715, 5.00, 1.426, 0.34231919445927034, 1.0, 0.045
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

lroot=str(OUT/"lcdm_")
lip=OUT/"lcdm.ini"; lip.write_text(lcdm_ini(lroot))
cp=subprocess.run(["./class",str(lip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lcdm=score(Path(lroot+"00_cl_lensed.dat"))
print("TRACKER_CS2_LCDM",lcdm,flush=True)

rows=[]
for i,cs2 in enumerate(CS2_VALUES,1):
    root=str(OUT/(f"c{i:02d}_"))
    ip=OUT/(f"c{i:02d}.ini"); ip.write_text(ini(cs2,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"cs2_fld":cs2,"returncode":cp.returncode,"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    pkp=Path(root+"00_pk.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update({
          "min_D_all":float(bg["kin (D)"].min()),
          "min_cs2_smg_all":float(bg["c_s^2"].min()),
          "max_cs2_smg_all":float(bg["c_s^2"].max()),
          "max_abs_noslip_all":float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
        })
    if cp.returncode==0 and clp.exists() and bgp.exists() and pkp.exists():
        s=score(clp)
        rec.update({
          "status":"OK","A_planck":s[0],"chi2_high":s[1],"chi2_lowT":s[2],
          "chi2_lowE":s[3],"chi2_lensing":s[4],"chi2_cal":s[5],
          "chi2_planck":s[6],"delta_planck":s[6]-lcdm[-1],
          "sigma8":sigma8(pkp),
        })
    else:
        rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    print("TRACKER_CS2_POINT",rec,flush=True)
    rows.append(rec)

df=pd.DataFrame(rows)
df.to_csv(OUT/"tracker_cs2_response.csv",index=False)
ok=df[df.status.eq("OK")].copy()
if ok.empty: raise SystemExit("No tracker sound-speed point completed")
base=ok.iloc[(ok.cs2_fld-0.003).abs().argsort()[:1]].iloc[0]
for c in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_planck","sigma8","A_planck"]:
    ok["d_"+c]=ok[c]-float(base[c])
ok=ok.sort_values("chi2_planck")
ok.to_csv(OUT/"tracker_cs2_ranked.csv",index=False)
print("TRACKER_CS2_BASE",base.to_dict(),flush=True)
print("TRACKER_CS2_RANKED",ok.to_dict("records"),flush=True)
