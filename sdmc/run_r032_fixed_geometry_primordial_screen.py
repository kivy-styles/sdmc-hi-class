#!/usr/bin/env python3
"""
Fixed-background primordial/reionization compensation screen for r032.

H0, omega_b, omega_cdm and all SDMC background/modified-gravity parameters are
held at the exact leader. We scan only (tau, n_s, Q), where
Q = ln(1e10 A_s) - 2 tau, to separate the primary-CMB amplitude direction from
the physical A_s/lensing amplitude. This directly tests whether Planck can be
improved without sacrificing the DESI-favored background geometry.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/fixed_primordial"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025
H0=70.8514; OB=.02239952; OC=.12444227328918850
OR=4.17998772e-5; OX=1.-(OB+OC+OR)/(H0/100.)**2
AF=.0715; ZC=5.; WIDTH=1.426; D0=.34231919445927034; DFLOOR=.045
LAM=17.925; ZT=17.775

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
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
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp},a

def score(path):
    harr,dls,a=load_cls(path)
    def pcs(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda x:pcs(float(x))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5]),a

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]; q=(k>1e-4)&(p>0)
    k=k[q]; p=p[q]; x=8*k; W=np.ones_like(x); m=x!=0
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

def ini(root,ns,tau,Q):
    lnAs=Q+2*tau; As=math.exp(lnAs)/1e10
    return textwrap.dedent(f"""\
    H0 = {H0}
    omega_b = {OB}
    omega_cdm = {OC}
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
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF},{ZC},{WIDTH},{D0:.17g},1.0,{DFLOOR}
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
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
base,lcl=score(Path(lr+"00_cl_lensed.dat"))
print("FIXPRIM_LCDM",base,flush=True)

# Include exact current point plus a dense grid in physically meaningful
# amplitude coordinates. Q controls A_s e^{-2 tau}; tau can then lift A_s and
# lensing without restoring the same primary-CMB amplitude.
design=[("current",.964,.0544,3.076-2*.0544)]
for tau in [.056,.060,.064,.068,.072]:
  for ns in [.960,.964,.968,.972]:
    for Q in [2.948,2.953,2.958,2.963,2.968]:
      design.append((f"g_t{tau:.3f}_n{ns:.3f}_q{Q:.3f}",ns,tau,Q))

rows=[]
for i,(tag,ns,tau,Q) in enumerate(design,1):
    lnAs=Q+2*tau; As=math.exp(lnAs)/1e10
    root=str(OUT/(f"{i:03d}_{tag}_")); ip=OUT/(f"{i:03d}_{tag}.ini")
    ip.write_text(ini(root,ns,tau,Q))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,index=i,n_s=ns,tau_reio=tau,Q=Q,ln10As=lnAs,A_s=As,returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat"); pkp=Path(root+"00_pk.dat")
    if bgp.exists():
      bg=table(bgp); rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
    if cp.returncode==0 and clp.exists() and pkp.exists():
      stable=rec.get("min_D",0)>0 and rec.get("min_cs2",0)>0 and rec.get("max_cs2",2)<=1
      rec["stable_subluminal"]=bool(stable)
      if stable:
        sc,a=score(clp); rec.update(sc); rec["sigma8"]=sigma8(pkp); rec["status"]="OK"
        for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
          rec["delta_"+k]=rec[k]-base[k]
        for L in [40,100,400,1000]:
          e=a[:,0].astype(int); le=lcl[:,0].astype(int); ii=np.where(e==L)[0]; jj=np.where(le==L)[0]
          if len(ii) and len(jj): rec[f"pp_ratio_L{L}"]=float(a[ii[0],5]/lcl[jj[0],5])
        for L in [200,500,1000,2000]:
          e=a[:,0].astype(int); le=lcl[:,0].astype(int); ii=np.where(e==L)[0]; jj=np.where(le==L)[0]
          if len(ii) and len(jj): rec[f"tt_ratio_l{L}"]=float(a[ii[0],1]/lcl[jj[0],1])
    if rec["status"]!="OK": rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    rows.append(rec); print("FIXPRIM_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"fixed_primordial_screen.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
print("FIXPRIM_BEST_TOTAL",ok.nsmallest(20,"chi2_total").to_dict("records"),flush=True)
print("FIXPRIM_BEST_LENS_CAL",ok.assign(lc=ok.chi2_lensing+ok.chi2_cal).nsmallest(20,"lc").to_dict("records"),flush=True)
