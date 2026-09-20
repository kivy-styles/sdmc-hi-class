#!/usr/bin/env python3
"""
Profile primordial amplitude/tilt on the best stable retimed No-Slip gravity
trajectory found by the physical Planck-gap scan:
  A_F=0.045, z_c=3.5, width=0.5.

The background/ruler sector remains the exact r032 structural point.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/latesharp_prim"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.8514; OB=0.02239952; OC=0.12444227328918850
AF=0.045; ZC=3.5; WIDTH=0.5
D0=0.34231919445927034; DFLOOR=0.045
LAM=17.925; ZT=17.775; OR=4.17998772e-5
OX=1.-(OB+OC+OR)/(H0/100.)**2
LNAS=[3.050,3.054,3.058,3.062,3.066,3.070,3.074,3.078]
NSS=[0.960,0.962,0.964,0.966]\nTAUS=[0.058,0.060,0.062,0.064,0.066]

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def pcs(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp)); cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda a:pcs(float(a))[-1],bounds=(0.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(op.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def sigma8(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]; q=(k>1e-4)&(p>0)
    k=k[q]; p=p[q]; x=8*k
    W=np.ones_like(x); m=x!=0
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

def ini(root,lnAs,ns,tau):
    return textwrap.dedent(f"""\
    H0 = {H0}
    omega_b = {OB}
    omega_cdm = {OC}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {math.exp(lnAs)/1e10:.17e}
    n_s = {ns:.16g}
    tau_reio = {tau:.16g}
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
base=score(Path(lr+"00_cl_lensed.dat"))
print("LATESHARP_PRIM_LCDM",base,flush=True)

rows=[]
for tau in TAUS:
  for ns in NSS:
    for lnAs in LNAS:
      tag=f"t{tau:.3f}_ns{ns:.3f}_la{lnAs:.3f}".replace(".","p")
      root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,lnAs,ns,tau))
      cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
      rec=dict(id=tag,tau_reio=tau,n_s=ns,ln10As=lnAs,Q=lnAs-2*tau,
               A_s=math.exp(lnAs)/1e10,returncode=cp.returncode,status="FAIL")
      bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat"); pkp=Path(root+"00_pk.dat")
      if bgp.exists():
          bg=table(bgp)
          rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                     max_cs2=float(bg["c_s^2"].max()),F0=float(bg.iloc[np.argmin(np.abs(bg.z.to_numpy()))]["M*^2_smg"]))
      if cp.returncode==0 and clp.exists() and pkp.exists():
          stable=rec.get("min_D",0)>0 and rec.get("min_cs2",0)>0 and rec.get("max_cs2",2)<=1
          rec["stable_subluminal"]=bool(stable)
          if stable:
              sc=score(clp); rec.update(sc); rec["sigma8"]=sigma8(pkp); rec["status"]="OK"
              for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
                  rec["delta_"+k]=rec[k]-base[k]
      if rec["status"]!="OK": rec["error"]=cp.stdout[-900:].replace("\n"," | ")
      rows.append(rec); print("LATESHARP_PRIMTAU_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"latesharp_gravity_primordial_tau_profile.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
print("LATESHARP_PRIMTAU_BEST_TOTAL",ok.nsmallest(20,"chi2_total").to_dict("records"),flush=True)
print("LATESHARP_PRIMTAU_BEST_HIGH",ok.nsmallest(20,"chi2_high").to_dict("records"),flush=True)
