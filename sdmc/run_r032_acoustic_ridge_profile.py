#!/usr/bin/env python3
"""
Profile the narrow SDMC acoustic ridge instead of sampling H0 independently.

For each (omega_b, omega_cdm) and each target acoustic scale ell_A, solve H0
with a background+thermodynamics-only CLASS call, then evaluate the full
Planck-lite + low-l + lensing score at that on-ridge point.

This is a screening profile at frozen r032 structural/gravity parameters.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import brentq, minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/acoustic_ridge"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
NS=0.964; LNAS=3.076; TAU=0.0544
AF=0.0715; ZC=5.; WIDTH=1.426; D0=0.34231919445927034; DFLOOR=0.045
LAM=17.925; ZT=17.775; OR=4.17998772e-5
OBS=[0.0220,0.0222,0.0224,0.0226,0.0228]
OCS=[0.1180,0.1200,0.1220,0.1244422732891885,0.1260]
TARGETS=[301.45,301.585,301.72]

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

def acoustic(root):
    bg=table(Path(root+"00_background.dat"))
    th=table(Path(root+"00_thermodynamics.dat"))
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return math.pi*dm/rs,zstar,rs

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
    op=minimize_scalar(lambda a:pcs(float(a))[-1],bounds=(0.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def sdmc_ini(root,H0,ob,oc,full):
    h=H0/100.; Ox=1.-(ob+oc+OR)/(h*h)
    common=f"""
    H0 = {H0:.16g}
    omega_b = {ob:.16g}
    omega_cdm = {oc:.16g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {math.exp(LNAS)/1e10:.17e}
    n_s = {NS}
    tau_reio = {TAU}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF},{ZC},{WIDTH},{D0:.17g},1.0,{DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {Ox:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
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
    """
    if full:
        common += """
        modes = s
        output = tCl,pCl,lCl
        lensing = yes
        l_max_scalars = 3000
        """
    return textwrap.dedent(common)

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

lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2500:])
base=score(Path(lr+"00_cl_lensed.dat"))
print("RIDGE_LCDM",base,flush=True)

counter=0
def eval_ell(H0,ob,oc,key):
    global counter
    counter+=1
    root=str(OUT/(f"probe_{key}_{counter:04d}_"))
    ip=OUT/(f"probe_{key}_{counter:04d}.ini")
    ip.write_text(sdmc_ini(root,H0,ob,oc,False))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
    if cp.returncode:
        raise RuntimeError(cp.stdout[-1600:])
    ea,_,_=acoustic(root)
    return ea

rows=[]
for ob in OBS:
  for oc in OCS:
    key=f"ob{ob:.5f}_oc{oc:.5f}".replace(".","p")
    samples=[]
    for H in np.linspace(66.,74.,17):
        try:
            samples.append((float(H),eval_ell(float(H),ob,oc,key)))
        except Exception:
            pass
    if len(samples)<2:
        print("RIDGE_BRACKET_FAIL",key,"fewer than two stable H0 samples",flush=True); continue
    for target in TARGETS:
      rec=dict(omega_b=ob,omega_cdm=oc,target_ell_A=target,status="NO_ROOT",
               stable_H0_min=min(h for h,e in samples),stable_H0_max=max(h for h,e in samples))
      brackets=[]
      for (h1,e1),(h2,e2) in zip(samples[:-1],samples[1:]):
          if (e1-target)*(e2-target)<=0 and e1!=e2:
              brackets.append((h1,h2))
      if not brackets:
        rows.append(rec); print("RIDGE_POINT",json.dumps(rec,sort_keys=True),flush=True); continue
      try:
        h1,h2=min(brackets,key=lambda q:abs(0.5*(q[0]+q[1])-70.8514))
        f=lambda H: eval_ell(H,ob,oc,key)-target
        Hroot=float(brentq(f,h1,h2,xtol=2e-5,rtol=1e-9,maxiter=32))
        root=str(OUT/(f"final_{key}_e{target:.3f}_".replace(".","p")))
        ip=OUT/(f"final_{key}_e{target:.3f}.ini".replace(".","p"))
        ip.write_text(sdmc_ini(root,Hroot,ob,oc,True))
        cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
        if cp.returncode: raise RuntimeError(cp.stdout[-1800:])
        bgfiles=sorted(Path(root).parent.glob(Path(root).name+"*background.dat"))
        thfiles=sorted(Path(root).parent.glob(Path(root).name+"*thermodynamics.dat"))
        clfiles=sorted(Path(root).parent.glob(Path(root).name+"*cl_lensed.dat"))
        if len(bgfiles)!=1 or len(thfiles)!=1 or len(clfiles)!=1:
            raise RuntimeError(f"final output discovery failed bg={bgfiles} th={thfiles} cl={clfiles}")
        bg=table(bgfiles[0]); th=table(thfiles[0])
        zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
        rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
        dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
        ea=math.pi*dm/rs; zs=zstar
        stable=(float(bg["kin (D)"].min())>0 and float(bg["c_s^2"].min())>0 and float(bg["c_s^2"].max())<=1)
        rec.update(status="OK" if stable else "UNSTABLE",H0=Hroot,ell_A=ea,z_star=zs,rs_star=rs,
                   Omega_m0=(ob+oc)/(Hroot/100.)**2,
                   min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
        if stable:
            sc=score(clfiles[0]); rec.update(sc)
            for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
                rec["delta_"+k]=rec[k]-base[k]
      except Exception as ex:
        rec["status"]="FAIL"; rec["error"]=str(ex)
      rows.append(rec); print("RIDGE_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"acoustic_ridge_profile.csv",index=False)
ok=df[df.status=="OK"].copy()
if ok.empty: raise SystemExit("no acoustic-ridge points")
print("RIDGE_BEST_TOTAL",ok.nsmallest(20,"chi2_total").to_dict("records"),flush=True)
print("RIDGE_BEST_HIGH",ok.nsmallest(20,"chi2_high").to_dict("records"),flush=True)
