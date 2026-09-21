#!/usr/bin/env python3
"""
Axis-split and acoustic-scale diagnostic for the equality hypothesis.

Tests at fixed geometry-refined SDMC structure:
  1) current geometry-refined densities;
  2) omega_m fixed to LCDM local021 while omega_b remains SDMC;
  3) omega_b fixed to LCDM local021 while total omega_m remains SDMC;
  4) both omega_b and omega_m fixed to LCDM local021;
  5) H0 scan at the combined equality-fixed densities.

For every successful model we compute:
  * Planck TTTEEE-lite + lowT + lowE + lensing + calibration score;
  * exact background matter-radiation equality crossing;
  * recombination visibility peak z_*;
  * acoustic scale ell_A = pi D_M(z_*) / r_s(z_*).
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/zeq_axis_acoustic"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0G=70.5653567390982
NS=0.9625227132590487
TAU=0.055202901571989066
LNAS=3.0598696043919773
AS=math.exp(LNAS)/1e10

OBG=0.022011983189284802
OCG=0.12404331885203719
WMG=OBG+OCG

OBL=0.0224063693780079
OCL=0.1182415105248335
WML=OBL+OCL

AF=0.06017362505197525
ZC=2.827464461401105
WIDTH=0.5514493708219379
D0=0.34231919445927034
DF=0.045
LAM=17.7
ZT=17.1

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

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    harr,dls=load_cls(path)
    def pc(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def crossing(z,y):
    order=np.argsort(z); z=np.asarray(z)[order]; y=np.asarray(y)[order]
    q=y-1.
    ii=np.where(q[:-1]*q[1:]<=0)[0]
    if not len(ii): return float("nan")
    i=ii[-1]
    if y[i+1]==y[i]: return float(z[i])
    return float(z[i]+(1-y[i])*(z[i+1]-z[i])/(y[i+1]-y[i]))

def derived(bg,th):
    b=bg.sort_values("z")
    matter=b["(.)rho_b"].to_numpy()+b["(.)rho_cdm"].to_numpy()
    rad=b["(.)rho_g"].to_numpy()+b["(.)rho_ur"].to_numpy()
    zeq=crossing(b["z"].to_numpy(),matter/rad)
    gcol="g [Mpc^-1]"
    imax=int(np.argmax(th[gcol].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,b["z"],b["comov. dist."]))
    rs=float(np.interp(zstar,b["z"],b["comov.snd.hrz."]))
    ellA=float(np.pi*DM/rs)
    return dict(z_eq_exact=zeq,z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=ellA)

def ini(root,H0,ob,oc):
    h=H0/100.
    ox=1.-(ob+oc+OR)/(h*h)
    return textwrap.dedent(f"""\
    H0 = {H0:.17g}
    omega_b = {ob:.17g}
    omega_cdm = {oc:.17g}
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
    parameters_smg = {AF},{ZC},{WIDTH},{D0},1.0,{DF}
    expansion_model = sdmc_full
    expansion_smg = {ox:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
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

def run(tag,H0,ob,oc):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,H0,ob,oc))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=tag,H0=H0,omega_b=ob,omega_cdm=oc,omega_m=ob+oc,status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and thp.exists() and clp.exists():
        bg=table(bgp); th=table(thp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()))
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(pscore(clp)); rec.update(derived(bg,th)); rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    print("ZEQAXIS_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[]
# 1. Current geometry-refined milestone.
rows.append(run("geom_current",H0G,OBG,OCG))
# 2. Equality-only: match LCDM omega_m, retain SDMC omega_b.
rows.append(run("wm_only",H0G,OBG,WML-OBG))
# 3. Baryon-only: match LCDM omega_b, retain SDMC total omega_m.
rows.append(run("ob_only",H0G,OBL,WMG-OBL))
# 4. Both equality-defining densities matched.
rows.append(run("wm_ob_both",H0G,OBL,OCL))

# 5. H0/acoustic-geometry compensation scan at matched omega_b, omega_m.
for H0 in [64.,65.,66.,67.,68.,68.5,69.,69.5,70.,70.5653567390982,71.,71.5,72.,73.,74.]:
    rows.append(run(("eq_h%06.3f"%H0).replace(".","p"),H0,OBL,OCL))

df=pd.DataFrame(rows)
base=df[df.id.eq("geom_current")].iloc[0]
for col in ["chi2_planck","chi2_high","chi2_lowT","chi2_lensing","ell_A"]:
    if col in df:
        df["delta_"+col+"_vs_geom"]=df[col]-float(base[col])
df.to_csv(OUT/"zeq_axis_acoustic.csv",index=False)
ok=df[df.status.eq("OK")].copy()
print("ZEQAXIS_BEST_PLANCK",ok.nsmallest(12,"chi2_planck").to_dict("records"),flush=True)
print("ZEQAXIS_CLOSEST_ELLA",ok.assign(della=np.abs(ok.ell_A-float(base.ell_A))).nsmallest(12,"della").to_dict("records"),flush=True)
summary={
 "geom":df[df.id.eq("geom_current")].iloc[0].to_dict(),
 "wm_only":df[df.id.eq("wm_only")].iloc[0].to_dict(),
 "ob_only":df[df.id.eq("ob_only")].iloc[0].to_dict(),
 "both":df[df.id.eq("wm_ob_both")].iloc[0].to_dict(),
 "best_eq_h0":ok[ok.id.str.startswith("eq_h")].nsmallest(1,"chi2_planck").iloc[0].to_dict()
}
(OUT/"zeq_axis_acoustic_summary.json").write_text(json.dumps(summary,indent=2))
print("ZEQAXIS_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
