#!/usr/bin/env python3
"""
Wide exact-leader structural stability ridge:
scan (AF, zc, width) with the corrected r032 background held fixed.
Stage A: background/stability only over a broad Sobol design.
Stage B: Planck-lite score the stable subluminal trajectories with the
smallest present-day Planck-mass excursion F0, plus the current point.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/sharp_refine"); OUT.mkdir(parents=True,exist_ok=True)
H0=70.8514; OB=0.02239952; OC=0.12444227328918850
AS=2.16715426219737e-9; NS=0.964; TAU=0.0544
OX=0.7073985893; D0=0.34231919445927034; DFLOOR=0.045
TCMB=2.7255; SIG=0.0025
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def ini(AF,zc,w,root,full=False):
    outputs="modes = s\noutput = tCl,pCl,lCl\nlensing = yes\nl_max_scalars = 3000\n" if full else ""
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
    parameters_smg = {AF:.17g}, {zc:.17g}, {w:.17g}, {D0:.17g}, 1.0, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX},17.925,17.775,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    {outputs}
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

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),dict(tt=tt,te=te,ee=ee,pp=pp)

def score(path):
    h,d=load_cls(path)
    def p(A):
        A=float(A)
        ch=float(high.chi_squared(h,A_planck=A))
        ct=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(d,**lp))
        cp=((A-1.)/SIG)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    o=minimize_scalar(lambda A:p(A)[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(o.x); return (A,)+p(A)

# Reference.
lr=str(OUT/"lcdm_"); li=OUT/"lcdm.ini"; li.write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(li)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lscore=score(lr+"00_cl_lensed.dat")

# Broad transformed design. width is sampled logarithmically; zc linearly.
sob=qmc.Sobol(d=3,scramble=True,seed=20260920)
u=sob.random_base2(m=7)  # 128
AF=0.025 + u[:,0]*(0.047-0.025)
ZC=3.0 + u[:,1]*(5.6-3.0)
W=np.exp(np.log(0.24)+u[:,2]*(np.log(0.60)-np.log(0.24)))
design=[
 ("current",0.0715,5.0,1.426),
 ("sharp_s195",0.036947266645729546,4.2729806587100025,0.3756495761351945),
 ("sharp_s107",0.0444651716761291,3.40291923135519,0.35471224024498715)
]
design += [(f"s{i+1:03d}",float(a),float(z),float(w)) for i,(a,z,w) in enumerate(zip(AF,ZC,W))]

rows=[]
for i,(tag,af,zc,w) in enumerate(design,1):
    root=str(OUT/f"bg_{i:03d}_"); ip=OUT/f"bg_{i:03d}.ini"; ip.write_text(ini(af,zc,w,root,False))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
    rec=dict(id=tag,AF=af,z_c=zc,width=w,status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat")
    if cp.returncode==0 and bgp.exists():
        bg=tab(bgp); z=bg["z"].to_numpy(); iz=int(np.argmin(np.abs(z)))
        rec.update(status="OK",F0=float(bg["M*^2_smg"].to_numpy()[iz]),
                   alphaM0=float(bg["M2_running_smg"].to_numpy()[iz]),
                   min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        rec["stable_subluminal"]=bool(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.0)
    else:
        m=re.search(r"minimum c_s\^2=([+\-0-9.eE]+) at a=([+\-0-9.eE]+)",cp.stdout)
        if m:
            rec["failed_min_cs2"]=float(m.group(1)); rec["failed_a"]=float(m.group(2))
        rec["error"]=cp.stdout[-800:].replace("\n"," | ")
    print("SHARPREFINE_BG",json.dumps(rec,sort_keys=True),flush=True); rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"sharp_refine_stability.csv",index=False)
stable=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
if stable.empty: raise SystemExit("no stable points")
cur=stable[stable.id=="current"].iloc[0]
print("SHARPREFINE_CURRENT",cur.to_dict(),flush=True)
print("SHARPREFINE_LOWEST_F0",stable.nsmallest(30,"F0").to_dict("records"),flush=True)

# Score a union: 24 smallest F0, 12 largest stability margin among F0<current,
# and current. This avoids selecting an almost-marginal trajectory solely by F0.
keep=set(["current"])
keep.update(stable.nsmallest(32,"F0").id.tolist())
low=stable[stable.F0 < float(cur.F0)]
if not low.empty:
    keep.update(low.nlargest(16,"min_cs2").id.tolist())
sel=stable[stable.id.isin(keep)].copy().sort_values("F0")

scores=[]
for j,r in sel.iterrows():
    tag=str(r.id); af=float(r.AF); zc=float(r.z_c); w=float(r.width)
    root=str(OUT/f"cl_{tag}_"); ip=OUT/f"cl_{tag}.ini"; ip.write_text(ini(af,zc,w,root,True))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=r.to_dict()
    if cp.returncode==0 and Path(root+"00_cl_lensed.dat").exists():
        p=score(root+"00_cl_lensed.dat")
        rec.update(A_planck=p[0],chi2_high=p[1],chi2_lowT=p[2],chi2_lowE=p[3],
                   chi2_lensing=p[4],chi2_cal=p[5],chi2_planck=p[6],delta_planck=p[6]-lscore[-1],
                   cmb_status="OK")
    else:
        rec.update(cmb_status="FAIL",cmb_error=cp.stdout[-1000:].replace("\n"," | "))
    print("SHARPREFINE_CMB",json.dumps(rec,sort_keys=True,default=str),flush=True); scores.append(rec)

sd=pd.DataFrame(scores); sd.to_csv(OUT/"sharp_refine_planck.csv",index=False)
ok=sd[sd.cmb_status=="OK"]
print("SHARPREFINE_BEST_PLANCK",ok.nsmallest(20,"delta_planck").to_dict("records"),flush=True)
print("SHARPREFINE_BEST_LOW_F0",ok.nsmallest(20,"F0").to_dict("records"),flush=True)
