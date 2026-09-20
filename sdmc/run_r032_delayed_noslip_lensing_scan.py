#!/usr/bin/env python3
"""
Broad delayed-No-Slip diagnostic on the corrected r032 background.

Scan (A_F,z_c,width) well below the historically explored z_c~5 corridor to
test the physical interpretation of the A_L diagnostic: leave F~1 through more
of the CMB-lensing kernel and move the Planck-mass release later.

Two early-tracker closures are evaluated separately (0.003 and cold 1e-8).
This is a diagnostic screen; promoted points require raw DESI and exact
candidate-specific covariant replay.
"""
from pathlib import Path
import json,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/delayed_noslip_local"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; SIG=.0025
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

def cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    h,d=cls(path)
    def f(A):
        A=float(A)
        ch=float(high.chi_squared(h,A_planck=A))
        ct=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(d,**lp)); cp=((A-1)/SIG)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    o=minimize_scalar(lambda x:f(x)[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(o.x); return (A,)+f(A)

def ini(AF,zc,w,cs2,root):
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
    parameters_smg = {AF:.17g},{zc:.17g},{w:.17g},0.34231919445927034,1.0,0.045
    expansion_model = sdmc_full
    expansion_smg = 0.7073985893,17.925,17.775,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
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

anchors=[
 ("current",.0715,5.,1.426),
 ("broad_best",0.05912758882623165,3.284172759018838,0.3490921955090016),
 ("broad_s32",0.06683718590065836,2.8873871725052593,0.8156898508314043),
 ("broad_s16",0.06350275671109556,4.872160436213017,0.4909840957261622)
]
lo=np.array([.050,2.35,.18]); hi=np.array([.0685,4.25,.72])
pts=qmc.scale(qmc.Sobol(d=3,scramble=True,seed=20260922).random_base2(m=6),lo,hi)
design=anchors+[(f"s{i+1:02d}",*v) for i,v in enumerate(pts)]
(OUT/"design.json").write_text(json.dumps({"low":lo.tolist(),"high":hi.tolist(),
 "stage":"local_refine","n_geometry":len(design),"tracker_cs2":[.003,1e-8]},indent=2))

rows=[]; k=0
for cphi in [.003,1e-8]:
  for tag,AF,zc,w in design:
    k+=1; root=str(OUT/(f"r{k:03d}_")); ip=OUT/(f"r{k:03d}.ini")
    ip.write_text(ini(float(AF),float(zc),float(w),cphi,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"id":tag,"tracker_cs2":cphi,"A_F":float(AF),"z_c":float(zc),"width":float(w),
         "returncode":cp.returncode,"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
      bg=table(bgp)
      i0=np.argmin(np.abs(bg.z.to_numpy()))
      rec.update({"min_D_all":float(bg["kin (D)"].min()),"min_cs2_all":float(bg["c_s^2"].min()),
       "max_cs2_all":float(bg["c_s^2"].max()),"F0":float(bg.iloc[i0]["M*^2_smg"]),
       "alphaM0":float(bg.iloc[i0]["M2_running_smg"]),
       "max_abs_noslip_all":float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))})
    stable=rec.get("min_D_all",-1)>0 and rec.get("min_cs2_all",-1)>0 and rec.get("max_cs2_all",2)<=1
    if cp.returncode==0 and clp.exists() and stable:
      v=score(clp)
      rec.update({"status":"OK","A_planck":v[0],"chi2_high":v[1],"chi2_lowT":v[2],
       "chi2_lowE":v[3],"chi2_lensing":v[4],"chi2_cal":v[5],"chi2_total":v[6]})
    else:
      rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    print("DELAYED_NOSLIP_POINT",json.dumps(rec,sort_keys=True),flush=True); rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"delayed_noslip_scan.csv",index=False)
ok=df[df.status.eq("OK")].copy()
for cphi in [.003,1e-8]:
  sub=ok[np.isclose(ok.tracker_cs2,cphi)].copy()
  if sub.empty: continue
  cur=sub[sub.id.eq("current")].iloc[0]
  for c in ["chi2_total","chi2_high","chi2_lensing","chi2_cal"]:
    sub["d_"+c]=sub[c]-float(cur[c])
  sub=sub.sort_values("chi2_total")
  sub.to_csv(OUT/f"delayed_noslip_ranked_{cphi:.0e}.csv",index=False)
  print("DELAYED_NOSLIP_BEST",cphi,sub.head(15).to_dict("records"),flush=True)
