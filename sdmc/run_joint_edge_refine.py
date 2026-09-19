#!/usr/bin/env python3
from pathlib import Path
import csv, math, re, subprocess, textwrap
import numpy as np
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/edge_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
WIDTH=1.426; D0=0.40; POWER=1.0; DFLOOR=1e-4
LCDM_PLANCK=1042.04059473969

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

ZCS=[4.80,4.85,4.90,4.95,5.00,5.05,5.10,5.15,5.20,5.25,5.30,5.35,5.40]
EPS=[1.0e-5,3.0e-5,6.0e-5,1.2e-4]

def ini(AF,zc,root,full):
    As=math.exp(3.076)/1e10
    outputs="""modes = s
output = tCl,pCl,lCl
lensing = yes
l_max_scalars = 3000
write background = yes
""" if full else """write background = yes
"""
    return textwrap.dedent(f"""\
    H0 = 70.8514
    omega_b = 0.02239952
    omega_cdm = 0.12444227328918850
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {As:.16e}
    n_s = 0.964
    tau_reio = 0.0544
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF:.12g}, {zc:.12g}, {WIDTH:.12g}, {D0}, {POWER}, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = 0.7073985893,20.0,21.10,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    {outputs}
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

def run(AF,zc,tag,full=False):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini")
    ip.write_text(ini(AF,zc,root,full))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240 if full else 60)
    (OUT/(tag+".log")).write_text(cp.stdout)
    return cp,root

def parse_bg(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
    harr,dls=load_cls(path)
    def comp(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda A:comp(A)[-1],bounds=(.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x)
    return (A,)+comp(A)

def succeeds(AF,zc,tag):
    cp,root=run(AF,zc,tag,False)
    return cp.returncode==0

thresholds=[]
rows=[]
for iz,zc in enumerate(ZCS):
    lo,hi=0.067,0.078
    # Guarantee monotonic bracket: low unstable, high stable.
    if succeeds(lo,zc,f"br_{iz}_lo"):
        raise RuntimeError(f"lower bracket unexpectedly stable at zc={zc}")
    if not succeeds(hi,zc,f"br_{iz}_hi"):
        raise RuntimeError(f"upper bracket unexpectedly unstable at zc={zc}")
    for k in range(13):
        mid=.5*(lo+hi)
        if succeeds(mid,zc,f"br_{iz}_{k}"):
            hi=mid
        else:
            lo=mid
    threshold=hi
    thresholds.append(dict(zc=zc,AF_threshold=threshold,AF_unstable=lo,
                           bracket_width=hi-lo))
    print("EDGE_THRESHOLD",zc,threshold,lo,hi-lo,flush=True)

    for je,eps in enumerate(EPS):
        AF=threshold+eps
        tag=f"e_{iz:02d}_{je:02d}"
        cp,root=run(AF,zc,tag,True)
        rec={"zc":zc,"AF":AF,"AF_threshold":threshold,"epsilon":eps,
             "width":WIDTH,"returncode":cp.returncode,"status":"FAIL"}
        bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
        if bgp.exists():
            d=parse_bg(bgp); sub=d["z"]<=100
            rec.update({
                "min_D_z100":float(np.min(d["kin (D)"][sub])),
                "min_cs2_z100":float(np.min(d["c_s^2"][sub])),
                "max_cs2_z100":float(np.max(d["c_s^2"][sub])),
                "min_cs2_all":float(np.min(d["c_s^2"])),
                "max_abs_noslip_z100":float(np.max(np.abs(d["braiding_smg"][sub]+2*d["M2_running_smg"][sub]))),
            })
        if cp.returncode==0 and bgp.exists() and clp.exists():
            sc=score(clp)
            rec.update({"status":"OK","A_planck":sc[0],"chi2_high":sc[1],
                        "chi2_lowT":sc[2],"chi2_lowE":sc[3],"chi2_lensing":sc[4],
                        "chi2_cal":sc[5],"chi2_planck":sc[6],
                        "delta_planck":sc[6]-LCDM_PLANCK})
            rec["stable_subluminal"]=bool(rec["min_D_z100"]>0 and rec["min_cs2_z100"]>0 and rec["max_cs2_z100"]<=1)
            print("EDGE_OK",rec,flush=True)
        else:
            rec["error"]=cp.stdout[-1000:]
            print("EDGE_FAIL",rec,flush=True)
        rows.append(rec)

pd_fields=sorted({k for r in rows for k in r})
with (OUT/"edge_refine.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=pd_fields); w.writeheader(); w.writerows(rows)
with (OUT/"edge_thresholds.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(thresholds[0])); w.writeheader(); w.writerows(thresholds)
ok=[r for r in rows if r.get("status")=="OK" and r.get("stable_subluminal")]
ok=sorted(ok,key=lambda r:r["chi2_planck"])
with (OUT/"edge_best.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=pd_fields); w.writeheader()
    for r in ok[:15]: w.writerow(r)
print("EDGE_BEST",ok[:15],flush=True)
