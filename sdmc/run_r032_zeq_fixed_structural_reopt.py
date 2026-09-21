#!/usr/bin/env python3
"""
Equality-scale constrained structural re-optimization.

Hypothesis under test:
The remaining geometry-refined SDMC primary-CMB residual is driven mainly by
the higher physical matter density / matter-radiation equality scale.

Protocol:
  * Fix omega_b and omega_m to the optimized LCDM local021 values.
  * Therefore omega_cdm is fixed to local021 as well.
  * Keep H0, n_s, tau, A_s at the geometry-refined SDMC values.
  * Re-optimize only the SDMC structural sector:
        (A_F, z_c, DeltaN, D_floor, lambda_e, z_t)
  * Keep D0 fixed at the geometry-refined value and preserve exact No-Slip.
  * Screen with Planck TTTEEE-lite + lowT + lowE + lensing + calibration prior.
  * Record Pantheon+, Union3 and DES-Y5 separately.
Finalists must later pass full Plik, raw DESI and exact-covariant replay.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/zeq_fixed_structural"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Geometry-refined ordinary parameters retained except equality-defining densities.
H0=70.5653567390982
NS=0.9625227132590487
TAU=0.055202901571989066
LNAS=3.0598696043919773
AS=math.exp(LNAS)/1e10

# Original geometry-refined densities (for an internal same-pipeline baseline).
OB_GEOM=0.022011983189284802
OC_GEOM=0.12404331885203719

# Optimized LCDM local021 equality-defining densities.
OB=0.0224063693780079
OC=0.1182415105248335
WM=OB+OC

# Structural geometry-refined center.
AF0=0.06017362505197525
ZC0=2.827464461401105
W0=0.5514493708219379
D0=0.34231919445927034
DF0=0.045
LAM0=17.7
ZT0=17.1

# Broad but targeted structural compensation box.
LOW=np.array([0.035, 2.15, 0.38, 0.025, 16.4, 15.8],float)
HIGH=np.array([0.090, 3.55, 0.78, 0.075, 18.8, 18.4],float)

print("ZEQ_CONVENTION_OMEGA_R",OR,flush=True)
print("ZEQ_GEOM_APPROX", (OB_GEOM+OC_GEOM)/OR-1., flush=True)
print("ZEQ_LCDM021_TARGET_APPROX", WM/OR-1., flush=True)
print("ZEQ_FRACTIONAL_SHIFT", ((OB_GEOM+OC_GEOM)/WM)-1., flush=True)

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

# SN setup
DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"covariance size mismatch {path}: {len(a)} vs {n*n}")
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False); return cf,one,float(one.dot(u))
def pchi(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm)); m=ppz>0.01
ppm,ppz,ppzh,C=ppm[m],ppz[m],ppzh[m],C[np.ix_(m,m)]; ppp=setup(C)

p=DATA/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unp=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1.+z)
    return 5.*np.log10((1.+zh)*(1.+z)*DA)

def sns(bg):
    return pchi(ppp,ppm-mu(bg,ppz,ppzh)),pchi(unp,unm-mu(bg,unz,unz)),pchi(dep,dem-mu(bg,dez,dezh))

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

def ini(root,ob,oc,af,zc,width,dfloor,lam,zt):
    h=H0/100.
    ox=1.-(ob+oc+OR)/(h*h)
    return textwrap.dedent(f"""\
    H0 = {H0}
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
    parameters_smg = {af:.17g},{zc:.17g},{width:.17g},{D0:.17g},1.0,{dfloor:.17g}
    expansion_model = sdmc_full
    expansion_smg = {ox:.17g},{lam:.17g},{zt:.17g},0.5,0.01105624999,0.25,0.01951933685,1.5
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

def run_model(tag,ob,oc,af,zc,width,dfloor,lam,zt):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini")
    ip.write_text(ini(root,ob,oc,af,zc,width,dfloor,lam,zt))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=tag,omega_b=ob,omega_cdm=oc,omega_m=ob+oc,
             A_F=af,z_c=zc,width=width,D_floor=dfloor,lambda_e=lam,z_t=zt,
             status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    if cp.returncode==0 and bgp.exists() and clp.exists():
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            bg=table(bgp); sc=pscore(clp); sn=sns(bg)
            rec.update(sc)
            rec.update(sn_pantheonplus=sn[0],sn_union3=sn[1],sn_desy5=sn[2])
            rec["status"]="OK"
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    return rec

# Fixed reference and current geometry-refined center in the same screen.
lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lbg=table(lr+"00_background.dat"); lp=pscore(lr+"00_cl_lensed.dat"); lsn=sns(lbg)
base=dict(planck=lp["chi2_planck"],pp=lsn[0],u3=lsn[1],dy=lsn[2])
print("ZEQFIX_LCDM_REFERENCE",json.dumps(base,sort_keys=True),flush=True)

geom=run_model("geom_current",OB_GEOM,OC_GEOM,AF0,ZC0,W0,DF0,LAM0,ZT0)
if geom["status"]!="OK": raise RuntimeError("geometry-refined current center failed")
print("ZEQFIX_GEOM_CURRENT",json.dumps(geom,sort_keys=True),flush=True)

eq0=run_model("eqfix_center",OB,OC,AF0,ZC0,W0,DF0,LAM0,ZT0)
print("ZEQFIX_EQUALITY_CENTER",json.dumps(eq0,sort_keys=True),flush=True)

# 64 deterministic Sobol structural points plus both explicit centers.
sam=qmc.Sobol(d=6,scramble=False)
u=sam.random_base2(m=6)
pts=qmc.scale(u,LOW,HIGH)
rows=[geom,eq0]
for i,p in enumerate(pts):
    af,zc,w,df,lam,zt=map(float,p)
    rec=run_model(f"sobol{i:03d}",OB,OC,af,zc,w,df,lam,zt)
    rows.append(rec)
    if rec["status"]=="OK":
        rec["delta_planck"]=rec["chi2_planck"]-base["planck"]
        rec["delta_vs_geom_planck"]=rec["chi2_planck"]-geom["chi2_planck"]
        rec["delta_high_vs_geom"]=rec["chi2_high"]-geom["chi2_high"]
        rec["delta_lowT_vs_geom"]=rec["chi2_lowT"]-geom["chi2_lowT"]
        rec["delta_lensing_vs_geom"]=rec["chi2_lensing"]-geom["chi2_lensing"]
        rec["delta_pp"]=rec["sn_pantheonplus"]-base["pp"]
        rec["delta_union3"]=rec["sn_union3"]-base["u3"]
        rec["delta_desy5"]=rec["sn_desy5"]-base["dy"]
    print("ZEQFIX_POINT",json.dumps(rec,sort_keys=True),flush=True)

# Add deltas for explicit center rows.
for rec in rows[:2]:
    if rec["status"]=="OK":
        rec["delta_planck"]=rec["chi2_planck"]-base["planck"]
        rec["delta_vs_geom_planck"]=rec["chi2_planck"]-geom["chi2_planck"]
        rec["delta_high_vs_geom"]=rec["chi2_high"]-geom["chi2_high"]
        rec["delta_lowT_vs_geom"]=rec["chi2_lowT"]-geom["chi2_lowT"]
        rec["delta_lensing_vs_geom"]=rec["chi2_lensing"]-geom["chi2_lensing"]
        rec["delta_pp"]=rec["sn_pantheonplus"]-base["pp"]
        rec["delta_union3"]=rec["sn_union3"]-base["u3"]
        rec["delta_desy5"]=rec["sn_desy5"]-base["dy"]

df=pd.DataFrame(rows)
df.to_csv(OUT/"zeq_fixed_structural_reopt.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
for col,label in [
    ("chi2_planck","PLANCK"),
    ("chi2_high","HIGH"),
    ("delta_vs_geom_planck","VS_GEOM"),
    ("delta_high_vs_geom","HIGH_VS_GEOM")
]:
    print("ZEQFIX_BEST_"+label,ok.nsmallest(10,col).to_dict("records"),flush=True)

summary={
    "omega_r_convention":OR,
    "z_eq_geom_approx":(OB_GEOM+OC_GEOM)/OR-1.,
    "z_eq_target_approx":WM/OR-1.,
    "geom_screen_chi2":geom["chi2_planck"],
    "eqfix_center_chi2":eq0.get("chi2_planck"),
    "best_id":str(ok.nsmallest(1,"chi2_planck").iloc[0]["id"]),
    "best_chi2_planck":float(ok["chi2_planck"].min()),
    "best_delta_vs_geom":float(ok["chi2_planck"].min()-geom["chi2_planck"]),
}
(OUT/"zeq_fixed_structural_summary.json").write_text(json.dumps(summary,indent=2))
print("ZEQFIX_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
