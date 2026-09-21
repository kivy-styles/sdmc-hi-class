#!/usr/bin/env python3
from pathlib import Path
import json, math, subprocess, textwrap, re
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from iminuit import Minuit
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE import TTTEEE
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/finalordinary_fullplik_line")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255
OR=4.17998772e-5
D0=0.34231919445927034

# Final Q-corrected SDMC structural point.
AF=0.04984375
ZC=2.98625
WIDTH=0.515
DFLOOR=0.041875
LAM=18.40625
ZT=16.496875

# Current final ordinary cosmology.
CUR=dict(
    H0=70.79608815124146,
    ob=0.022083219194622913,
    oc=0.12299536722293603,
    ns=0.9625227132590487,
    tau=0.055202901571989066,
    lnAs=3.054869604391977,
)
CUR["Q"]=CUR["lnAs"]-2*CUR["tau"]

# LCDM local021 target coordinates, used only as a direction in ordinary space.
LOC=dict(
    H0=68.56858744695782,
    ob=0.022406369378007947,
    oc=0.11824151052483357,
    ns=0.9649593164240942,
    tau=0.044985382026527077,
    lnAs=math.log(1e10*2.05109266435559e-9),
)
LOC["Q"]=LOC["lnAs"]-2*LOC["tau"]

ELL0=301.6798271349195
counter=0

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def acoustic(bg,th):
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return dict(z_star=zstar,rs_star=rs,DM_star=dm,ell_A=math.pi*dm/rs)

def ini(H0,ob,oc,ns,lnAs,tau,root,full):
    h=H0/100.
    Ox=1.-(ob+oc+OR)/(h*h)
    As=math.exp(lnAs)/1e10
    extra="""modes = s
output = tCl,pCl,lCl
lensing = yes
l_max_scalars = 3000
""" if full else ""
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {ob:.15g}
    omega_cdm = {oc:.15g}
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
    parameters_smg = {AF:.17g}, {ZC:.17g}, {WIDTH:.17g}, {D0:.17g}, 1.0, {DFLOOR:.17g}
    expansion_model = sdmc_full
    expansion_smg = {Ox:.17g}, {LAM:.17g}, {ZT:.17g}, 0.5, 0.01105624999, 0.25, 0.01951933685, 1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    a_ini_over_a_today_default = 1.e-8
    a_ini_test_qs_smg = 1.e-8
    pert_ic_ini_z_ref_smg = 1.e7
    a_min_stability_test_smg = 1.e-8
    output_background_smg = 3
    {extra}
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

def ell_for(H0,ob,oc,ns,lnAs,tau,tag):
    global counter
    counter+=1
    root=str(OUT/f"lock_{tag}_{counter:04d}_")
    ip=OUT/f"lock_{tag}_{counter:04d}.ini"
    ip.write_text(ini(H0,ob,oc,ns,lnAs,tau,root,False))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=160)
    if cp.returncode: return None
    try:
        return acoustic(table(root+"00_background.dat"),table(root+"00_thermodynamics.dat"))["ell_A"]
    except Exception:
        return None

def solve_H0(ob,oc,ns,lnAs,tau,target,tag):
    for lo,hi in [(67.0,74.0),(65.0,76.0),(63.0,78.0)]:
        elo=ell_for(lo,ob,oc,ns,lnAs,tau,tag)
        ehi=ell_for(hi,ob,oc,ns,lnAs,tau,tag)
        if elo is None or ehi is None: continue
        if (elo-target)*(ehi-target)<=0:
            def f(H):
                e=ell_for(H,ob,oc,ns,lnAs,tau,tag)
                if e is None: raise ValueError("acoustic root eval failed")
                return e-target
            return float(brentq(f,lo,hi,xtol=1e-4,rtol=1e-8,maxiter=32))
    return None

# Ordinary-space line plus acoustic-phase probes.
specs=[
    ("current",0.0,0.0),
    ("f0125",0.125,0.0),
    ("f025",0.25,0.0),
    ("f0375",0.375,0.0),
    ("f050",0.50,0.0),
    ("f0625",0.625,0.0),
    ("f075",0.75,0.0),
    ("f025_ellm03",0.25,-0.03),
    ("f025_ellp03",0.25,+0.03),
    ("f050_ellm03",0.50,-0.03),
    ("f050_ellp03",0.50,+0.03),
]
CANDS=[]
for tag,f,dell in specs:
    ob=CUR["ob"]+f*(LOC["ob"]-CUR["ob"])
    oc=CUR["oc"]+f*(LOC["oc"]-CUR["oc"])
    ns=CUR["ns"]+f*(LOC["ns"]-CUR["ns"])
    tau=CUR["tau"]+f*(LOC["tau"]-CUR["tau"])
    Q=CUR["Q"]+f*(LOC["Q"]-CUR["Q"])
    lnAs=Q+2*tau
    target=ELL0+dell
    H0=solve_H0(ob,oc,ns,lnAs,tau,target,tag)
    if H0 is None:
        print("FINALORDINARY_ACOUSTIC_FAIL",tag,flush=True)
        continue
    CANDS.append(dict(id=tag,fraction=f,dell=dell,H0=H0,ob=ob,oc=oc,ns=ns,tau=tau,Q=Q,lnAs=lnAs,target_ell_A=target))

# Generate spectra and stability metadata.
meta={}
for c in CANDS:
    root=str(OUT/(c["id"]+"_"))
    ip=OUT/(c["id"]+".ini")
    ip.write_text(ini(c["H0"],c["ob"],c["oc"],c["ns"],c["lnAs"],c["tau"],root,True))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    if cp.returncode:
        print("FINALORDINARY_CLASS_FAIL",c["id"],cp.stdout[-2400:],flush=True)
        continue
    bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat")
    ac=acoustic(bg,th)
    rec=dict(c); rec.update(ac)
    rec.update(
        min_D=float(bg["kin (D)"].min()),
        min_cs2=float(bg["c_s^2"].min()),
        max_cs2=float(bg["c_s^2"].max()),
        max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
    )
    rec["stable_subluminal"]=bool(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.)
    print("FINALORDINARY_STABILITY",json.dumps(rec,sort_keys=True),flush=True)
    if rec["stable_subluminal"]:
        meta[c["id"]]=rec

# Standard LCDM spectrum for an absolute reference.
lcdm_root=str(OUT/"lcdm_")
lcdm_ini=OUT/"lcdm.ini"
lcdm_ini.write_text(textwrap.dedent(f"""\
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
output = tCl,pCl,lCl
modes = s
lensing = yes
l_max_scalars = 3000
root = {lcdm_root}
format = class
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
spectra_verbose = 0
lensing_verbose = 0
output_verbose = 0
"""))
cp=subprocess.run(["./class",str(lcdm_ini)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
if cp.returncode: raise RuntimeError(cp.stdout[-2400:])

high=TTTEEE(packages_path="planck_full")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

fixed={
  "cib_index":-1.3,"galf_TE_index":-2.4,"galf_EE_index":-2.4,
  "A_sbpx_100_100_TT":1.0,"A_sbpx_143_143_TT":1.0,"A_sbpx_143_217_TT":1.0,"A_sbpx_217_217_TT":1.0,
  "galf_EE_A_100":0.055,"galf_EE_A_100_143":0.040,"galf_EE_A_100_217":0.094,
  "galf_EE_A_143":0.086,"galf_EE_A_143_217":0.21,"galf_EE_A_217":0.70,
  "A_cnoise_e2e_100_100_EE":1.0,"A_cnoise_e2e_143_143_EE":1.0,"A_cnoise_e2e_217_217_EE":1.0,
  "A_sbpx_100_100_EE":1.0,"A_sbpx_100_143_EE":1.0,"A_sbpx_100_217_EE":1.0,
  "A_sbpx_143_143_EE":1.0,"A_sbpx_143_217_EE":1.0,"A_sbpx_217_217_EE":1.0,
  "A_pol":1.0,"calib_100P":1.021,"calib_143P":0.966,"calib_217P":1.040,
}
names=[
  "A_planck","calib_100T","calib_217T","A_cib_217","xi_sz_cib","A_sz","ksz_norm",
  "gal545_A_100","gal545_A_143","gal545_A_143_217","gal545_A_217",
  "ps_A_100_100","ps_A_143_143","ps_A_143_217","ps_A_217_217",
  "galf_TE_A_100","galf_TE_A_100_143","galf_TE_A_100_217","galf_TE_A_143","galf_TE_A_143_217","galf_TE_A_217"
]
default=np.array([1.0,1.0002,0.99805,67.,0.05,7.,3.,8.6,10.6,23.5,91.9,257.,47.,40.,104.,.130,.130,.46,.207,.69,1.938],float)
seed_final=np.array([
1.0012377970239048,0.999692,0.99825,51.5,0.0,7.38,0.0001,
8.80,11.06,19.5,94.1,261.9,47.3,40.3,117.3,.1153,.1356,.4825,.2298,.6739,2.1288],float)
bounds=[
(.97,1.03),(.9946,1.0058),(.99285,1.00325),(0,200),(0,1),(0,10),(0,10),
(0,24.6),(0,26.6),(0,91.5),(0,251.9),(0,400),(0,400),(0,400),(0,400),
(0,0.466),(0,0.418),(0,1.18),(0,0.783),(0,1.41),(0,6.258)]
gauss={
"A_planck":(1.,.0025),"calib_100T":(1.0002,.0007),"calib_217T":(.99805,.00065),
"gal545_A_100":(8.6,2.),"gal545_A_143":(10.6,2.),"gal545_A_143_217":(23.5,8.5),"gal545_A_217":(91.9,20.),
"galf_TE_A_100":(.130,.042),"galf_TE_A_100_143":(.130,.036),"galf_TE_A_100_217":(.46,.09),
"galf_TE_A_143":(.207,.072),"galf_TE_A_143_217":(.69,.09),"galf_TE_A_217":(1.938,.54)}

if set(high.expected_params)!=(set(fixed)|set(names)):
    raise RuntimeError("full Plik parameter set mismatch")

def load(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    dl={k:np.zeros(n) for k in ["tt","ee","bb","te","pp","tp","ep"]}
    dl["tt"][ell]=a[:,1]*conv; dl["ee"][ell]=a[:,2]*conv; dl["te"][ell]=a[:,3]*conv; dl["bb"][ell]=a[:,4]*conv
    ll=ell.astype(float)*(ell.astype(float)+1.)
    dl["pp"][ell]=a[:,5]*ll
    dl["tp"][ell]=a[:,6]*np.sqrt(ll)*(TCMB*1e6)
    dl["ep"][ell]=a[:,7]*np.sqrt(ll)*(TCMB*1e6)
    cl={k:np.zeros(n) for k in ["tt","ee","bb","te"]}
    fac=np.zeros(n); q=ell>=2; fac[ell[q]]=2*np.pi/ll[q]
    for k in cl: cl[k][ell]=dl[k][ell]*fac[ell]
    return cl,dl

def prior_chi2(p):
    v=sum(((p[n]-mu)/sig)**2 for n,(mu,sig) in gauss.items())
    v+=((p["ksz_norm"]+1.6*p["A_sz"]-9.5)/3.)**2
    return float(v)

def make_objective(path):
    cl,dl=load(path)
    def fun(*x):
        p=dict(fixed); p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p)); A=p["A_planck"]
        lt=float(lowT.log_likelihood(dl["tt"],calib=A))
        le=float(lowE.log_likelihood(dl["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        if not all(np.isfinite([lh,lt,le,ll])): return 1e100
        return float(-2*(lh+lt+le+ll)+prior_chi2(p))
    def pieces(x):
        p=dict(fixed); p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p)); A=p["A_planck"]
        lt=float(lowT.log_likelihood(dl["tt"],calib=A))
        le=float(lowE.log_likelihood(dl["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        return p,lh,lt,le,ll
    return fun,pieces

def profile(label,path,seeds):
    fun,pieces=make_objective(path)
    best=None
    for iseed,seed in enumerate(seeds):
        m=Minuit(fun,*seed,name=names); m.errordef=1.; m.strategy=1; m.tol=.1
        for n,b in zip(names,bounds): m.limits[n]=b
        for n in names: m.errors[n]=max(1e-5,.02*max(abs(m.values[n]),1.))
        m.migrad(ncall=5500)
        if not m.fmin.is_valid:
            m.simplex(ncall=2200); m.migrad(ncall=5500)
        row=(float(m.fval),np.array([m.values[n] for n in names],float),bool(m.fmin.is_valid),int(m.nfcn))
        print("FINALORDINARY_FULLPLIK_SEED",label,iseed,row[0],row[2],row[3],flush=True)
        if best is None or row[0]<best[0]: best=row
    fval,x,valid,nfcn=best
    p,lh,lt,le,ll=pieces(x)
    row=dict(model=label,chi2_profile_total=fval,chi2_high_full=-2*lh,chi2_lowT=-2*lt,
             chi2_lowE=-2*le,chi2_lensing=-2*ll,chi2_nuisance_priors=prior_chi2(p),
             A_planck=p["A_planck"],valid=valid,nfcn=nfcn)
    row.update({n:p[n] for n in names if n!="A_planck"})
    print("FINALORDINARY_FULLPLIK_RESULT",json.dumps(row,sort_keys=True),flush=True)
    return row,x

lcdm,lx=profile("LCDM",OUT/"lcdm_00_cl_lensed.dat",[default,seed_final])
rows=[lcdm]
prev=seed_final
for c in CANDS:
    if c["id"] not in meta: continue
    row,cx=profile(c["id"],OUT/(c["id"]+"_00_cl_lensed.dat"),[prev,seed_final])
    prev=cx
    row.update(meta[c["id"]])
    row["delta_chi2_vs_lcdm"]=row["chi2_profile_total"]-lcdm["chi2_profile_total"]
    rows.append(row)
    print("FINALORDINARY_FULLPLIK_DELTA",c["id"],row["delta_chi2_vs_lcdm"],flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"finalordinary_fullplik_line.csv",index=False)
valid=[r for r in rows[1:] if r.get("valid",False)]
best=min(valid,key=lambda r:r["chi2_profile_total"])
base=next((r for r in valid if r["model"]=="current"),None)
summary=dict(
    lcdm_chi2=lcdm["chi2_profile_total"],
    best_model=best["model"],
    best_delta_vs_lcdm=best["delta_chi2_vs_lcdm"],
    current_delta_vs_lcdm=base["delta_chi2_vs_lcdm"] if base else None,
    improvement_vs_current=(best["chi2_profile_total"]-base["chi2_profile_total"]) if base else None,
)
(OUT/"finalordinary_fullplik_line_summary.json").write_text(json.dumps(summary,indent=2))
print("FINALORDINARY_FULLPLIK_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
