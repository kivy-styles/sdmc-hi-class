#!/usr/bin/env python3
"""
Acoustically matched continuation along the SDMC -> LCDM local021 density axis.

f=0  : geometry-refined SDMC densities
f=1  : optimized LCDM local021 omega_b and omega_cdm
f<0  : continuation in the opposite (higher omega_m) direction

At each f, solve H0 so that the SDMC acoustic scale ell_A matches the
geometry-refined target. This removes peak-position drift and isolates the
effect of matter-radiation equality / baryon loading on the CMB shape.

SDMC structural and primordial coordinates are held at the current
geometry-refined expansion-best values.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import brentq, minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/zeq_acoustic_continuation"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

H0G=70.5653567390982
OBG=0.022011983189284802
OCG=0.12404331885203719
OBL=0.0224063693780079
OCL=0.1182415105248335

NS=0.9625227132590487
TAU=0.055202901571989066
LNAS=3.0598696043919773
AS=math.exp(LNAS)/1e10

AF=0.06017362505197525
ZC=2.827464461401105
WIDTH=0.5514493708219379
D0=0.34231919445927034
DF=0.045
LAM=17.7
ZT=17.1

FS=np.array([-0.50,-0.40,-0.30,-0.20,-0.10,0.0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.0])

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

def derived(bg,th):
    b=bg.sort_values("z")
    matter=b["(.)rho_b"].to_numpy()+b["(.)rho_cdm"].to_numpy()
    rad=b["(.)rho_g"].to_numpy()+b["(.)rho_ur"].to_numpy()
    z=b["z"].to_numpy(); ratio=matter/rad; q=ratio-1.
    ii=np.where(q[:-1]*q[1:]<=0)[0]
    zeq=float("nan")
    if len(ii):
        i=ii[-1]
        zeq=float(z[i]+(1-ratio[i])*(z[i+1]-z[i])/(ratio[i+1]-ratio[i]))
    imax=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,b["z"],b["comov. dist."]))
    rs=float(np.interp(zstar,b["z"],b["comov.snd.hrz."]))
    return dict(z_eq_exact=zeq,z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=float(np.pi*DM/rs))

def ini(root,H0,ob,oc,full):
    h=H0/100.; ox=1.-(ob+oc+OR)/(h*h)
    extra = """
    modes=s
    output=tCl,pCl,lCl
    lensing=yes
    l_max_scalars=3000
    """ if full else ""
    return textwrap.dedent(f"""\
    H0={H0:.17g}
    omega_b={ob:.17g}
    omega_cdm={oc:.17g}
    N_ncdm=0
    N_ur=3.046
    T_cmb=2.7255
    YHe=0.2453
    A_s={AS:.17e}
    n_s={NS:.17g}
    tau_reio={TAU:.17g}
    Omega_Lambda=0
    Omega_fld=3.1443554e-8
    fluid_equation_of_state=SDMC_TRACKER
    cs2_fld=0.003
    use_ppf=no
    Omega_smg=-1
    gravity_model=sdmc_v3_independent_kinetic
    parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
    expansion_model=sdmc_full
    expansion_smg={ox:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg=zero
    method_qs_smg=fully_dynamic
    output_background_smg=3
    {extra}
    write background=yes
    write thermodynamics=yes
    root={root}
    format=class
    input_verbose=0
    background_verbose=0
    thermodynamics_verbose=0
    perturbations_verbose=0
    spectra_verbose=0
    lensing_verbose=0
    output_verbose=0
    """)

def run_bg(tag,H0,ob,oc,full=False):
    root=str(OUT/(tag+"_"))
    ip=OUT/(tag+".ini"); ip.write_text(ini(root,H0,ob,oc,full))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    if cp.returncode:
        raise RuntimeError(f"{tag} failed: {cp.stdout[-1200:]}")
    bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat")
    der=derived(bg,th)
    stab=dict(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
    return root,der,stab

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
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

# Baseline target acoustic scale.
root0,der0,stab0=run_bg("baseline_target",H0G,OBG,OCG,full=True)
TARGET=der0["ell_A"]
base_score=pscore(root0+"00_cl_lensed.dat")
print("ZEQCONT_TARGET",json.dumps({**der0,**stab0,**base_score},sort_keys=True),flush=True)

rows=[]
for idx,f in enumerate(FS):
    ob=OBG+f*(OBL-OBG)
    oc=OCG+f*(OCL-OCG)

    cache={}
    def g(H0):
        key=round(float(H0),8)
        if key not in cache:
            _,d,_=run_bg(f"root_f{idx:02d}_{str(key).replace('.','p')}",float(H0),ob,oc,full=False)
            cache[key]=d["ell_A"]-TARGET
        return cache[key]

    # Wide bracket; ell_A decreases monotonically with H0 in this neighborhood.
    lo,hi=67.0,77.0
    flo,fhi=g(lo),g(hi)
    if flo*fhi>0:
        raise RuntimeError(f"no acoustic root for f={f}: g(lo)={flo}, g(hi)={fhi}")
    H0=float(brentq(g,lo,hi,xtol=2e-5,rtol=1e-10,maxiter=30))

    root,der,stab=run_bg(f"final_f{idx:02d}",H0,ob,oc,full=True)
    score=pscore(root+"00_cl_lensed.dat")
    rec=dict(f=float(f),H0=H0,omega_b=ob,omega_cdm=oc,omega_m=ob+oc,
             **der,**stab,**score,
             delta_planck_vs_geom=score["chi2_planck"]-base_score["chi2_planck"],
             delta_high_vs_geom=score["chi2_high"]-base_score["chi2_high"],
             delta_lowT_vs_geom=score["chi2_lowT"]-base_score["chi2_lowT"],
             delta_lensing_vs_geom=score["chi2_lensing"]-base_score["chi2_lensing"],
             delta_ellA=der["ell_A"]-TARGET)
    rows.append(rec)
    print("ZEQCONT_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"zeq_acoustic_continuation.csv",index=False)
best=df.nsmallest(5,"chi2_planck").to_dict("records")
summary={"target_ellA":TARGET,"best":best,"f_min":float(df.loc[df.chi2_planck.idxmin(),"f"]),
         "best_delta_vs_geom":float(df.delta_planck_vs_geom.min())}
(OUT/"zeq_acoustic_continuation_summary.json").write_text(json.dumps(summary,indent=2))
print("ZEQCONT_BEST",json.dumps(best,sort_keys=True),flush=True)
print("ZEQCONT_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
