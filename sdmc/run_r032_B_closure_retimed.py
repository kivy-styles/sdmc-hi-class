#!/usr/bin/env python3
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import brentq,minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/b_closure"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.8514; OB=0.02239952; OC=0.12444227328918850
OR=4.17998772e-5; Q0=-0.55; TARGET=301.471
LAM=17.925; ZT=17.775; DN=0.5; TAUA=.25; TAUB=1.5
D0=.34231919445927034; DFLOOR=.045
h=H0/100.; OX=1.-(OB+OC+OR)/(h*h); OM=(OB+OC)/(h*h); OR0=OR/(h*h)
A_LATE=1.+Q0-.5*(3.*OM+4.*OR0)

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def acoustic(root):
    bg=tab(Path(root+"00_background.dat")); th=tab(Path(root+"00_thermodynamics.dat"))
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return math.pi*dm/rs,zstar,rs,dm

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
    op=minimize_scalar(lambda a:pcs(float(a))[-1],bounds=(.97,1.03),method="bounded",options={"xatol":1e-10})
    A=float(op.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],chi2_lensing=p[3],chi2_cal=p[4],chi2_total=p[5])

def ini(root,B,AF,ZC,W,ns,lnAs,full=True,skip=False):
    out="" if not full else """
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
    """
    sk="skip_stability_tests_smg = yes" if skip else ""
    return textwrap.dedent(f"""\
    H0 = {H0}
    omega_b = {OB}
    omega_cdm = {OC}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {math.exp(lnAs)/1e10:.17e}
    n_s = {ns}
    tau_reio = 0.0544
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF},{ZC},{W},{D0},1.0,{DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{LAM},{ZT},{DN},{A_LATE:.17g},{TAUA},{B:.17g},{TAUB}
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    {sk}
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
    {out}
    """)

counter=0
def evalB(B):
    global counter
    counter+=1
    root=str(OUT/f"probe_{counter:03d}_")
    ip=OUT/f"probe_{counter:03d}.ini"
    ip.write_text(ini(root,B,.0715,5.,1.426,.964,3.076,False,True))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
    if cp.returncode: raise RuntimeError(cp.stdout[-1200:])
    return acoustic(root)[0]

samples=[]
for B in np.linspace(-.005,.055,25):
    try: samples.append((float(B),evalB(float(B))-TARGET))
    except Exception: pass
br=[]
for (b1,f1),(b2,f2) in zip(samples[:-1],samples[1:]):
    if f1*f2<=0 and f1!=f2: br.append((b1,b2))
if not br: raise SystemExit("no B root "+repr(samples))
b1,b2=min(br,key=lambda q:abs(.5*(q[0]+q[1])-.01951933685))
B=float(brentq(lambda b:evalB(float(b))-TARGET,b1,b2,xtol=1e-8,rtol=1e-10,maxiter=50))
print("B_CLOSURE_SOLVED",json.dumps({"A_late":A_LATE,"B_old":.01951933685,"B_new":B,"target_ellA":TARGET}),flush=True)

# LCDM reference
lcdm=textwrap.dedent("""\
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
root = output/b_closure/lcdm_
format = class
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
spectra_verbose = 0
lensing_verbose = 0
output_verbose = 0
""")
(OUT/"lcdm.ini").write_text(lcdm)
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
base=score(OUT/"lcdm_00_cl_lensed.dat")
print("B_CLOSURE_LCDM",base,flush=True)

cases=[
 ("original_oldB",.01951933685,.0715,5.,1.426,.964,3.076),
 ("original_newB",B,.0715,5.,1.426,.964,3.076),
 ("retimed_prim_oldB",.01951933685,.060,3.5,.8,.962,3.062),
 ("retimed_prim_newB",B,.060,3.5,.8,.962,3.062),
]
rows=[]
for tag,b,af,zc,w,ns,la in cases:
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,b,af,zc,w,ns,la,True,False))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,B=b,A_F=af,z_c=zc,width=w,n_s=ns,ln10As=la,returncode=cp.returncode,status="FAIL")
    if cp.returncode==0:
        bg=tab(Path(root+"00_background.dat")); ea,zs,rs,dm=acoustic(root)
        rec.update(ell_A=ea,rs_star=rs,DM_star=dm,min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
        if rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1:
            sc=score(Path(root+"00_cl_lensed.dat")); rec.update(sc); rec["delta_planck"]=sc["chi2_total"]-base["chi2_total"]; rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    rows.append(rec); print("B_CLOSURE_CASE",json.dumps(rec,sort_keys=True),flush=True)
pd.DataFrame(rows).to_csv(OUT/"b_closure_cases.csv",index=False)
