#!/usr/bin/env python3
"""
Self-consistent ordinary-cosmology screen for SDMC r032.

Unlike the earlier frozen-A/B six-dimensional screen, this enforces the
construction conditions at every trial point:
  q0 = -0.55  -> derived A_late,
  ell_A = 301.471 -> B_late solved numerically.

The structural r032 point (lambda_e,z_t,Dfloor) and the current No-Slip
gravity trajectory are held fixed in this stage.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/selfconsistent6d"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
Q0=-0.55; TARGET_ELLA=301.471
LAMBDA=17.925; ZT=17.775; DN=0.5; TAUA=0.25; TAUB=1.5
AF=0.0715; ZC=5.; WIDTH=1.426; D0=0.34231919445927034; DFLOOR=0.045
OMEGA_R_PHYS=4.17998772e-5
CURRENT=np.array([70.8514,0.02239952,0.12444227328918850,0.964,3.076,0.0544],float)
REFERENCE=np.array([67.36,0.02237,0.1200,0.9649,math.log(1e10*2.10e-9),0.0544],float)
NAMES=["H0","omega_b","omega_cdm","n_s","ln10As","tau_reio"]
MARGIN=np.array([0.55,0.00050,0.0040,0.015,0.040,0.012],float)
LOW=np.minimum(CURRENT,REFERENCE)-MARGIN
HIGH=np.maximum(CURRENT,REFERENCE)+MARGIN

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

def acoustic(root):
    bg=table(Path(root+"00_background.dat")); th=table(Path(root+"00_thermodynamics.dat"))
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

def planck_score(path):
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
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

# SN likelihoods
DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"bad covariance {path}")
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False); return cf,one,float(one.dot(u))
def proj(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
pm=pp["m_b_corr"].to_numpy(float); pz=pp["zHD"].to_numpy(float); pzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(pm)); m=pz>0.01
pm,pz,pzh,C=pm[m],pz[m],pzh[m],C[np.ix_(m,m)]; ppack=setup(C)
p=DATA/"Union3/lcparam_full.txt"; lines=p.read_text().splitlines(); cols=lines[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
um=un[cm["mb"]].to_numpy(float); uz=un[cm["zcmb"]].to_numpy(float); upack=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(um)))
de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dm=de[cm["mu"]].to_numpy(float); dz=de[cm["zhd"]].to_numpy(float); dzh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dpack=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(dm))+np.diag(derr*derr))
def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1.+z)
    return 5*np.log10((1.+zh)*(1.+z)*DA)
def sn_scores(bg):
    return proj(ppack,pm-mu(bg,pz,pzh)),proj(upack,um-mu(bg,uz,uz)),proj(dpack,dm-mu(bg,dz,dzh))

def derive_A(H0,ob,oc):
    h=H0/100.; Om=(ob+oc)/(h*h); Or=OMEGA_R_PHYS/(h*h)
    return 1.+Q0-0.5*(3.*Om+4.*Or)

def ini_text(x,root,A_late,B_late,full,skip_stability=False):
    H0,ob,oc,ns,lnAs,tau=map(float,x); h=H0/100.
    Ox=1.-(ob+oc+OMEGA_R_PHYS)/(h*h)
    sskip="skip_stability_tests_smg = yes" if skip_stability else ""
    outblk="" if not full else """
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
    """
    return textwrap.dedent(f"""\
    H0 = {H0:.16g}
    omega_b = {ob:.16g}
    omega_cdm = {oc:.16g}
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
    expansion_smg = {Ox:.17g},{LAMBDA},{ZT},{DN},{A_late:.17g},{TAUA},{B_late:.17g},{TAUB}
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    {sskip}
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
    {outblk}
    """)

probe_counter=0
def eval_ella(x,A_late,B,key):
    global probe_counter
    probe_counter+=1
    root=str(OUT/(f"probe_{key}_{probe_counter:05d}_"))
    ip=OUT/(f"probe_{key}_{probe_counter:05d}.ini")
    ip.write_text(ini_text(x,root,A_late,B,False,True))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=120)
    if cp.returncode: raise RuntimeError(cp.stdout[-1200:])
    return acoustic(root)[0]

def solve_B(x,A_late,key):
    vals=[]
    for B in np.linspace(-0.005,0.055,13):
        try: vals.append((float(B),eval_ella(x,A_late,float(B),key)-TARGET_ELLA))
        except Exception: pass
    brackets=[]
    for (b1,f1),(b2,f2) in zip(vals[:-1],vals[1:]):
        if f1==0: return b1
        if f1*f2<=0: brackets.append((b1,b2))
    if not brackets:
        raise RuntimeError(f"no B bracket; samples={vals}")
    b1,b2=min(brackets,key=lambda q:abs(0.5*(q[0]+q[1])-0.0195))
    return float(brentq(lambda b:eval_ella(x,A_late,float(b),key)-TARGET_ELLA,b1,b2,xtol=2e-7,rtol=1e-9,maxiter=32))

# Fixed LCDM reference for deltas.
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
lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lps=planck_score(Path(lr+"00_cl_lensed.dat")); lbg=table(Path(lr+"00_background.dat")); lsn=sn_scores(lbg)
BASE=dict(planck=lps["chi2_planck"],pp=lsn[0],union3=lsn[1],desy5=lsn[2])
print("SELF6D_LCDM",BASE,flush=True)

anchors=[("current",CURRENT),("reference_coords",REFERENCE),
         ("line25",CURRENT+.25*(REFERENCE-CURRENT)),("line50",CURRENT+.5*(REFERENCE-CURRENT)),
         ("line75",CURRENT+.75*(REFERENCE-CURRENT))]
sob=qmc.Sobol(d=6,scramble=True,seed=20260920)
xs=qmc.scale(sob.random_base2(m=5),LOW,HIGH)
design=anchors+[(f"sobol{i+1:03d}",x) for i,x in enumerate(xs)]

rows=[]
for i,(tag,x) in enumerate(design,1):
    rec={"id":tag,"index":i,**{n:float(v) for n,v in zip(NAMES,x)},"status":"FAIL"}
    A_late=derive_A(*x[:3]); rec["A_late"]=A_late
    try:
        B_late=solve_B(x,A_late,tag); rec["B_late"]=B_late
        root=str(OUT/(f"{i:03d}_{tag}_")); ip=OUT/(f"{i:03d}_{tag}.ini")
        ip.write_text(ini_text(x,root,A_late,B_late,True,False))
        cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
        rec["returncode"]=cp.returncode
        if cp.returncode: raise RuntimeError(cp.stdout[-1400:])
        bg=table(Path(root+"00_background.dat")); ea,zs,rs,dmstar=acoustic(root)
        rec.update(ell_A=ea,z_star=zs,rs_star=rs,DM_star=dmstar,
                   min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            ps=planck_score(Path(root+"00_cl_lensed.dat")); sn=sn_scores(bg); rec.update(ps)
            rec.update(delta_planck=ps["chi2_planck"]-BASE["planck"],
                       delta_pp=sn[0]-BASE["pp"],delta_union3=sn[1]-BASE["union3"],delta_desy5=sn[2]-BASE["desy5"])
            rec["joint_pp"]=rec["delta_planck"]+rec["delta_pp"]
            rec["joint_union3"]=rec["delta_planck"]+rec["delta_union3"]
            rec["joint_desy5"]=rec["delta_planck"]+rec["delta_desy5"]
            rec["status"]="OK"
    except Exception as ex:
        rec["error"]=str(ex)
    rows.append(rec); print("SELF6D_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"selfconsistent_ordinary6d.csv",index=False)
ok=df[df.status=="OK"].copy()
if ok.empty: raise SystemExit("no successful self-consistent points")
print("SELF6D_BEST_PLANCK",ok.nsmallest(15,"delta_planck").to_dict("records"),flush=True)
print("SELF6D_BEST_PP",ok.nsmallest(15,"joint_pp").to_dict("records"),flush=True)
print("SELF6D_BEST_UNION3",ok.nsmallest(15,"joint_union3").to_dict("records"),flush=True)
print("SELF6D_BEST_DESY5",ok.nsmallest(15,"joint_desy5").to_dict("records"),flush=True)
