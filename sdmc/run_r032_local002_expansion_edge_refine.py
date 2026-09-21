#!/usr/bin/env python3
"""
Re-optimize the frozen r032 expansion-sector coordinates (lambda_e, z_t)
around the exact local002 cosmology/No-Slip structure.

This tests whether the remaining primary-CMB gap is caused by carrying the
old expansion optimum into a new ordinary cosmology.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/local002_expansion_edge_refine"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5
H0=70.5653567390982; OB=0.022011983189284802; OC=0.12404331885203719
NS=0.9625227132590487; TAU=0.055202901571989066; LNAS=3.0598696043919773
AS=math.exp(LNAS)/1e10
AF=0.06017362505197525; ZC=2.827464461401105; WIDTH=0.5514493708219379
D0=0.34231919445927034; DFLOOR=.045
OX=1.-(OB+OC+OR)/(H0/100.)**2

LAMS=[17.40,17.50,17.60,17.70,17.80]
ZTS =[15.90,16.20,16.50,16.80,17.10]

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
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
    return pchi(ppp,ppm-mu(bg,ppz,ppzh)),pchi(unp,unm-mu(bg,unz,unz)),pchi(dep,dem-mu(bg,dez,dezh))

def lcdm_ini(root):
    return textwrap.dedent(f"""\
    H0=67.36
    omega_b=0.02237
    omega_cdm=0.1200
    N_ncdm=0
    N_ur=3.046
    T_cmb=2.7255
    YHe=0.2453
    A_s=2.10e-9
    n_s=0.9649
    tau_reio=0.0544
    modes=s
    output=tCl,pCl,lCl
    lensing=yes
    l_max_scalars=3000
    write background=yes
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

def ini(root,lam,zt):
    return textwrap.dedent(f"""\
    H0={H0}
    omega_b={OB}
    omega_cdm={OC}
    N_ncdm=0
    N_ur=3.046
    T_cmb=2.7255
    YHe=0.2453
    A_s={AS:.17e}
    n_s={NS}
    tau_reio={TAU}
    Omega_Lambda=0
    Omega_fld=3.1443554e-8
    fluid_equation_of_state=SDMC_TRACKER
    cs2_fld=0.003
    use_ppf=no
    Omega_smg=-1
    gravity_model=sdmc_v3_independent_kinetic
    parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DFLOOR}
    expansion_model=sdmc_full
    expansion_smg={OX:.17g},{lam},{zt},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg=zero
    method_qs_smg=fully_dynamic
    output_background_smg=3
    modes=s
    output=tCl,pCl,lCl
    lensing=yes
    l_max_scalars=3000
    write background=yes
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

lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lbg=table(lr+"00_background.dat"); lp=pscore(lr+"00_cl_lensed.dat"); lsn=sns(lbg)
base=dict(planck=lp["chi2_planck"],pp=lsn[0],u3=lsn[1],dy=lsn[2])

rows=[]
for lam in LAMS:
  for zt in ZTS:
    tag=f"l{lam:.3f}_z{zt:.3f}".replace(".","p")
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,lam,zt))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,lambda_e=lam,z_t=zt,status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
      bg=table(bgp)
      rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                 max_cs2=float(bg["c_s^2"].max()))
    if cp.returncode==0 and bgp.exists() and clp.exists():
      stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
      rec["stable_subluminal"]=bool(stable)
      if stable:
        bg=table(bgp); sc=pscore(clp); sn=sns(bg)
        rec.update(sc)
        rec.update(delta_planck=sc["chi2_planck"]-base["planck"],
                   delta_pp=sn[0]-base["pp"],delta_union3=sn[1]-base["u3"],delta_desy5=sn[2]-base["dy"])
        rec["joint_pp"]=rec["delta_planck"]+rec["delta_pp"]
        rec["joint_union3"]=rec["delta_planck"]+rec["delta_union3"]
        rec["joint_desy5"]=rec["delta_planck"]+rec["delta_desy5"]
        rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-800:].replace("\n"," | ")
    rows.append(rec); print("EXPREF_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"local002_expansion_edge_refine.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
for col,label in [("delta_planck","PLANCK"),("joint_pp","PP"),("joint_union3","UNION3"),("joint_desy5","DESY5")]:
    print("EXPREF_BEST_"+label,ok.nsmallest(8,col).to_dict("records"),flush=True)
