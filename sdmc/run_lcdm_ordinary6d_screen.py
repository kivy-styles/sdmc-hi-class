#!/usr/bin/env python3
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

OUT=Path("output/lcdm6d"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
CURRENT=np.array([70.8514,0.02239952,0.12444227328918850,0.964,3.076,0.0544],float)
REFERENCE=np.array([67.36,0.02237,0.1200,0.9649,math.log(1e10*2.10e-9),0.0544],float)
NAMES=["H0","omega_b","omega_cdm","n_s","ln10As","tau_reio"]
MARGIN=np.array([0.55,0.00050,0.0040,0.015,0.040,0.012],float)
LOW=np.minimum(CURRENT,REFERENCE)-MARGIN; HIGH=np.maximum(CURRENT,REFERENCE)+MARGIN
NUR=3.046
OMEGA_GAMMA=2.47282e-5
OMEGA_R=OMEGA_GAMMA*(1.+0.22710731766*NUR)

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

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def planck_score(path):
    harr,dls=load_cls(path)
    def pieces(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp)); cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda A:pieces(float(A))[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x); return (A,)+pieces(A)

DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"covariance mismatch {path}")
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False); return cf,one,float(one.dot(u))
def pchi(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
pm=pp["m_b_corr"].to_numpy(float); pz=pp["zHD"].to_numpy(float); pzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(pm)); m=pz>0.01
pm,pz,pzh,C=pm[m],pz[m],pzh[m],C[np.ix_(m,m)]; ppack=setup(C)

p=DATA/"Union3/lcparam_full.txt"; lines=p.read_text().splitlines(); cols=lines[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
um=un[cm["mb"]].to_numpy(float); uz=un[cm["zcmb"]].to_numpy(float); uzh=uz.copy()
upack=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(um)))

de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dm=de[cm["mu"]].to_numpy(float); dz=de[cm["zhd"]].to_numpy(float); dzh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dpack=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(dm))+np.diag(derr*derr))

def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1.+z)
    return 5.*np.log10((1.+zh)*(1.+z)*DA)
def sn(bg): return (pchi(ppack,pm-mu(bg,pz,pzh)),pchi(upack,um-mu(bg,uz,uzh)),pchi(dpack,dm-mu(bg,dz,dzh)))

def ini(x,root):
    H0,ob,oc,ns,lnAs,tau=map(float,x); h=H0/100.; Ol=1.-(ob+oc+OMEGA_R)/(h*h)
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {ob:.15g}
    omega_cdm = {oc:.15g}
    N_ncdm = 0
    N_ur = {NUR}
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {math.exp(lnAs)/1e10:.17e}
    n_s = {ns:.15g}
    tau_reio = {tau:.15g}
    Omega_k = 0
    Omega_Lambda = {Ol:.17g}
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

anchors=[("current",CURRENT),("reference_coords",REFERENCE),
         ("line25",CURRENT+.25*(REFERENCE-CURRENT)),("line50",CURRENT+.5*(REFERENCE-CURRENT)),
         ("line75",CURRENT+.75*(REFERENCE-CURRENT))]
sob=qmc.Sobol(d=6,scramble=True,seed=20260920)
xs=qmc.scale(sob.random_base2(m=6),LOW,HIGH)
design=anchors+[(f"sobol{i+1:03d}",x) for i,x in enumerate(xs)]
(OUT/"lcdm6d_design.json").write_text(json.dumps({"names":NAMES,"low":LOW.tolist(),"high":HIGH.tolist(),"n":len(design),
    "omega_r_derived":OMEGA_R},indent=2))

rows=[]
for i,(tag,x) in enumerate(design,1):
    root=str(OUT/(f"{i:03d}_{tag}_")); ip=OUT/(f"{i:03d}_{tag}.ini"); ip.write_text(ini(x,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"id":tag,"index":i,**{n:float(v) for n,v in zip(NAMES,x)},"status":"FAIL","returncode":cp.returncode}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        ps=planck_score(clp); s=sn(table(bgp)); rec.update({
            "status":"OK","A_planck":ps[0],"chi2_planck":ps[-1],
            "chi2_pantheonplus":s[0],"chi2_union3":s[1],"chi2_desy5":s[2],
            "joint_pp":ps[-1]+s[0],"joint_union3":ps[-1]+s[1],"joint_desy5":ps[-1]+s[2]})
    else: rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    print("LCDM6D_POINT",json.dumps(rec,sort_keys=True),flush=True); rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"lcdm6d_screen.csv",index=False)
ok=df[df.status=="OK"].copy()
keep=set(["current","reference_coords","line25","line50","line75"])
for col in ["chi2_planck","joint_pp","joint_union3","joint_desy5"]:
    for x in ok.nsmallest(8,col)["id"]: keep.add(x)
short=ok[ok.id.isin(keep)].copy()
short["best_any"]=short[["chi2_planck","joint_pp","joint_union3","joint_desy5"]].min(axis=1)
short.sort_values("best_any").to_csv(OUT/"lcdm6d_shortlist.csv",index=False)
print("LCDM6D_BEST_PLANCK",ok.nsmallest(10,"chi2_planck")[NAMES+["id","chi2_planck"]].to_dict("records"),flush=True)
print("LCDM6D_BEST_PP",ok.nsmallest(10,"joint_pp")[NAMES+["id","joint_pp"]].to_dict("records"),flush=True)
print("LCDM6D_BEST_UNION3",ok.nsmallest(10,"joint_union3")[NAMES+["id","joint_union3"]].to_dict("records"),flush=True)
print("LCDM6D_BEST_DESY5",ok.nsmallest(10,"joint_desy5")[NAMES+["id","joint_desy5"]].to_dict("records"),flush=True)
