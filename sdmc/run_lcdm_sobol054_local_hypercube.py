#!/usr/bin/env python3
"""
Fair-control six-dimensional ordinary-cosmology screen for flat LCDM.

Uses the exact same deterministic design box and Sobol seed as the SDMC
ordinary-cosmology screen so the two model families can be compared before
the expensive raw DESI/full-Plik promotion.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/lcdm_sobol054_local")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255
CAL_SIGMA=0.0025
OMEGA_R_PHYS=4.17998772e-5

CURRENT=np.array([70.8514,0.02239952,0.12444227328918850,0.964,3.076,0.0544],float)
REFERENCE=np.array([67.36,0.02237,0.1200,0.9649,math.log(1e10*2.10e-9),0.0544],float)
NAMES=["H0","omega_b","omega_cdm","n_s","ln10As","tau_reio"]
MARGIN=np.array([0.55,0.00050,0.0040,0.015,0.040,0.012],float)
LOW=np.minimum(CURRENT,REFERENCE)-MARGIN
HIGH=np.maximum(CURRENT,REFERENCE)+MARGIN

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
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def planck_score(path):
    harr,dls=load_cls(path)
    def pieces(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x)
    return (A,)+pieces(A)

DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:
        a=a[1:]
    if len(a)!=n*n:
        raise RuntimeError(f"covariance mismatch {path}: {len(a)} vs {n*n}")
    return a.reshape(n,n)

def project_setup(C):
    cf=cho_factor(C,lower=True,check_finite=False)
    one=np.ones(C.shape[0]); u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def projected_chi2(pack,resid):
    cf,one,denom=pack
    y=cho_solve(cf,resid,check_finite=False)
    return float(resid.dot(y)-(one.dot(y))**2/denom)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
pp_mag=pp["m_b_corr"].to_numpy(float); pp_z=pp["zHD"].to_numpy(float); pp_zh=pp["zHEL"].to_numpy(float)
pp_C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(pp_mag))
mask=pp_z>0.01
pp_mag,pp_z,pp_zh=pp_mag[mask],pp_z[mask],pp_zh[mask]
pp_C=pp_C[np.ix_(mask,mask)]; pp_pack=project_setup(pp_C)

p=DATA/"Union3/lcparam_full.txt"
lines=p.read_text().splitlines(); cols=lines[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
cmap={c.lower():c for c in un.columns}
un_mag=un[cmap["mb"]].to_numpy(float); un_z=un[cmap["zcmb"]].to_numpy(float); un_zh=un_z.copy()
un_C=read_cov(DATA/"Union3/mag_covmat.txt",len(un_mag)); un_pack=project_setup(un_C)

de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv")
cmap={c.lower():c for c in de.columns}
de_mag=de[cmap["mu"]].to_numpy(float); de_z=de[cmap["zhd"]].to_numpy(float); de_zh=de[cmap["zhel"]].to_numpy(float)
de_err=de[cmap["muerr_final"]].to_numpy(float)
de_C=read_cov(DATA/"DESY5/covsys_000.txt",len(de_mag))+np.diag(de_err*de_err)
de_pack=project_setup(de_C)

def sn_mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1.+z)
    return 5.*np.log10((1.+zh)*(1.+z)*DA)

def sn_scores(bg):
    return (
        projected_chi2(pp_pack,pp_mag-sn_mu(bg,pp_z,pp_zh)),
        projected_chi2(un_pack,un_mag-sn_mu(bg,un_z,un_zh)),
        projected_chi2(de_pack,de_mag-sn_mu(bg,de_z,de_zh)),
    )

def sigma8_from_pk(path):
    a=np.loadtxt(path); k=a[:,0]; p=a[:,1]
    q=(k>1e-4)&(p>0)&np.isfinite(p); k=k[q]; p=p[q]
    x=8*k; W=np.ones_like(x); m=x!=0
    W[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    y=k**3*p*W**2/(2*np.pi**2)
    return float(np.sqrt(np.trapezoid(y,x=np.log(k))))

def ini(x,root):
    H0,ob,oc,ns,lnAs,tau=map(float,x); h=H0/100.; As=math.exp(lnAs)/1e10
    OL=1.-(ob+oc+OMEGA_R_PHYS)/(h*h)
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
    Omega_k = 0
    Omega_Lambda = {OL:.17g}
    modes = s
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
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

SOBOL054=np.array([68.78113517636451,0.022678412813829034,0.11795410220236471,
                         0.9621755753970704,3.0204865057301538,0.0451812363281846],float)
LOCAL_LOW=np.array([68.35,0.02240,0.11720,0.9585,3.008,0.0395],float)
LOCAL_HIGH=np.array([69.15,0.02292,0.11875,0.9655,3.034,0.0510],float)
anchors=[
    ("sobol054_center",SOBOL054),
    ("reference",REFERENCE),
]
sob=qmc.Sobol(d=6,scramble=True,seed=202609212)
xs=qmc.scale(sob.random_base2(m=6),LOCAL_LOW,LOCAL_HIGH)
design=anchors+[(f"local{i+1:03d}",x) for i,x in enumerate(xs)]
print("LCDM_LOCAL_DESIGN_SIZE",len(design),flush=True)
(OUT/"lcdm_sobol054_local_design.json").write_text(json.dumps({
    "names":NAMES,"current":CURRENT.tolist(),"reference":REFERENCE.tolist(),
    "low":LOCAL_LOW.tolist(),"high":LOCAL_HIGH.tolist(),"n_design":len(design)
},indent=2))

rows=[]
for ic,(tag,x) in enumerate(design,1):
    H0,ob,oc,ns,lnAs,tau=map(float,x); h=H0/100.
    root=str(OUT/(f"{ic:03d}_{tag}_")); ip=OUT/(f"{ic:03d}_{tag}.ini"); ip.write_text(ini(x,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"id":tag,"index":ic,"H0":H0,"omega_b":ob,"omega_cdm":oc,"n_s":ns,
         "ln10As":lnAs,"A_s":math.exp(lnAs)/1e10,"tau_reio":tau,
         "Omega_m0":(ob+oc)/(h*h),"returncode":cp.returncode,"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat"); pkp=Path(root+"00_pk.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists() and pkp.exists():
        bg=table(bgp); ps=planck_score(clp); sn=sn_scores(bg)
        rec.update({
            "status":"OK","A_planck":ps[0],
            "chi2_high":ps[1],"chi2_lowT":ps[2],"chi2_lowE":ps[3],
            "chi2_lensing":ps[4],"chi2_cal":ps[5],"chi2_planck":ps[6],
            "chi2_pantheonplus":sn[0],"chi2_union3":sn[1],"chi2_desy5":sn[2],
            "sigma8":sigma8_from_pk(pkp),
        })
    else:
        rec["error"]=cp.stdout[-1400:].replace("\n"," | ")
    print("LCDM_LOCAL_POINT",json.dumps(rec,sort_keys=True),flush=True)
    rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"lcdm_sobol054_local_screen.csv",index=False)
ok=df[df.status.eq("OK")].copy()
if ok.empty: raise SystemExit("No LCDM screen candidate succeeded")
for col,sncol in [("planck_pp","chi2_pantheonplus"),("planck_union3","chi2_union3"),("planck_desy5","chi2_desy5")]:
    ok[col]=ok["chi2_planck"]+ok[sncol]
keep=set(["sobol054_center","reference"])
for col in ["chi2_planck","planck_pp","planck_union3","planck_desy5"]:
    keep.update(ok.nsmallest(8,col)["id"].astype(str))
short=ok[ok.id.astype(str).isin(keep)].copy()
short["screen_best_any"]=short[["chi2_planck","planck_pp","planck_union3","planck_desy5"]].min(axis=1)
short=short.sort_values("screen_best_any")
short.to_csv(OUT/"lcdm_sobol054_local_shortlist.csv",index=False)
print("LCDM_LOCAL_BEST_PLANCK",ok.nsmallest(10,"chi2_planck")[NAMES+["id","chi2_planck","A_planck","sigma8"]].to_dict("records"),flush=True)
print("LCDM_LOCAL_BEST_PP",ok.nsmallest(10,"planck_pp")[NAMES+["id","planck_pp","chi2_planck","chi2_pantheonplus"]].to_dict("records"),flush=True)
print("LCDM_LOCAL_BEST_UNION3",ok.nsmallest(10,"planck_union3")[NAMES+["id","planck_union3","chi2_planck","chi2_union3"]].to_dict("records"),flush=True)
print("LCDM_LOCAL_BEST_DESY5",ok.nsmallest(10,"planck_desy5")[NAMES+["id","planck_desy5","chi2_planck","chi2_desy5"]].to_dict("records"),flush=True)
