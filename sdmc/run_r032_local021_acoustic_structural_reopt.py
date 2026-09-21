#!/usr/bin/env python3
"""
Acoustically compensated local021-density structural re-optimization.

Goal:
Test the gap-closing route suggested by the equality diagnosis:
  * fix omega_b and omega_m to optimized LCDM local021;
  * allow H0 to move around the acoustic-compensation value ~73;
  * re-optimize only SDMC structure
      (A_F, z_c, DeltaN, D_floor, lambda_e, z_t);
  * keep n_s, A_s and tau at the geometry-refined SDMC values.

This is a Planck screening stage. Any promising point must later pass:
  full Plik, raw DESI, exact-covariant reconstruction, and the SN alternatives.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from scipy.linalg import cho_factor, cho_solve
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/local021_acoustic_structural"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Geometry-refined ordinary/primordial coordinates.
H0G=70.5653567390982
OBG=0.022011983189284802
OCG=0.12404331885203719
NS=0.9625227132590487
TAU=0.055202901571989066
LNAS=3.0598696043919773
AS=math.exp(LNAS)/1e10

# Optimized LCDM local021 matter/baryon densities.
OB=0.0224063693780079
OC=0.1182415105248335

# Geometry-refined structural center.
AF0=0.06017362505197525
ZC0=2.827464461401105
W0=0.5514493708219379
D0=0.34231919445927034
DF0=0.045
LAM0=17.7
ZT0=17.1

# [H0, AF, zc, width, Dfloor, lambda_e, zt]
LOW=np.array([72.25,0.035,2.15,0.38,0.025,16.4,15.8],float)
HIGH=np.array([73.75,0.090,3.55,0.78,0.075,18.8,18.4],float)

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

def crossing(z,y):
    order=np.argsort(z); z=np.asarray(z)[order]; y=np.asarray(y)[order]
    q=y-1.; ii=np.where(q[:-1]*q[1:]<=0)[0]
    if not len(ii): return float("nan")
    i=ii[-1]
    return float(z[i]+(1-y[i])*(z[i+1]-z[i])/(y[i+1]-y[i]))

def derived(bg,th):
    b=bg.sort_values("z")
    matter=b["(.)rho_b"].to_numpy()+b["(.)rho_cdm"].to_numpy()
    rad=b["(.)rho_g"].to_numpy()+b["(.)rho_ur"].to_numpy()
    zeq=crossing(b["z"].to_numpy(),matter/rad)
    imax=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,b["z"],b["comov. dist."]))
    rs=float(np.interp(zstar,b["z"],b["comov.snd.hrz."]))
    return dict(z_eq_exact=zeq,z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=float(np.pi*DM/rs))

# SN setup, only to record geometry penalties for shortlisted screens.
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
    b=bg.sort_values("z")
    DM=np.interp(z,b["z"],b["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
    b=bg.sort_values("z")
    return pchi(ppp,ppm-mu(b,ppz,ppzh)),pchi(unp,unm-mu(b,unz,unz)),pchi(dep,dem-mu(b,dez,dezh))

def ini(root,H0,ob,oc,af,zc,width,dfloor,lam,zt):
    h=H0/100.; ox=1.-(ob+oc+OR)/(h*h)
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
    parameters_smg={af:.17g},{zc:.17g},{width:.17g},{D0:.17g},1.0,{dfloor:.17g}
    expansion_model=sdmc_full
    expansion_smg={ox:.17g},{lam:.17g},{zt:.17g},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg=zero
    method_qs_smg=fully_dynamic
    output_background_smg=3
    modes=s
    output=tCl,pCl,lCl
    lensing=yes
    l_max_scalars=3000
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

def run(tag,H0,ob,oc,af,zc,width,dfloor,lam,zt):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,H0,ob,oc,af,zc,width,dfloor,lam,zt))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=dict(id=tag,H0=H0,omega_b=ob,omega_cdm=oc,omega_m=ob+oc,A_F=af,z_c=zc,width=width,
             D_floor=dfloor,lambda_e=lam,z_t=zt,status="FAIL",returncode=cp.returncode)
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and thp.exists() and clp.exists():
        bg=table(bgp); th=table(thp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            rec.update(pscore(clp)); rec.update(derived(bg,th))
            sn=sns(bg); rec.update(sn_pantheonplus=sn[0],sn_union3=sn[1],sn_desy5=sn[2])
            rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    print("ACCOMP_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

rows=[]
# Geometry-refined baseline and acoustic-compensated center.
geom=run("geom_current",H0G,OBG,OCG,AF0,ZC0,W0,DF0,LAM0,ZT0); rows.append(geom)
center=run("acoustic_center",73.0,OB,OC,AF0,ZC0,W0,DF0,LAM0,ZT0); rows.append(center)

sam=qmc.Sobol(d=7,scramble=False)
pts=qmc.scale(sam.random_base2(m=6),LOW,HIGH)
for i,p in enumerate(pts):
    H0,af,zc,w,df,lam,zt=map(float,p)
    rows.append(run(f"sobol{i:03d}",H0,OB,OC,af,zc,w,df,lam,zt))

df=pd.DataFrame(rows)
base=geom
for rec in rows:
    if rec.get("status")=="OK":
        rec["delta_planck_vs_geom"]=rec["chi2_planck"]-base["chi2_planck"]
        rec["delta_high_vs_geom"]=rec["chi2_high"]-base["chi2_high"]
        rec["delta_lowT_vs_geom"]=rec["chi2_lowT"]-base["chi2_lowT"]
        rec["delta_lensing_vs_geom"]=rec["chi2_lensing"]-base["chi2_lensing"]
        rec["delta_ellA_vs_geom"]=rec["ell_A"]-base["ell_A"]
        rec["delta_pp_vs_geom"]=rec["sn_pantheonplus"]-base["sn_pantheonplus"]
        rec["delta_u3_vs_geom"]=rec["sn_union3"]-base["sn_union3"]
        rec["delta_dy_vs_geom"]=rec["sn_desy5"]-base["sn_desy5"]
df=pd.DataFrame(rows)
df.to_csv(OUT/"local021_acoustic_structural_reopt.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
print("ACCOMP_BEST_PLANCK",ok.nsmallest(12,"chi2_planck").to_dict("records"),flush=True)
print("ACCOMP_BEST_HIGH",ok.nsmallest(12,"chi2_high").to_dict("records"),flush=True)
print("ACCOMP_CLOSEST_ELLA",ok.assign(della=np.abs(ok.ell_A-float(geom["ell_A"]))).nsmallest(12,"della").to_dict("records"),flush=True)
best=ok.nsmallest(1,"chi2_planck").iloc[0].to_dict()
summary={
  "geom_chi2":float(geom["chi2_planck"]),
  "center_chi2":float(center["chi2_planck"]),
  "best":best,
  "best_delta_vs_geom":float(best["chi2_planck"]-geom["chi2_planck"]),
  "geom_ellA":float(geom["ell_A"]),
  "target_z_eq":float(center["z_eq_exact"])
}
(OUT/"local021_acoustic_structural_summary.json").write_text(json.dumps(summary,indent=2))
print("ACCOMP_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
