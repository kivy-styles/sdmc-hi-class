#!/usr/bin/env python3
"""
Local simultaneous SDMC compensation search.

Vary the ordinary cosmology, primordial amplitude/reionization and physical
No-Slip Planck-mass transition together in one deterministic Sobol design.
The frozen r032 coordinates lambda_e=17.925, z_t=17.775 and Dfloor=0.045 stay
fixed. Every sample must pass D>0, c_s^2>0 and c_s^2<=1. Pantheon+, Union3 and
DES-Y5 are evaluated separately and are never summed together.

This is a screening likelihood (Plik-lite + lowT + lowE + lensing). Finalists
must be replayed with full Plik, raw DESI full shape and exact covariant G_i.
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

OUT=Path("output/simultaneous_physical_local"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025
D0=0.34231919445927034
DFLOOR=.045
LAM=17.925; ZT=17.775
OR=4.17998772e-5

# [H0, omega_b, omega_cdm, n_s, tau, ln10As, A_F, z_c, width]
LOW=np.array([69.55,0.02208,0.1238,0.9610,0.0500,3.050,0.054,2.65,0.48],float)
HIGH=np.array([70.95,0.02247,0.1262,0.9670,0.0605,3.070,0.071,3.55,0.78],float)
ANCHORS=[
 ("old_phys",np.array([70.8514,0.02239952,0.12444227328918850,0.964,0.0544,3.060,
                       0.0661676420550793,3.0664408076554537,0.5638187950477004])),
 ("ac_phys",np.array([69.7095092787552,0.02211778594478965,0.1258596659367904,
                      0.9637720508351921,0.0544,3.060,
                      0.0573794336710125,3.0902165911160413,0.5979282500222326])),
 ("old_original_amp",np.array([70.8514,0.02239952,0.12444227328918850,0.964,0.0544,3.076,
                               0.0661676420550793,3.0664408076554537,0.5638187950477004])),
 ("mid_phys",np.array([70.28,0.0222587,0.1251510,0.963886,0.0544,3.060,
                       0.0618,3.08,0.58])),
]
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

def acoustic(bg,th):
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i]); cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross: return {}
    zd=float(cross[0]); rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    rd=float(np.interp(zd,bg.z,bg["comov.snd.hrz."])); dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return dict(z_star=zstar,z_drag=zd,rs_star=rs,rd=rd,ell_A=math.pi*dm/rs)

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
ppmag=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppmag)); m=ppz>0.01
ppmag,ppz,ppzh,C=ppmag[m],ppz[m],ppzh[m],C[np.ix_(m,m)]; pppack=setup(C)
p=DATA/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unmag=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unzh=unz.copy()
unpack=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(unmag)))
de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
demag=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
deerr=de[cm["muerr_final"]].to_numpy(float)
depack=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(demag))+np.diag(deerr*deerr))
def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
    return pchi(pppack,ppmag-mu(bg,ppz,ppzh)),pchi(unpack,unmag-mu(bg,unz,unzh)),pchi(depack,demag-mu(bg,dez,dezh))

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

def ini(root,v):
    H0,ob,oc,ns,tau,lnAs,AF,ZC,WIDTH=map(float,v)
    h=H0/100.; OX=1.-(ob+oc+OR)/(h*h)
    return textwrap.dedent(f"""\\
    H0 = {H0:.15g}
    omega_b = {ob:.15g}
    omega_cdm = {oc:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {math.exp(lnAs)/1e10:.17e}
    n_s = {ns:.15g}
    tau_reio = {tau:.15g}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF:.17g},{ZC:.17g},{WIDTH:.17g},{D0:.17g},1.0,{DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
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

lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lbg=table(lr+"00_background.dat"); lp=pscore(lr+"00_cl_lensed.dat"); lsn=sns(lbg)
base=dict(planck=lp["chi2_planck"],pp=lsn[0],union3=lsn[1],desy5=lsn[2])
print("PGC_LCDM",base,flush=True)

sob=qmc.Sobol(d=9,scramble=True,seed=20260921)
design=ANCHORS+[(f"sobol{i+1:03d}",v) for i,v in enumerate(qmc.scale(sob.random_base2(m=6),LOW,HIGH))]

rows=[]
for i,(tag,v) in enumerate(design,1):
    H0,ob,oc,ns,tau,lnAs,AF,ZC,WIDTH=map(float,v)
    root=str(OUT/(f"{i:03d}_{tag}_")); ip=OUT/(f"{i:03d}_{tag}.ini"); ip.write_text(ini(root,v))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec=dict(id=tag,index=i,H0=H0,omega_b=ob,omega_cdm=oc,n_s=ns,tau_reio=tau,
             ln10As=lnAs,A_s=math.exp(lnAs)/1e10,A_F=AF,z_c=ZC,width=WIDTH,
             returncode=cp.returncode,status="FAIL")
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()),
                   max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    if cp.returncode==0 and bgp.exists() and thp.exists() and clp.exists():
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.
        rec["stable_subluminal"]=bool(stable)
        if stable:
            bg=table(bgp); th=table(thp); sc=pscore(clp); sn=sns(bg)
            rec.update(sc); rec.update(acoustic(bg,th))
            rec.update(delta_planck=sc["chi2_planck"]-base["planck"],
                       delta_pp=sn[0]-base["pp"],delta_union3=sn[1]-base["union3"],
                       delta_desy5=sn[2]-base["desy5"])
            rec["joint_pp"]=rec["delta_planck"]+rec["delta_pp"]
            rec["joint_union3"]=rec["delta_planck"]+rec["delta_union3"]
            rec["joint_desy5"]=rec["delta_planck"]+rec["delta_desy5"]
            rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    rows.append(rec); print("SPL_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"simultaneous_physical_local.csv",index=False)
ok=df[(df.status=="OK")&(df.stable_subluminal==True)].copy()
if ok.empty: raise SystemExit("no stable simultaneous points")
for col,label in [("delta_planck","PLANCK"),("joint_pp","PP"),("joint_union3","UNION3"),("joint_desy5","DESY5")]:
    print("SPL_BEST_"+label,ok.nsmallest(12,col).to_dict("records"),flush=True)
