#!/usr/bin/env python3
"""
Acoustic-degeneracy screen for the frozen r032 SDMC structural point.

Instead of sampling H0 blindly, each trial solves H0 so that ell_A stays on a
narrow CMB acoustic-scale surface. The sampled coordinates are
(omega_b, omega_cdm, n_s, tau, Q=ln(1e10 As)-2 tau, delta_ellA).

This is a screening likelihood only: Plik-lite + lowT + lowE + lensing, with
Pantheon+, Union3 and DES-Y5 scored separately. Finalists require exact
covariant reconstruction, full Plik nuisance profiling and raw DESI full shape.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, brentq
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/acoustic6d"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
AF=0.0715; ZC=5.00; WIDTH=1.426
D0=0.34231919445927034; POWER=1.0; DFLOOR=0.045
LAMBDA_E=17.925; ZT=17.775; OMEGA_R_PHYS=4.17998772e-5

CURRENT=dict(H0=70.8514,ob=0.02239952,oc=0.12444227328918850,
             ns=0.964,lnAs=3.076,tau=0.0544)
Q0=CURRENT["lnAs"]-2*CURRENT["tau"]

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

def acoustic(bg,th):
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy()
    cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i]); cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross: raise RuntimeError("no drag crossing")
    zd=float(cross[0])
    rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    rd=float(np.interp(zd,bg.z,bg["comov.snd.hrz."]))
    dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return dict(z_star=zstar,z_drag=zd,rs_star=rs,rd=rd,ell_A=math.pi*dm/rs)

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
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),
                        method="bounded",options={"xatol":1e-10})
    A=float(opt.x)
    return (A,)+pieces(A)

DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"covariance size mismatch {path}")
    return a.reshape(n,n)
def psetup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))
def pchi(pack,resid):
    cf,one,den=pack; y=cho_solve(cf,resid,check_finite=False)
    return float(resid.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppmag=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppmag)); m=ppz>0.01
ppmag,ppz,ppzh,C=ppmag[m],ppz[m],ppzh[m],C[np.ix_(m,m)]; pppack=psetup(C)
p=DATA/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unmag=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unzh=unz.copy()
unpack=psetup(read_cov(DATA/"Union3/mag_covmat.txt",len(unmag)))
de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
demag=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
deerr=de[cm["muerr_final"]].to_numpy(float)
depack=psetup(read_cov(DATA/"DESY5/covsys_000.txt",len(demag))+np.diag(deerr*deerr))
def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
    return pchi(pppack,ppmag-mu(bg,ppz,ppzh)),pchi(unpack,unmag-mu(bg,unz,unzh)),pchi(depack,demag-mu(bg,dez,dezh))

def common_ini(H0,ob,oc,ns,lnAs,tau,root,full):
    h=H0/100.; Ox=1.-(ob+oc+OMEGA_R_PHYS)/(h*h); As=math.exp(lnAs)/1e10
    extra = """modes = s
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
    parameters_smg = {AF}, {ZC}, {WIDTH}, {D0:.17g}, {POWER}, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {Ox:.17g}, {LAMBDA_E}, {ZT}, 0.5, 0.01105624999, 0.25, 0.01951933685, 1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
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

# Reference LCDM: sets the acoustic target and score zero.
lr=str(OUT/"lcdm_"); li=OUT/"lcdm.ini"; li.write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(li)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lbg=table(lr+"00_background.dat"); lth=table(lr+"00_thermodynamics.dat")
ELL0=acoustic(lbg,lth)["ell_A"]; lp=planck_score(lr+"00_cl_lensed.dat"); lsn=sns(lbg)
LCDM=dict(planck=lp[-1],pp=lsn[0],union3=lsn[1],desy5=lsn[2],ell_A=ELL0)
print("ACOUSTIC6D_LCDM",LCDM,flush=True)

counter=0
def ell_for(H0,ob,oc,ns,lnAs,tau,tag):
    global counter
    counter+=1
    root=str(OUT/f"root_{tag}_{counter:05d}_"); ip=OUT/f"root_{tag}_{counter:05d}.ini"
    ip.write_text(common_ini(H0,ob,oc,ns,lnAs,tau,root,False))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=100)
    if cp.returncode: return None
    try: return acoustic(table(root+"00_background.dat"),table(root+"00_thermodynamics.dat"))["ell_A"]
    except Exception: return None

def solve_H0(ob,oc,ns,lnAs,tau,target,tag):
    lo,hi=68.0,73.0
    elo=ell_for(lo,ob,oc,ns,lnAs,tau,tag); ehi=ell_for(hi,ob,oc,ns,lnAs,tau,tag)
    if elo is None or ehi is None: return None
    flo=elo-target; fhi=ehi-target
    if flo*fhi>0:
        lo,hi=66.0,75.0
        elo=ell_for(lo,ob,oc,ns,lnAs,tau,tag); ehi=ell_for(hi,ob,oc,ns,lnAs,tau,tag)
        if elo is None or ehi is None or (elo-target)*(ehi-target)>0: return None
    def f(H):
        e=ell_for(H,ob,oc,ns,lnAs,tau,tag)
        if e is None: raise ValueError("root eval failed")
        return e-target
    return float(brentq(f,lo,hi,xtol=2e-4,rtol=2e-7,maxiter=20))

# Coordinates: ob, oc, ns, tau, Q, delta ell_A.
LOW=np.array([0.02210,0.1210,0.957,0.048,Q0-0.025,-0.10])
HIGH=np.array([0.02270,0.1275,0.976,0.070,Q0+0.010,+0.10])
# Current delta ell is measured once, so the exact current point is included on
# the same coordinate system rather than privileged outside the solve.
ecur=ell_for(CURRENT["H0"],CURRENT["ob"],CURRENT["oc"],CURRENT["ns"],CURRENT["lnAs"],CURRENT["tau"],"cur")
dcur=float(ecur-ELL0)
anchors=[
 ("current",np.array([CURRENT["ob"],CURRENT["oc"],CURRENT["ns"],CURRENT["tau"],Q0,dcur])),
 ("tau60_ns968",np.array([CURRENT["ob"],CURRENT["oc"],0.968,0.060,Q0-0.005,dcur])),
 ("tau65_ns970",np.array([CURRENT["ob"],CURRENT["oc"],0.970,0.065,Q0-0.008,dcur])),
 ("tau65_ns972",np.array([CURRENT["ob"],CURRENT["oc"],0.972,0.065,Q0-0.012,dcur])),
 ("tau60_ns972",np.array([CURRENT["ob"],CURRENT["oc"],0.972,0.060,Q0-0.010,dcur])),
]
sob=qmc.Sobol(d=6,scramble=True,seed=20260920)
design=anchors+[(f"sobol{i+1:03d}",x) for i,x in enumerate(qmc.scale(sob.random_base2(m=6),LOW,HIGH))]
rows=[]
for i,(tag,x) in enumerate(design,1):
    ob,oc,ns,tau,Q,dell=map(float,x); lnAs=Q+2*tau; target=ELL0+dell
    H0=solve_H0(ob,oc,ns,lnAs,tau,target,tag)
    rec=dict(id=tag,index=i,omega_b=ob,omega_cdm=oc,n_s=ns,tau_reio=tau,Q=Q,ln10As=lnAs,target_ell_A=target,H0=H0 if H0 else np.nan,status="FAIL")
    if H0 is None:
        rec["error"]="no acoustic root"; print("ACOUSTIC6D_POINT",json.dumps(rec),flush=True); rows.append(rec); continue
    root=str(OUT/f"full_{i:03d}_{tag}_"); ip=OUT/f"full_{i:03d}_{tag}.ini"
    ip.write_text(common_ini(H0,ob,oc,ns,lnAs,tau,root,True))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    if cp.returncode==0 and Path(root+"00_cl_lensed.dat").exists():
        bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat")
        ac=acoustic(bg,th); ps=planck_score(root+"00_cl_lensed.dat"); sn=sns(bg)
        rec.update(ac)
        rec.update(dict(status="OK",A_planck=ps[0],chi2_high=ps[1],chi2_lowT=ps[2],chi2_lowE=ps[3],chi2_lensing=ps[4],chi2_cal=ps[5],chi2_planck=ps[6],
                        delta_planck=ps[6]-LCDM["planck"],delta_pp=sn[0]-LCDM["pp"],delta_union3=sn[1]-LCDM["union3"],delta_desy5=sn[2]-LCDM["desy5"],
                        min_D_all=float(bg["kin (D)"].min()),min_cs2_all=float(bg["c_s^2"].min()),max_cs2_all=float(bg["c_s^2"].max())))
        rec["stable_subluminal"]=bool(rec["min_D_all"]>0 and rec["min_cs2_all"]>0 and rec["max_cs2_all"]<=1)
        rec["joint_pp"]=rec["delta_planck"]+rec["delta_pp"]; rec["joint_union3"]=rec["delta_planck"]+rec["delta_union3"]; rec["joint_desy5"]=rec["delta_planck"]+rec["delta_desy5"]
    else: rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    print("ACOUSTIC6D_POINT",json.dumps(rec,sort_keys=True),flush=True); rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"acoustic6d_screen.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
if ok.empty: raise SystemExit("no stable acoustic candidates")
keep=set(["current","tau60_ns968","tau65_ns970","tau65_ns972","tau60_ns972"])
for col in ["delta_planck","joint_pp","joint_union3","joint_desy5","chi2_lensing","chi2_cal"]:
    keep.update(ok.nsmallest(10,col)["id"].tolist())
short=ok[ok.id.isin(keep)].copy()
short["best_screen"]=short[["delta_planck","joint_pp","joint_union3","joint_desy5"]].min(axis=1)
short=short.sort_values("best_screen"); short.to_csv(OUT/"acoustic6d_shortlist.csv",index=False)
print("ACOUSTIC6D_BEST_PLANCK",ok.nsmallest(12,"delta_planck").to_dict("records"),flush=True)
print("ACOUSTIC6D_BEST_PP",ok.nsmallest(12,"joint_pp").to_dict("records"),flush=True)
print("ACOUSTIC6D_BEST_UNION3",ok.nsmallest(12,"joint_union3").to_dict("records"),flush=True)
print("ACOUSTIC6D_BEST_DESY5",ok.nsmallest(12,"joint_desy5").to_dict("records"),flush=True)
print("ACOUSTIC6D_BEST_LENS",ok.nsmallest(12,"chi2_lensing").to_dict("records"),flush=True)
print("ACOUSTIC6D_BEST_CAL",ok.nsmallest(12,"chi2_cal").to_dict("records"),flush=True)
