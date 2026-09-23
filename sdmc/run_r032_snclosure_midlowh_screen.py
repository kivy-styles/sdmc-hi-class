#!/usr/bin/env python3
"""
Intermediate low-H SN-aware closure screen around the first direct-SN-closure region.

Stage 1: 256 Sobol backgrounds in (A_late, B_late, H0), exact SN covariances,
         acoustic-scale proxy, DESI DR1 Gaussian BAO, exact stability columns.
Stage 2: full C_l + Planck native-lite likelihood on the best 24 backgrounds.

The target is not "best SN alone". It is:
  (i) retain a plausible Planck+DESI advantage over optimized LCDM local021;
  (ii) make at least two of Pantheon+, Union3, DES-Y5 cross local021 after
       adding that same Planck+DESI advantage.

This is a screen. Shortlist winners must still be promoted to exact covariant
full-Plik + raw DESI full shape + all three SN likelihoods.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap, shutil
import numpy as np
import pandas as pd
from scipy.stats import qmc
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve

from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure_midlowh"); OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"tmp"; TMP.mkdir(exist_ok=True)
SNROOT=Path("sn_data")
BAOROOT=Path("bao_data")

# Successful edge024-Q structural / primordial anchor.
OB=0.022083219194622913
OC=0.12299536722293603
AS=2.119721074504234e-9
NS=0.9625227132590487
TAU=0.055202901571989066
AF=0.024290704212870225
ZC=3.6096407580714724
WIDTH=0.3631213944496205
D0=0.34231919445927034
DF=0.04602317962943721
LAM=18.40625
ZT=16.173189924377947
DNT=0.5
TAUA=0.25
TAUB=1.5
OR=4.17998772e-5

A0=0.00792232021316886
B0=0.012283301661722363
H00=69.71275747716427
EDGE024_PD_FAIR=-0.21524906298265822
EDGE024_BAO_CHI=13.251270395581628
LOCAL_BAO_CHI=13.400290718544086
LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}
RD_FIXED=145.8653484988842
ZSTAR=1089.9
TCMB=2.7255; CAL_SIGMA=.0025

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov mismatch {path}: {len(a)} != {n*n}")
    return a.reshape(n,n)

def setup_cov(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def chi_profile(pack,r):
    cf,one,den=pack
    y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

# SN covariance products
pp=pd.read_csv(SNROOT/'PantheonPlus/Pantheon+SH0ES.dat',sep=r'\s+',header=0,engine='python')
ppm=pp['m_b_corr'].to_numpy(float); ppz=pp['zHD'].to_numpy(float); ppzh=pp['zHEL'].to_numpy(float)
Cpp=read_cov(SNROOT/'PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov',len(ppm))
m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]
pppack=setup_cov(Cpp)

p=SNROOT/'Union3/lcparam_full.txt'
cols=p.read_text().splitlines()[0].removeprefix('#').split()
un=pd.read_csv(p,sep=r'\s+',comment='#',names=cols,engine='python')
cm={c.lower():c for c in un.columns}
unm=un[cm['mb']].to_numpy(float); unz=un[cm['zcmb']].to_numpy(float)
unpack=setup_cov(read_cov(SNROOT/'Union3/mag_covmat.txt',len(unm)))

de=pd.read_csv(SNROOT/'DESY5/DES-SN5YR_HD.csv')
cm={c.lower():c for c in de.columns}
dem=de[cm['mu']].to_numpy(float); dez=de[cm['zhd']].to_numpy(float); dezh=de[cm['zhel']].to_numpy(float)
derr=de[cm['muerr_final']].to_numpy(float)
depack=setup_cov(read_cov(SNROOT/'DESY5/covsys_000.txt',len(dem))+np.diag(derr*derr))

# DESI DR1 BAO
bao_rows=[]
for ln in (BAOROOT/'desi_2024_gaussian_bao_ALL_GCcomb_mean.txt').read_text().splitlines():
    if not ln.strip() or ln.startswith('#'): continue
    z,v,q=ln.split(); bao_rows.append((float(z),float(v),q))
bao_obs=np.array([x[1] for x in bao_rows])
bao_cov=np.loadtxt(BAOROOT/'desi_2024_gaussian_bao_ALL_GCcomb_cov.txt')
bao_inv=np.linalg.inv(bao_cov)

def interp(bg,col,z):
    zz=bg.z.to_numpy(); yy=bg[col].to_numpy(); o=np.argsort(zz)
    return float(np.interp(z,zz[o],yy[o]))

def mu_bg(bg,z,zh):
    zz=bg.z.to_numpy(); dm=bg['comov. dist.'].to_numpy(); o=np.argsort(zz)
    DM=np.interp(z,zz[o],dm[o]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def sn_scores(bg):
    return {
      "pantheonplus":chi_profile(pppack,ppm-mu_bg(bg,ppz,ppzh)),
      "union3":chi_profile(unpack,unm-mu_bg(bg,unz,unz)),
      "desy5":chi_profile(depack,dem-mu_bg(bg,dez,dezh)),
    }

def bao_score(bg):
    pred=[]
    for z,_,q in bao_rows:
        DM=interp(bg,'comov. dist.',z)
        H=interp(bg,'H [1/Mpc]',z)
        DH=1/H; DV=(z*DM*DM*DH)**(1/3)
        pred.append({"DM_over_rs":DM/RD_FIXED,
                     "DH_over_rs":DH/RD_FIXED,
                     "DV_over_rs":DV/RD_FIXED}[q])
    d=np.array(pred)-bao_obs
    return float(d@bao_inv@d)

def ini_text(root,A,B,H0,cls=False):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    output = """modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
""" if cls else ""
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
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
expansion_smg={ox:.17g},{LAM},{ZT},{DNT},{A:.17g},{TAUA},{B:.17g},{TAUB}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
{output}write background=yes
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

def run_background(label,A,B,H0):
    root=str(TMP/(label+"_")); ip=TMP/(label+".ini")
    ip.write_text(ini_text(root,A,B,H0,False))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=120)
    rec={"id":label,"A":float(A),"B":float(B),"H0":float(H0),
         "status":"FAIL","returncode":cp.returncode}
    p=Path(root+"00_background.dat")
    if cp.returncode==0 and p.exists():
        bg=table(p)
        try:
            rec["min_D"]=float(bg["kin (D)"].min())
            rec["min_cs2"]=float(bg["c_s^2"].min())
            rec["max_cs2"]=float(bg["c_s^2"].max())
            rec["max_abs_noslip"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
            stable=(rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1)
            rec["stable_subluminal"]=bool(stable)
            if stable:
                sn=sn_scores(bg)
                rec.update({f"sn_{k}":v for k,v in sn.items()})
                for k,v in sn.items(): rec[f"sn_delta_{k}"]=v-LOCAL_SN[k]
                rec["bao_chi2"]=bao_score(bg)
                rec["bao_delta_local"]=rec["bao_chi2"]-LOCAL_BAO_CHI
                DM=interp(bg,'comov. dist.',ZSTAR)
                rs=interp(bg,'comov.snd.hrz.',ZSTAR)
                rec["ellA_proxy"]=math.pi*DM/rs
                rec["status"]="OK"
        except Exception as e:
            rec["error"]=repr(e)
    if rec["status"]!="OK":
        rec["error"]=rec.get("error",cp.stdout[-800:].replace("\n"," | "))
    # clean per-point files
    for q in TMP.glob(label+"_*"): 
        try: q.unlink()
        except: pass
    try: ip.unlink()
    except: pass
    return rec

# Center establishes the acoustic proxy exactly in this screening code.
center=run_background("center",A0,B0,H00)
if center["status"]!="OK": raise RuntimeError(f"center failed: {center}")
ELL0=center["ellA_proxy"]
print("SNCLOSE_CENTER_BACKGROUND",json.dumps(center,sort_keys=True),flush=True)

sob=qmc.Sobol(d=3,scramble=True,seed=6006107)
u=sob.random_base2(m=8)  # 256
# Focused extension through the lower-H boundary exposed by sobol006.
Alo,Ahi=0.0035,0.0105
Blo,Bhi=0.0100,0.0145
Hlo,Hhi=69.20,69.76
rows=[]
for i,x in enumerate(u):
    A=Alo+(Ahi-Alo)*x[0]
    B=Blo+(Bhi-Blo)*x[1]
    H0=Hlo+(Hhi-Hlo)*x[2]
    r=run_background(f"sobol{i:03d}",A,B,H0)
    if r["status"]=="OK":
        r["ellA_frac"]=(r["ellA_proxy"]/ELL0-1.0)
        # Direct target: at least two exact SN covariance scores below matched LCDM.
        sn2=sorted([r["sn_delta_pantheonplus"],r["sn_delta_union3"],r["sn_delta_desy5"]])[1]
        r["second_best_sn_delta"]=sn2
        # Rank direct SN closure while protecting acoustic and BAO continuity.
        r["stage1_score"]=sn2 + 1600*abs(r["ellA_frac"]) + 0.10*max(0,r["bao_chi2"]-EDGE024_BAO_CHI)
    rows.append(r)
    if i%16==0: print("SNCLOSE_PROGRESS",i,json.dumps(r,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"snclosure_stage1_all.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
# Acoustic/BAO guard; relax only if too few.
guard=ok[(ok.ellA_frac.abs()<0.0020) & (ok.bao_chi2<20.0)].copy()
if len(guard)<32: guard=ok.copy()
short=guard.nsmallest(32,"stage1_score").copy()
# Always include exact center for consistent Planck-lite normalization.
centerrow=pd.DataFrame([{**center,"ellA_frac":0.0,
                         "sn_delta_pantheonplus":center["sn_pantheonplus"]-LOCAL_SN["pantheonplus"],
                         "sn_delta_union3":center["sn_union3"]-LOCAL_SN["union3"],
                         "sn_delta_desy5":center["sn_desy5"]-LOCAL_SN["desy5"],
                         "bao_delta_local":center["bao_chi2"]-LOCAL_BAO_CHI,
                         "second_best_sn_delta":sorted([
                           center["sn_pantheonplus"]-LOCAL_SN["pantheonplus"],
                           center["sn_union3"]-LOCAL_SN["union3"],
                           center["sn_desy5"]-LOCAL_SN["desy5"]])[1],
                         "stage1_score":999.}])
short=pd.concat([centerrow,short],ignore_index=True)
short.to_csv(OUT/"snclosure_stage1_shortlist.csv",index=False)

# Planck native-lite stage
high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

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
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-9})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def planck_for(row):
    label="cl_"+str(row["id"]); root=str(TMP/(label+"_")); ip=TMP/(label+".ini")
    ip.write_text(ini_text(root,float(row.A),float(row.B),float(row.H0),True))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec=row.to_dict(); rec["planck_status"]="FAIL"
    cl=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and cl.exists():
        rec.update(pscore(cl)); rec["planck_status"]="OK"
    else:
        rec["planck_error"]=cp.stdout[-800:].replace("\n"," | ")
    for q in TMP.glob(label+"_*"):
        try:q.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    print("SNCLOSE_PLANCK_POINT",json.dumps(rec,sort_keys=True),flush=True)
    return rec

p_rows=[planck_for(r) for _,r in short.iterrows()]
pdf=pd.DataFrame(p_rows)
centerP=float(pdf.loc[pdf.id=="center","chi2_planck"].iloc[0])
pdf["delta_planck_vs_center"]=pdf.chi2_planck-centerP
# Conservative DESI movement proxy: one quarter of the BAO chi2 movement.
pdf["pd_proxy"]=EDGE024_PD_FAIR+pdf.delta_planck_vs_center+0.25*(pdf.bao_chi2-EDGE024_BAO_CHI)
pdf["second_best_sn_delta"]=pdf[["sn_delta_pantheonplus","sn_delta_union3","sn_delta_desy5"]].apply(
    lambda r: sorted(map(float,r))[1],axis=1)
pdf["goal_score"]=np.maximum(pdf.pd_proxy,pdf.second_best_sn_delta)
pdf=pdf.sort_values(["goal_score","second_best_sn_delta","pd_proxy"])
pdf.to_csv(OUT/"snclosure_planck_shortlist.csv",index=False)
best=pdf.head(8).to_dict("records")
summary={
 "center_pd_fair":EDGE024_PD_FAIR,
 "center_planck_lite":centerP,
 "center_ellA_proxy":ELL0,
 "n_stage1_ok":int(len(ok)),
 "n_planck":int((pdf.planck_status=="OK").sum()),
 "best":best,
 "goal":"pd_proxy<0 and second_best_sn_delta<0 (at least two exact SN datasets beat LCDM); then exact full-Plik + raw DESI promotion required"
}
(OUT/"snclosure_summary.json").write_text(json.dumps(summary,indent=2))
print("SNCLOSE_BEST",json.dumps(best,sort_keys=True),flush=True)
print("SNCLOSE_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
