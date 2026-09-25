#!/usr/bin/env python3
"""
Uniform downstream observational battery for all successful Part IX exact candidates.

Gates:
  * Pantheon+, Union3, DES-Y5 shape-only SN covariance
  * DESI DR1 2024 Gaussian BAO (12-vector covariance)
  * PRIMAT-calibrated local BBN response screen (Yp, D/H; Li reported separately)
  * galaxy weak-lensing compressed S8 comparisons
  * exact cosmic age t0 and t(z)
  * old passive-galaxy age-margin stress tests

Every result is compared against the same optimized LCDM local021 control.
No statistically overlapping gates are summed together.
"""
from pathlib import Path
import json, math, re
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.integrate import quad
import camb

C_KMS = 299792.458
MPC_KM = 3.0856775814913673e19
SEC_GYR = 365.25*86400.0*1e9

ROOT = Path("inputs")
SNROOT = Path("sn_data")
BAOROOT = Path("bao_data")
OUT = Path("output/uniform_battery")
OUT.mkdir(parents=True, exist_ok=True)

PLANCK_DELTA_LOCAL021 = 3.0095900971796254
CANDS = {
 "bo002_exact": dict(dir="bo002", omega_b=0.022083219194622913, lambda_e=18.40625,
                     pd_fair=0.0),
}
LOCAL021 = dict(H0=68.56858744695782,
                omega_b=0.022406369378007947,
                omega_cdm=0.11824151052483357,
                A_s=2.05109266435559e-9,
                n_s=0.9649593164240942,
                tau=0.044985382026527077,
                N_ur=3.046, Tcmb=2.7255, YHe=0.2453)

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith('#') and re.search(r'1\s*:',l)][-1].lstrip('#').strip()
    ms=list(re.finditer(r'(\d+)\s*:\s*',hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values('z')

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"covariance mismatch {path}: {len(a)} vs {n*n}")
    return a.reshape(n,n)

def setup_cov(C):
    cf=cho_factor(C,lower=True,check_finite=False)
    one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def chi_profile_intercept(pack,r):
    cf,one,den=pack
    y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

def rd_from(bg,th):
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); crosses=[]
    for i in range(len(zz)-1):
        if (kb[i]-1)*(kb[i+1]-1)<=0 and kb[i]!=kb[i+1]:
            q=(1-kb[i])/(kb[i+1]-kb[i]); crosses.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not crosses: raise RuntimeError("no kappa_b=1 drag crossing")
    zd=float(crosses[0])
    rd=float(np.interp(zd,bg.z,bg['comov.snd.hrz.']))
    return rd,zd

# ---------- SN data ----------
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

def mu_bg(bg,z,zh):
    DM=np.interp(z,bg.z,bg['comov. dist.'])
    DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def sn_scores_bg(bg):
    return {
      "pantheonplus":chi_profile_intercept(pppack,ppm-mu_bg(bg,ppz,ppzh)),
      "union3":chi_profile_intercept(unpack,unm-mu_bg(bg,unz,unz)),
      "desy5":chi_profile_intercept(depack,dem-mu_bg(bg,dez,dezh)),
    }

# ---------- local021 CAMB ----------
pars=camb.CAMBparams()
pars.set_cosmology(H0=LOCAL021["H0"],ombh2=LOCAL021["omega_b"],omch2=LOCAL021["omega_cdm"],
                   mnu=0.0,num_massive_neutrinos=0,nnu=LOCAL021["N_ur"],
                   TCMB=LOCAL021["Tcmb"],YHe=LOCAL021["YHe"],tau=LOCAL021["tau"])
pars.InitPower.set_params(As=LOCAL021["A_s"],ns=LOCAL021["n_s"])
pars.set_matter_power(redshifts=[0.0],kmax=3.0)
pars.NonLinear = camb.model.NonLinear_none
lres=camb.get_results(pars)
ldp=lres.get_derived_params()
local_rd=float(ldp['rdrag'])

def mu_local(z,zh):
    DA=np.array([lres.angular_diameter_distance(float(x)) for x in np.atleast_1d(z)])
    return 5*np.log10((1+np.atleast_1d(zh))*(1+np.atleast_1d(z))*DA)

local_sn={
 "pantheonplus":chi_profile_intercept(pppack,ppm-mu_local(ppz,ppzh)),
 "union3":chi_profile_intercept(unpack,unm-mu_local(unz,unz)),
 "desy5":chi_profile_intercept(depack,dem-mu_local(dez,dezh)),
}

# ---------- DESI DR1 Gaussian BAO ----------
mean_path=BAOROOT/'desi_2024_gaussian_bao_ALL_GCcomb_mean.txt'
bao_rows=[]
for ln in mean_path.read_text().splitlines():
    if not ln.strip() or ln.startswith('#'): continue
    z,v,q=ln.split()
    bao_rows.append((float(z),float(v),q))
bao_cov=np.loadtxt(BAOROOT/'desi_2024_gaussian_bao_ALL_GCcomb_cov.txt')
bao_inv=np.linalg.inv(bao_cov)
bao_obs=np.array([x[1] for x in bao_rows])

def bao_pred_bg(bg,rd):
    vals=[]
    for z,_,q in bao_rows:
        DM=float(np.interp(z,bg.z,bg['comov. dist.']))
        H=float(np.interp(z,bg.z,bg['H [1/Mpc]']))
        DH=1.0/H
        DV=(z*DM*DM*DH)**(1/3)
        vals.append({"DM_over_rs":DM/rd,"DH_over_rs":DH/rd,"DV_over_rs":DV/rd}[q])
    return np.array(vals)

def bao_pred_local(rd):
    vals=[]
    for z,_,q in bao_rows:
        DM=float(lres.comoving_radial_distance(z))
        DH=C_KMS/float(lres.hubble_parameter(z))
        DV=(z*DM*DM*DH)**(1/3)
        vals.append({"DM_over_rs":DM/rd,"DH_over_rs":DH/rd,"DV_over_rs":DV/rd}[q])
    return np.array(vals)

def bao_chi(pred):
    d=pred-bao_obs
    return float(d@bao_inv@d)
local_bao_pred=bao_pred_local(local_rd)
local_bao_chi=bao_chi(local_bao_pred)

# ---------- BBN calibrated local-sensitivity screen ----------
# Calibrated to the Part VIII PRIMAT-local benchmark at omega_b_ref and lambda_e=20.
ob_ref=0.02239952
Y_ref=0.246990
D_ref=2.43955 # in 1e-5 units
Li_ref=5.54707 # in 1e-10 units
dneff_ref=0.062049
dY_ref=0.247820-0.246990
dD_ref=2.45993-2.43955
dLi_ref=5.51555-5.54707
# transparent abundance data adopted for this controlled screen
Y_obs,Y_sig=0.2450,0.0034
D_obs,D_sig=2.527,0.038
Li_obs=1.60

def bbn_standard(ob):
    eta=273.9*ob; eta0=273.9*ob_ref
    Y=Y_ref+0.0016*(eta-eta0)
    D=D_ref*(ob_ref/ob)**1.6
    Li=Li_ref*(ob/ob_ref)**2.0
    return Y,D,Li

def bbn_tracker(ob,lambda_e=None):
    Y,D,Li=bbn_standard(ob)
    if lambda_e is None:
        dneff=0.0
    else:
        fr=4.0/lambda_e**2
        dneff=(43/7)*fr/(1-fr)
        scale=dneff/dneff_ref
        Y += dY_ref*scale
        D += dD_ref*scale*(D/D_ref)
        Li += dLi_ref*scale*(Li/Li_ref)
    chi=((Y-Y_obs)/Y_sig)**2+((D-D_obs)/D_sig)**2
    return dict(Yp=Y,DH_1e5=D,LiH_1e10=Li,Li_over_plateau=Li/Li_obs,
                delta_Neff_bg=dneff,chi2_Yp_DH=float(chi))

local_bbn=bbn_tracker(LOCAL021["omega_b"],None)

# ---------- sigma8 / S8 weak-lensing compression ----------
def sigma8_from_pk(path):
    a=np.loadtxt(path); k=a[:,0]; P=a[:,1]; x=8.0*k
    W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
    sig2=np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))
    return float(np.sqrt(sig2))

# CAMB local021 matter power, same top-hat convention
kh,zz,pk=lres.get_matter_power_spectrum(minkh=1e-5,maxkh=3.0,npoints=600)
x=8*kh; W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
local_sigma8=float(np.sqrt(np.trapezoid(kh**3*pk[0]/(2*np.pi**2)*W**2,np.log(kh))))
h=LOCAL021["H0"]/100
local_Om=(LOCAL021["omega_b"]+LOCAL021["omega_cdm"])/h**2
local_S8=local_sigma8*math.sqrt(local_Om/0.3)

WL = {
 "DES_Y3_3x2": dict(mu=0.776,lo=0.017,hi=0.017),
 "KiDS_Legacy": dict(mu=0.815,lo=0.021,hi=0.016),
 "HSC_Y3_DESIcal": dict(mu=0.805,lo=0.018,hi=0.018),
}
def wl_scores(S8):
    out={}
    for k,d in WL.items():
        sig=d["hi"] if S8>=d["mu"] else d["lo"]
        pull=(S8-d["mu"])/sig
        out[k]=dict(S8_obs=d["mu"],sigma_used=sig,pull=float(pull),chi2=float(pull*pull))
    return out
local_wl=wl_scores(local_S8)

# ---------- local021 age ----------
Tc=LOCAL021["Tcmb"]; h=LOCAL021["H0"]/100
Og=2.4728e-5/h**2*(Tc/2.7255)**4
Or=Og*(1+0.22710731766*LOCAL021["N_ur"])
Om=local_Om
Ol=1-Om-Or
H0_s=LOCAL021["H0"]/MPC_KM
def local_age(z):
    # Integrate in ln(a), not directly to an enormous z upper bound.
    # dt = d ln(a) / H(a); this is numerically stable through radiation era.
    xmax=math.log(1.0/(1.0+z))
    xmin=-32.0
    def fx(x):
        a=math.exp(x)
        E=math.sqrt(Or*a**-4+Om*a**-3+Ol)
        return 1.0/(H0_s*E)
    sec=quad(fx,xmin,xmax,epsabs=0,epsrel=2e-9,limit=500)[0]
    return sec/SEC_GYR

AGE_Z=[0.0,0.5,1.0,1.432,1.552,2.0,3.0,6.0,10.0,20.0]
local_ages={str(z):local_age(z) for z in AGE_Z}

old_galaxies = {
 "53W069_stress": dict(z=1.432,age_Gyr=4.0,note="upper end of historical 3-4 Gyr stellar-population estimate"),
 "53W069_mid": dict(z=1.432,age_Gyr=3.5,note="midpoint of historical 3-4 Gyr estimate"),
 "53W091_revised": dict(z=1.552,age_Gyr=1.8,note="upper end of later 1.4-1.8 Gyr synthesis estimate"),
 "53W091_historical": dict(z=1.552,age_Gyr=3.5,note="historical minimum-age stress claim"),
}
def oldgal_scores(agefunc):
    out={}
    for name,d in old_galaxies.items():
        tu=float(agefunc(d["z"]))
        out[name]=dict(z=d["z"],stellar_age_Gyr=d["age_Gyr"],universe_age_Gyr=tu,
                       margin_Gyr=tu-d["age_Gyr"],ratio=d["age_Gyr"]/tu,note=d["note"])
    return out
local_old=oldgal_scores(local_age)

# ---------- candidate loop ----------
rows=[]; details={}
for label,meta in CANDS.items():
    d=ROOT/meta["dir"]
    bg=table(d/'desi_cov_00_background.dat')
    th=table(d/'desi_cov_00_thermodynamics.dat')
    rd,zd=rd_from(bg,th)
    sn=sn_scores_bg(bg)
    pred=bao_pred_bg(bg,rd)
    bchi=bao_chi(pred)
    bbn=bbn_tracker(meta["omega_b"],meta["lambda_e"])
    sig8=sigma8_from_pk(d/'desi_cov_00_z1_pk.dat')
    i0=int(np.argmin(np.abs(bg.z.to_numpy())))
    Om0=float(bg.iloc[i0]['Omega_m(z)'])
    H0=float(bg.iloc[i0]['H [1/Mpc]']*C_KMS)
    S8=sig8*math.sqrt(Om0/0.3)
    wl=wl_scores(S8)
    def age_bg(z):
        zz=bg.z.to_numpy(); tt=bg['proper time [Gyr]'].to_numpy()
        order=np.argsort(zz)
        return float(np.interp(z,zz[order],tt[order]))
    ages={str(z):age_bg(z) for z in AGE_Z}
    old=oldgal_scores(age_bg)
    details[label]=dict(H0=H0,omega_b=meta["omega_b"],lambda_e=meta["lambda_e"],
                        rd_Mpc=rd,z_drag=zd,sn=sn,bao_chi2=bchi,bao_vector=pred.tolist(),
                        bbn=bbn,sigma8=sig8,Omega_m0=Om0,S8=S8,weak_lensing=wl,
                        ages_Gyr=ages,old_galaxies=old,pd_fair=meta["pd_fair"])
    dbao=bchi-local_bao_chi
    dpp=sn["pantheonplus"]-local_sn["pantheonplus"]
    du3=sn["union3"]-local_sn["union3"]
    dd5=sn["desy5"]-local_sn["desy5"]
    pbao=PLANCK_DELTA_LOCAL021+dbao
    row=dict(candidate=label,H0=H0,rd_Mpc=rd,
             exact_planck_delta_vs_local021=PLANCK_DELTA_LOCAL021,
             exact_planck_plus_desi_bao_delta=pbao,
             exact_planck_bao_pp_delta=pbao+dpp,
             exact_planck_bao_u3_delta=pbao+du3,
             exact_planck_bao_desy5_delta=pbao+dd5,
             pd_fair=meta["pd_fair"],
             sn_pp_delta_vs_local021=dpp,
             sn_u3_delta_vs_local021=du3,
             sn_desy5_delta_vs_local021=dd5,
             joint_pd_pp_fair=meta["pd_fair"]+sn["pantheonplus"]-local_sn["pantheonplus"],
             joint_pd_u3_fair=meta["pd_fair"]+sn["union3"]-local_sn["union3"],
             joint_pd_desy5_fair=meta["pd_fair"]+sn["desy5"]-local_sn["desy5"],
             bao_chi2=bchi,bao_delta_vs_local021=dbao,
             bbn_chi2=bbn["chi2_Yp_DH"],bbn_delta_vs_local021=bbn["chi2_Yp_DH"]-local_bbn["chi2_Yp_DH"],
             bbn_deltaNeff=bbn["delta_Neff_bg"],bbn_Li_factor=bbn["Li_over_plateau"],
             sigma8=sig8,Omega_m0=Om0,S8=S8,
             wl_DESY3_chi2=wl["DES_Y3_3x2"]["chi2"],
             wl_DESY3_delta_vs_local021=wl["DES_Y3_3x2"]["chi2"]-local_wl["DES_Y3_3x2"]["chi2"],
             wl_KiDS_chi2=wl["KiDS_Legacy"]["chi2"],
             wl_KiDS_delta_vs_local021=wl["KiDS_Legacy"]["chi2"]-local_wl["KiDS_Legacy"]["chi2"],
             wl_HSC_chi2=wl["HSC_Y3_DESIcal"]["chi2"],
             wl_HSC_delta_vs_local021=wl["HSC_Y3_DESIcal"]["chi2"]-local_wl["HSC_Y3_DESIcal"]["chi2"],
             age_t0_Gyr=ages["0.0"],age_t0_delta_vs_local021=ages["0.0"]-local_ages["0.0"],
             age_z1432_Gyr=ages["1.432"],age_z1552_Gyr=ages["1.552"],
             W069_4Gyr_margin=old["53W069_stress"]["margin_Gyr"],
             W091_35Gyr_margin=old["53W091_historical"]["margin_Gyr"])
    rows.append(row)

summary=pd.DataFrame(rows).sort_values("candidate")
summary.to_csv(OUT/"successful_candidates_uniform_battery.csv",index=False)
(OUT/"successful_candidates_uniform_battery.json").write_text(json.dumps(details,indent=2))

local=dict(H0=LOCAL021["H0"],rd_Mpc=local_rd,sn=local_sn,bao_chi2=local_bao_chi,
           bao_vector=local_bao_pred.tolist(),bbn=local_bbn,sigma8=local_sigma8,
           Omega_m0=local_Om,S8=local_S8,weak_lensing=local_wl,ages_Gyr=local_ages,
           old_galaxies=local_old,pd_total=-31.152232235955807)
(OUT/"local021_uniform_control.json").write_text(json.dumps(local,indent=2))

print("BO002_EXACT_BATTERY_LOCAL021",json.dumps(local,sort_keys=True),flush=True)
for _,r in summary.iterrows():
    print("UNIFORM_BATTERY_CANDIDATE",json.dumps(r.to_dict(),sort_keys=True),flush=True)

# A compact verdict ledger: no cross-gate summation because FS/BAO overlap and WL is compressed.
verdict=[]
for r in rows:
    verdict.append(dict(
      candidate=r["candidate"],
      Planck_DESI_crosses_local021=(r["exact_planck_plus_desi_bao_delta"]<0),
      Pantheon_joint_crosses_local021=(r["exact_planck_bao_pp_delta"]<0),
      Union3_joint_crosses_local021=(r["exact_planck_bao_u3_delta"]<0),
      DESY5_joint_crosses_local021=(r["exact_planck_bao_desy5_delta"]<0),
      BAO_better_than_local021=(r["bao_delta_vs_local021"]<0),
      BBN_screen_better_than_local021=(r["bbn_delta_vs_local021"]<0),
      KiDS_compressed_better_than_local021=(r["wl_KiDS_delta_vs_local021"]<0),
      age_older_than_local021=(r["age_t0_delta_vs_local021"]>0),
      passes_53W069_4Gyr_age=(r["W069_4Gyr_margin"]>0),
      passes_53W091_35Gyr_age=(r["W091_35Gyr_margin"]>0),
    ))
pd.DataFrame(verdict).to_csv(OUT/"uniform_gate_verdict.csv",index=False)
print("UNIFORM_BATTERY_DONE",flush=True)
