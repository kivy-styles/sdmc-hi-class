#!/usr/bin/env python3
"""
Full-SDMC multi-observable candidate generator around the exact S117 / AF11 basin.

This is a screening workflow, not the final exact promotion. It preserves the
GitHub Part IX comparison convention and ranks candidates before promotion to:
  * full Planck Plik TT/TE/EE + low-l + lensing nuisance profile,
  * official DESI DR1 P0/P2/P4 full-shape,
  * exact Pantheon+ / Union3 / DES-Y5 covariance,
  * calibrated Pantheon+SH0ES ladder,
  * uniform BAO / BBN / weak-lensing / age battery,
  * conditional JWST/FRESCO nonlinear-collapse likelihood-equivalent test.

The JWST block uses the newly derived local mapper branch as a conditional
nonlinear completion. The minimal identification g3 = A_F is kept explicit;
it is not silently promoted to a fundamental identity.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/fulljoint_corridor_refine")
OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
BAOROOT=Path("bao_data")

TCMB=2.7255
CAL_SIGMA=0.0025
OR=4.17998772e-5
C_KMS=299792.458

# Exact benchmark / control ledger carried forward from Part IX.
S117_LITE=1020.917005539309
LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}
LOCAL_BAO=13.400290718544086
LOCAL_BBN=5.78515
LOCAL_S8=0.80746
LOCAL_AGE=13.75343
LADDER_H0=73.49
LADDER_SIG=1.048

# Fixed full-SDMC transition constants retained from the successful basin.
D0=0.34231919445927034
POWER=0.45
LAM=18.40625
ZT=16.10
DNT=0.5
TAUA=0.25
TAUB=1.5

# Targeted anchors: exact S117 plus the two best lower-S8 points from the first broad scan.
S117=dict(id="S117",H0=69.71482083084993,ob=0.022083219194622913,oc=0.12299536722293603,
      As=2.117602412937069e-9,ns=0.9625227132590487,tau=0.055202901571989066,
      AF=0.02114389337040484,zc=3.715736017236486,width=0.35696151830255984,
      Dfloor=0.052829781055450435,Alate=0.024822032167576252,Blate=0.01264704344701022)
FJ7=dict(id="fj007",H0=69.9257593356073,ob=0.022026450499221684,oc=0.12225877693593502,
      As=2.0788681025616825e-9,ns=0.9653349212072789,tau=0.04896082428097725,
      AF=0.02155169515684247,zc=3.7798996351845564,width=0.33331606168299915,
      Dfloor=0.04623167278245091,Alate=0.022394198394566776,Blate=0.012352856796979905)
FJ28=dict(id="fj028",H0=70.10553486775606,ob=0.022147933021215722,oc=0.12157599684279412,
      As=2.0819216857664285e-9,ns=0.962911429583095,tau=0.052764227174222474,
      AF=0.02101232836768031,zc=3.8788957125507295,width=0.3907623110339046,
      Dfloor=0.04894181532226503,Alate=0.018275821322947743,Blate=0.011977504090964795)
KEYS=["H0","ob","oc","As","ns","tau","AF","zc","width","Dfloor","Alate","Blate"]
def blend(a,b,t,label):
    d={k:(a[k]+t*(b[k]-a[k])) for k in KEYS}; d["id"]=label; return d
SEEDS=[S117,FJ7,FJ28]
for t in (0.20,0.35,0.50,0.65,0.80):
    SEEDS.append(blend(S117,FJ7,t,f"s7_t{int(round(t*100)):02d}"))
for t in (0.20,0.35,0.50,0.65):
    SEEDS.append(blend(S117,FJ28,t,f"s28_t{int(round(t*100)):02d}"))

# Focus the Sobol cloud on the low-S8 corridor rather than resampling the entire Part IX basin.
LOW=np.array([69.66,0.02200,0.12145,2.055e-9,0.9610,0.0460,
              0.0198,3.55,0.315,0.0455,0.0180,0.01175])
HIGH=np.array([70.12,0.02216,0.12310,2.119e-9,0.9680,0.0565,
               0.0222,3.93,0.395,0.0555,0.0252,0.01335])

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def class_table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"covariance mismatch {path}: {len(a)} != {n*n}")
    return a.reshape(n,n)

def setup_cov(C):
    cf=cho_factor(C,lower=True,check_finite=False)
    one=np.ones(C.shape[0]); u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one@u)

def prof(pack,r):
    cf,one,den=pack
    y=cho_solve(cf,r,check_finite=False)
    return float(r@y-(one@y)**2/den)

# SN data
pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm))
m=ppz>0.01; ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]
pppack=setup_cov(Cpp)
p=SNROOT/"Union3/lcparam_full.txt"
cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
cm={c.lower():c for c in un.columns}; unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float)
unpack=setup_cov(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))
de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
depack=setup_cov(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

# BAO
bao_rows=[]
for ln in (BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
    if ln.strip() and not ln.startswith("#"):
        z,v,q=ln.split(); bao_rows.append((float(z),float(v),q))
bao_cov=np.loadtxt(BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
bao_inv=np.linalg.inv(bao_cov); bao_obs=np.array([x[1] for x in bao_rows])

def mu_bg(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."])
    DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def sn_scores(bg):
    return {
      "pantheonplus":prof(pppack,ppm-mu_bg(bg,ppz,ppzh)),
      "union3":prof(unpack,unm-mu_bg(bg,unz,unz)),
      "desy5":prof(depack,dem-mu_bg(bg,dez,dezh))
    }

def rd_from(bg,th):
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); crosses=[]
    for i in range(len(zz)-1):
        if (kb[i]-1)*(kb[i+1]-1)<=0 and kb[i]!=kb[i+1]:
            q=(1-kb[i])/(kb[i+1]-kb[i]); crosses.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not crosses: raise RuntimeError("no drag crossing")
    zd=float(crosses[0]); rd=float(np.interp(zd,bg.z,bg["comov.snd.hrz."]))
    return rd,zd

def bao_score(bg,rd):
    pred=[]
    for z,_,q in bao_rows:
        DM=float(np.interp(z,bg.z,bg["comov. dist."]))
        H=float(np.interp(z,bg.z,bg["H [1/Mpc]"]))
        DH=1/H; DV=(z*DM*DM*DH)**(1/3)
        pred.append({"DM_over_rs":DM/rd,"DH_over_rs":DH/rd,"DV_over_rs":DV/rd}[q])
    d=np.array(pred)-bao_obs
    return float(d@bao_inv@d)

# Same transparent BBN screen as the Part IX uniform battery.
ob_ref=0.02239952; Y_ref=0.246990; D_ref=2.43955; Li_ref=5.54707; dneff_ref=0.062049
dY_ref=0.247820-0.246990; dD_ref=2.45993-2.43955; dLi_ref=5.51555-5.54707
Y_obs,Y_sig=0.2450,0.0034; D_obs,D_sig=2.527,0.038; Li_obs=1.60

def bbn(ob):
    eta=273.9*ob; eta0=273.9*ob_ref
    Y=Y_ref+0.0016*(eta-eta0); D=D_ref*(ob_ref/ob)**1.6; Li=Li_ref*(ob/ob_ref)**2
    fr=4.0/LAM**2; dneff=(43/7)*fr/(1-fr); scale=dneff/dneff_ref
    Y+=dY_ref*scale; D+=dD_ref*scale*(D/D_ref); Li+=dLi_ref*scale*(Li/Li_ref)
    chi=((Y-Y_obs)/Y_sig)**2+((D-D_obs)/D_sig)**2
    return dict(chi2=float(chi),Yp=float(Y),DH_1e5=float(D),Li_over_plateau=float(Li/Li_obs),delta_Neff=float(dneff))

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); ppv=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); ppv[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":ppv}

def planck_score(path):
    harr,dls=load_cls(path)
    def f(A):
        v=float(high.chi_squared(harr,A_planck=A))
        v+=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        v+=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        v+=float(-2*lens.log_likelihood(dls,**lp))
        v+=((A-1)/CAL_SIGMA)**2
        return v
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def sigma8_from_pk(path):
    a=np.loadtxt(path); k=a[:,0]; P=a[:,1]; x=8*k
    W=np.where(np.abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
    return float(np.sqrt(np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))))

def wl_chi(S8,mu,sig):
    return float(((S8-mu)/sig)**2)

# JWST is kept independent of the linear hi_class coefficient A_F.
# The nonlinear mapper coefficient g3 is a separate candidate UV coefficient.
JWST_G3_CENTRAL=0.01573642223
JWST_G3_EPS20_95=0.03496982717777778
JWST_G3_EPS20_68=0.04429511442518519
def jwst_independent():
    return dict(g3_central=JWST_G3_CENTRAL,
                g3_epsilon20_95=JWST_G3_EPS20_95,
                g3_epsilon20_68=JWST_G3_EPS20_68,
                status="independent nonlinear mapper target; does not alter linear screen")

def ini(root,c):
    h=c["H0"]/100.; ox=1.-(c["ob"]+c["oc"]+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={c['H0']:.17g}
omega_b={c['ob']:.17g}
omega_cdm={c['oc']:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={c['As']:.17e}
n_s={c['ns']:.17g}
tau_reio={c['tau']:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={c['AF']:.17g},{c['zc']:.17g},{c['width']:.17g},{D0:.17g},{POWER:.17g},{c['Dfloor']:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM:.17g},{ZT:.17g},{DNT:.17g},{c['Alate']:.17g},{TAUA:.17g},{c['Blate']:.17g},{TAUB:.17g}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
modes=s
output=tCl,pCl,lCl,mPk
lensing=yes
l_max_scalars=3000
P_k_max_h/Mpc=3.0
z_pk=0
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

def run(c):
    label=c["id"]; root=str(OUT/(label+"_")); ip=OUT/(label+".ini"); ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=420)
    rec=dict(c); rec.update(status="FAIL",returncode=cp.returncode)
    try:
        bg=class_table(root+"00_background.dat"); th=class_table(root+"00_thermodynamics.dat")
        D=bg["kin (D)"].to_numpy(); cs=bg["c_s^2"].to_numpy()
        rec["D_min"]=float(D.min()); rec["cs2_min"]=float(cs.min()); rec["cs2_max"]=float(cs.max())
        rec["noslip_max"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        stable=cp.returncode==0 and rec["D_min"]>0 and rec["cs2_min"]>0 and rec["cs2_max"]<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            pc,Ap=planck_score(root+"00_cl_lensed.dat")
            pks=sorted(OUT.glob(label+"_*pk.dat"))
            sig8=sigma8_from_pk(pks[0])
            i0=int(np.argmin(np.abs(bg.z.to_numpy())))
            Om0=float(bg.iloc[i0]["Omega_m(z)"]); S8=sig8*math.sqrt(Om0/0.3)
            sn=sn_scores(bg); rd,zd=rd_from(bg,th); bc=bao_score(bg,rd); bb=bbn(c["ob"])
            age0=float(np.interp(0.0,bg.z,bg["proper time [Gyr]"]))
            j=jwst_independent()
            rec.update(status="OK",chi2_planck_lite=pc,A_planck=Ap,delta_planck_vs_S117=pc-S117_LITE,
                       sn_pp_delta=sn["pantheonplus"]-LOCAL_SN["pantheonplus"],
                       sn_u3_delta=sn["union3"]-LOCAL_SN["union3"],
                       sn_d5_delta=sn["desy5"]-LOCAL_SN["desy5"],
                       sn_worst_delta=max(sn["pantheonplus"]-LOCAL_SN["pantheonplus"],sn["union3"]-LOCAL_SN["union3"],sn["desy5"]-LOCAL_SN["desy5"]),
                       rd_Mpc=rd,z_drag=zd,bao_chi2=bc,bao_delta=bc-LOCAL_BAO,
                       bbn_chi2=bb["chi2"],bbn_delta=bb["chi2"]-LOCAL_BBN,bbn_Li_factor=bb["Li_over_plateau"],
                       sigma8=sig8,Omega_m0=Om0,S8=S8,
                       wl_DESY3_delta=wl_chi(S8,0.776,0.017)-wl_chi(LOCAL_S8,0.776,0.017),
                       wl_HSC_delta=wl_chi(S8,0.805,0.018)-wl_chi(LOCAL_S8,0.805,0.018),
                       wl_KiDS_delta=wl_chi(S8,0.815,0.020)-wl_chi(LOCAL_S8,0.815,0.020),
                       age_t0_Gyr=age0,age_delta=age0-LOCAL_AGE,
                       shoes_proxy_H0_pull=(c["H0"]-LADDER_H0)/LADDER_SIG,
                       shoes_proxy_chi2=((c["H0"]-LADDER_H0)/LADDER_SIG)**2,
                       **{f"jwst_{k}":v for k,v in j.items()})
            # Screening objective: exact promotions, not this scalar, decide the scientific winner.
            planck_pen=max(0.,rec["delta_planck_vs_S117"]-0.10)
            sn_pen=max(0.,rec["sn_worst_delta"])
            bao_pen=max(0.,rec["bao_delta"])
            wl_pen=max(0.,(S8-0.815)/0.015)**2
            age_pen=max(0.,(13.45-age0)/0.1)**2
            jw_pen=0.0
            shoes_pen=rec["shoes_proxy_chi2"]/15.
            rec["screen_objective"]=float(planck_pen+0.45*sn_pen+0.15*bao_pen+0.30*wl_pen+0.10*age_pen+0.10*jw_pen+0.08*shoes_pen)
    except Exception as e:
        rec["error"]=repr(e)
    if rec["status"]!="OK":
        rec["error"]=rec.get("error",cp.stdout[-1200:].replace("\n"," | "))
    print("FULLJOINT_SCREEN_POINT",json.dumps(rec,sort_keys=True),flush=True)
    # Keep only summaries to make artifact small.
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return rec

rows=[]
for c in SEEDS: rows.append(run(c))

sob=qmc.Sobol(d=len(KEYS),scramble=True,seed=20261007)
cloud=qmc.scale(sob.random_base2(m=7),LOW,HIGH) # 128 targeted corridor candidates
for i,x in enumerate(cloud):
    c={k:float(v) for k,v in zip(KEYS,x)}
    c["id"]=f"fj{i:03d}"
    rows.append(run(c))

df=pd.DataFrame(rows)
df.to_csv(OUT/"fulljoint_screen.csv",index=False)
ok=df[df.status=="OK"].copy()
if len(ok):
    # Hard pre-promotion gates keep points close enough to the established exact basin.
    gate=ok[(ok.chi2_planck_lite<=S117_LITE+2.0)&(ok.sn_worst_delta<=1.2)&(ok.bao_delta<=0.5)&(ok.age_t0_Gyr>=13.45)&(ok.S8<=0.835)]
    if len(gate)<8: gate=ok.sort_values(["S8","screen_objective"]).head(20)
    promote=gate.sort_values(["screen_objective","S8","sn_worst_delta"]).head(10)
else:
    promote=ok
promote.to_csv(OUT/"promotion_shortlist.csv",index=False)
promote_cols=["id"]+KEYS+["chi2_planck_lite","delta_planck_vs_S117","sn_pp_delta","sn_u3_delta","sn_d5_delta",
                         "bao_delta","bbn_delta","S8","wl_DESY3_delta","age_t0_Gyr","shoes_proxy_chi2",
                         "jwst_A_M","jwst_epsilon95_equiv","screen_objective"]
summary={"n_total":len(rows),"n_ok":int(len(ok)),"n_shortlist":int(len(promote)),
         "promotion_rule":"Top stable points near S117/AF11 are not declared winners until exact full Plik + official raw DESI + exact SN + calibrated SH0ES are run.",
         "jwst_status":"independent nonlinear mapper target; g3 is not identified with hi_class A_F",
         "shortlist":promote[[c for c in promote_cols if c in promote.columns]].to_dict("records")}
(OUT/"fulljoint_screen_summary.json").write_text(json.dumps(summary,indent=2))
print("FULLJOINT_SCREEN_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
