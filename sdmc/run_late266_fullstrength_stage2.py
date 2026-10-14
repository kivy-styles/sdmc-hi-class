#!/usr/bin/env python3
"""
Full-SDMC second-stage candidate generator around the late266 successor basin.

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

OUT=Path("output/late266_fullstrength_stage2")
OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
BAOROOT=Path("bao_data")

TCMB=2.7255
CAL_SIGMA=0.0025
OR=4.17998772e-5
C_KMS=299792.458

# Exact benchmark / control ledger carried forward from Part IX.
S036_LITE=1021.1344170797494
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
POWER=1.0
LAM=18.40625
DNT=0.5
TAUA=0.25
TAUB=1.5

# Anchors from Manuscript-A late266 and the first-stage successor screen.
SEEDS=[
 dict(id="late266",H0=69.45160505326464,ob=0.02208511574370519,oc=0.12298509428428875,
      As=2.1197359302086e-9,ns=0.962555355281866,tau=0.055202901571989066,
      AF=0.023604633340554453,zc=3.4328776987879466,width=0.36879721635160295,
      Dfloor=0.050275813910968366,zt=16.173189924377947,
      Alate=0.013036667098291217,Blate=0.010544837576337158),
 dict(id="s036",H0=69.45160505326464,ob=0.02208511574370519,oc=0.12298509428428875,
      As=2.1153164253417825e-9,ns=0.963824442274234,tau=0.05415934741962701,
      AF=0.020326851338950964,zc=3.6185491493767987,width=0.37530905244079676,
      Dfloor=0.05087322041005618,zt=16.173189924377947,
      Alate=0.013036667098291217,Blate=0.010544837576337158),
 dict(id="s024",H0=69.45160505326464,ob=0.02208511574370519,oc=0.12298509428428875,
      As=2.10569837679848e-9,ns=0.9622208507576201,tau=0.051880733001045884,
      AF=0.019676369993974624,zc=3.6132100268690777,width=0.3528031259137139,
      Dfloor=0.05608912437767512,zt=16.173189924377947,
      Alate=0.013036667098291217,Blate=0.010544837576337158)
]

# Tight full-SDMC box.  Ordinary parameters are allowed to move only slightly
# because Part IX already showed that broad ordinary shifts worsen full Plik.
LOW=np.array([69.30,0.02203,0.12250,2.095e-9,0.9615,0.0490,
              0.0194,3.45,0.340,0.0470,16.05,0.0108,0.0095])
HIGH=np.array([69.75,0.02214,0.12340,2.125e-9,0.9655,0.0570,
               0.0212,3.75,0.405,0.0570,16.28,0.0145,0.0118])
KEYS=["H0","ob","oc","As","ns","tau","AF","zc","width","Dfloor","zt","Alate","Blate"]

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

# Conditional JWST/FRESCO likelihood-equivalent transfer from the standalone continuation.
# Minimal local-map unification hypothesis: g3 = A_F.
AF_REF=0.01573642223
Q_REF=0.054
AM_GRID=np.array([1.000,1.100,1.163,1.200,1.300,1.360,1.457])
E95_GRID=np.array([0.630,0.455,0.365,0.325,0.245,0.205,0.165])
E68_GRID=np.array([0.790,0.575,0.465,0.415,0.310,0.270,0.210])
def jwst(AF):
    q=Q_REF*(AF/AF_REF)
    q=float(np.clip(q,0,0.18))
    am=1+3*q
    e95=float(np.interp(am,AM_GRID,E95_GRID,left=E95_GRID[0],right=E95_GRID[-1]))
    e68=float(np.interp(am,AM_GRID,E68_GRID,left=E68_GRID[0],right=E68_GRID[-1]))
    b5=float(np.exp(25/2*(1-am**-2))); b6=float(np.exp(36/2*(1-am**-2)))
    return dict(g3_eq_AF=float(AF),Qeff=q,A_M=am,epsilon95_equiv=e95,epsilon68_equiv=e68,
                boost_5sigma=b5,boost_6sigma=b6,status="conditional likelihood-equivalent, not catalogue MCMC")

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
expansion_smg={ox:.17g},{LAM:.17g},{c['zt']:.17g},{DNT:.17g},{c['Alate']:.17g},{TAUA:.17g},{c['Blate']:.17g},{TAUB:.17g}
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
            j=jwst(c["AF"])
            rec.update(status="OK",chi2_planck_lite=pc,A_planck=Ap,delta_planck_vs_s036=pc-S036_LITE,
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
            planck_pen=rec["delta_planck_vs_s036"]
            sn_pen=max(0.,rec["sn_worst_delta"])
            bao_pen=max(0.,rec["bao_delta"])
            wl_pen=max(0.,(S8-0.815)/0.015)**2
            age_pen=max(0.,(13.45-age0)/0.1)**2
            jw_pen=0.0  # nonlinear JWST mapper coefficient is independent in this stage
            shoes_pen=rec["shoes_proxy_chi2"]/15.
            rec["screen_objective"]=float(planck_pen+2.0*sn_pen+0.15*bao_pen+0.20*wl_pen+0.10*age_pen+0.08*shoes_pen)
    except Exception as e:
        rec["error"]=repr(e)
    if rec["status"]!="OK":
        rec["error"]=rec.get("error",cp.stdout[-1200:].replace("\n"," | "))
    print("LATE266_STAGE2_POINT",json.dumps(rec,sort_keys=True),flush=True)
    # Keep only summaries to make artifact small.
    for p in OUT.glob(label+"_*"):
        try:p.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return rec

rows=[]
for c in SEEDS: rows.append(run(c))

sob=qmc.Sobol(d=len(KEYS),scramble=True,seed=20261006)
cloud=qmc.scale(sob.random_base2(m=6),LOW,HIGH) # 64 new candidates
for i,x in enumerate(cloud):
    c={k:float(v) for k,v in zip(KEYS,x)}
    c["id"]=f"fj{i:03d}"
    rows.append(run(c))

df=pd.DataFrame(rows)
df.to_csv(OUT/"late266_stage2_screen.csv",index=False)
ok=df[df.status=="OK"].copy()
if len(ok):
    # Hard pre-promotion gates keep points close enough to the established exact basin.
    gate=ok[(ok.chi2_planck_lite<=S036_LITE+0.10)&
            (ok.sn_worst_delta<=0.0)&
            (ok.bao_delta<=0.6)&
            (ok.age_t0_Gyr>=13.45)&
            (ok.cs2_min>=0.006)]
    if len(gate)<8: gate=ok.sort_values("screen_objective").head(16)
    promote=gate.sort_values(["screen_objective","S8","sn_worst_delta"]).head(8)
else:
    promote=ok
promote.to_csv(OUT/"late266_stage2_promotion.csv",index=False)
promote_cols=["id"]+KEYS+["chi2_planck_lite","delta_planck_vs_s036","sn_pp_delta","sn_u3_delta","sn_d5_delta",
                         "bao_delta","bbn_delta","S8","wl_DESY3_delta","age_t0_Gyr","shoes_proxy_chi2",
                         "jwst_A_M","jwst_epsilon95_equiv","screen_objective"]
summary={"n_total":len(rows),"n_ok":int(len(ok)),"n_shortlist":int(len(promote)),
         "promotion_rule":"Stable SN-closed points that improve the s036 screen are not winners until exact full Plik + official raw DESI + exact 3-SN + calibrated SH0ES are run.",
         "jwst_status":"reported diagnostically only; nonlinear mapper g3 is not tied to A_F in the stage-2 ranking",
         "shortlist":promote[[c for c in promote_cols if c in promote.columns]].to_dict("records")}
(OUT/"late266_stage2_summary.json").write_text(json.dumps(summary,indent=2))
print("LATE266_STAGE2_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
