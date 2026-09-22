#!/usr/bin/env python3
"""
Bayesian optimization of the late SDMC deformation using the actual target
likelihoods needed for the Part IX closure criterion.

Frozen structure: exact edge024-Q structural basin.
Varied late geometry: (A_late, tau_A, B_late, tau_B, H0).

Each evaluation computes:
  * exact parameterized hi_class background + stability
  * Planck native-lite TTTEEE+lowT+lowE+lensing screen
  * official DESI DR1 P0/P2/P4 raw full-shape profile
  * Pantheon+, Union3, DES-Y5 covariance likelihoods

The score is max(P+D fair estimate, second-smallest P+D+SN fair estimate).
Negative score means Planck+DESI closes local021 and at least 2/3 SN sets do too.
Exact covariant full-Plik/raw-DESI promotion is still required for any winner.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap, shutil
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar
from skopt import Optimizer
from skopt.space import Real

from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

from desi_fullshape_profile_utils import profile_candidate

OUT=Path("output"); RES=OUT/"snclosure_bo_rawdesi"; RES.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Frozen structural-reopt sobol009 and ordinary sector.
OB=0.022083219194622913
OC=0.12299536722293603
NS=0.9625227132590487
TAU=0.055202901571989066
Q=2.942463801247999
AS=math.exp(Q+2*TAU)/1e10
AF=0.03200255395658314
ZC=4.007449422683567
WIDTH=0.34799921004101636
D0=0.34231919445927034
DF=0.04607192634791136
LAM=18.40625
ZT=16.173189924377947
DNT=0.5

# Fair-ledger calibration to exact edge024-Q.
EDGE_PLANCK_LITE=1022.3837030216768
EDGE_PD_FAIR=-1.016724593277559
EDGE_DESI_DELTA=-13.44182058
LCDM_DESI_TOTAL=340.5324033804545
LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}

# ----- Planck-lite -----
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
        x=float(high.chi_squared(harr,A_planck=A))
        x+=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        x+=float(-2*lens.log_likelihood(dls,**lp))
        x+=((A-1.)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(pc,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

# ----- SN -----
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov mismatch {path}")
    return a.reshape(n,n)

def setup_cov(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def chi_profile(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(SNROOT/'PantheonPlus/Pantheon+SH0ES.dat',sep=r'\s+',header=0,engine='python')
ppm=pp['m_b_corr'].to_numpy(float); ppz=pp['zHD'].to_numpy(float); ppzh=pp['zHEL'].to_numpy(float)
Cpp=read_cov(SNROOT/'PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov',len(ppm))
m=ppz>0.01; ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]
pppack=setup_cov(Cpp)

p=SNROOT/'Union3/lcparam_full.txt'
cols=p.read_text().splitlines()[0].removeprefix('#').split()
un=pd.read_csv(p,sep=r'\s+',comment='#',names=cols,engine='python')
cm={c.lower():c for c in un.columns}
unm=un[cm['mb']].to_numpy(float); unz=un[cm['zcmb']].to_numpy(float)
unpack=setup_cov(read_cov(SNROOT/'Union3/mag_covmat.txt',len(unm)))

de=pd.read_csv(SNROOT/'DESY5/DES-SN5YR_HD.csv'); cm={c.lower():c for c in de.columns}
dem=de[cm['mu']].to_numpy(float); dez=de[cm['zhd']].to_numpy(float); dezh=de[cm['zhel']].to_numpy(float)
derr=de[cm['muerr_final']].to_numpy(float)
depack=setup_cov(read_cov(SNROOT/'DESY5/covsys_000.txt',len(dem))+np.diag(derr*derr))

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith('#') and re.search(r'1\s*:',l)][-1].lstrip('#').strip()
    ms=list(re.finditer(r'(\d+)\s*:\s*',hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

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

# ----- candidate CLASS run -----
ZPK=(OUT/"desi_zpk.txt").read_text().strip()

def ini(root,A,tauA,B,tauB,H0):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
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
parameters_smg={AF:.17g},{ZC:.17g},{WIDTH:.17g},{D0:.17g},1.0,{DF:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},{DNT},{A:.17g},{tauA:.17g},{B:.17g},{tauB:.17g}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
gauge=synchronous
modes=s
output=tCl,pCl,lCl,mPk,dTk,vTk
lensing=yes
l_max_scalars=3000
extra_metric_transfer_functions=yes
matter_source_in_current_gauge=no
P_k_max_h/Mpc=2.0
z_pk={ZPK}
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

history=[]
def cleanup(prefix):
    for p in OUT.glob(prefix+"_*"):
        try:p.unlink()
        except:pass
    ip=OUT/(prefix+".ini")
    try:ip.unlink()
    except:pass

def evaluate(x,tag):
    A,tauA,B,tauB,H0=map(float,x)
    prefix=f"bo_{tag}"
    root=str(OUT/(prefix+"_"))
    ip=OUT/(prefix+".ini"); ip.write_text(ini(root,A,tauA,B,tauB,H0))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=420)
    rec=dict(id=tag,A=A,tauA=tauA,B=B,tauB=tauB,H0=H0,status="FAIL",objective=100.0)
    bgp=OUT/f"{prefix}_00_background.dat"; clp=OUT/f"{prefix}_00_cl_lensed.dat"
    if cp.returncode==0 and bgp.exists() and clp.exists():
        try:
            bg=table(bgp)
            rec.update(min_D=float(bg["kin (D)"].min()),
                       min_cs2=float(bg["c_s^2"].min()),
                       max_cs2=float(bg["c_s^2"].max()),
                       max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
            stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1
            rec["stable_subluminal"]=bool(stable)
            if stable:
                pchi,Ap=pscore(clp); rec["chi2_planck_lite"]=pchi; rec["A_planck"]=Ap
                sn=sn_scores(bg)
                rec.update({f"sn_{k}":v for k,v in sn.items()})
                for k,v in sn.items():rec[f"sn_delta_{k}"]=v-LOCAL_SN[k]
                rows,desi_total,m=profile_candidate(prefix,H0,OB,OC,NS,label=f"BO_{tag}")
                desi_delta=float(desi_total-LCDM_DESI_TOTAL)
                rec["desi_chi2"]=float(desi_total); rec["desi_delta"]=desi_delta
                rec["delta_planck_vs_edge"]=pchi-EDGE_PLANCK_LITE
                rec["delta_desi_vs_edge"]=desi_delta-EDGE_DESI_DELTA
                pdfair=EDGE_PD_FAIR+rec["delta_planck_vs_edge"]+rec["delta_desi_vs_edge"]
                rec["pd_fair_est"]=float(pdfair)
                joints=[pdfair+rec["sn_delta_pantheonplus"],pdfair+rec["sn_delta_union3"],pdfair+rec["sn_delta_desy5"]]
                rec["joint_pp_est"],rec["joint_u3_est"],rec["joint_d5_est"]=map(float,joints)
                rec["second_joint_est"]=float(sorted(joints)[1])
                rec["n_sn_closed_est"]=int(sum(v<0 for v in joints))
                rec["objective"]=float(max(pdfair,rec["second_joint_est"]))
                rec["status"]="OK"
        except Exception as e:
            rec["error"]=repr(e)
    if rec["status"]!="OK":
        rec["error"]=rec.get("error",cp.stdout[-800:].replace("\n"," | "))
    print("BO_RAWDESI_POINT",json.dumps(rec,sort_keys=True),flush=True)
    history.append(rec)
    pd.DataFrame(history).to_csv(RES/"bo_history.csv",index=False)
    cleanup(prefix)
    return rec["objective"]

space=[
 Real(0.005,0.035,name="A"),
 Real(0.12,0.55,name="tauA"),
 Real(0.006,0.028,name="B"),
 Real(0.75,2.60,name="tauB"),
 Real(69.4,71.10,name="H0"),
]
opt=Optimizer(space,base_estimator="GP",acq_func="EI",acq_optimizer="sampling",
              random_state=923,initial_point_generator="sobol",n_initial_points=8)

# Seed physically important known points.
seeds=[
 [0.01105624999,0.25,0.01951933685,1.5,70.79187943335671], # edge024 geometry
 [0.018589897081255913,0.25,0.013815921754576268,1.5,69.85232012766193], # sobol223
 [0.0208165963771753,0.25,0.012095115539617837,1.5,69.57700432301499], # excellent SN local refine
]
for i,x in enumerate(seeds):
    y=evaluate(x,f"seed{i:02d}"); opt.tell(x,y)

N_TOTAL=28
for i in range(N_TOTAL-len(seeds)):
    x=opt.ask()
    y=evaluate(x,f"bo{i:02d}")
    opt.tell(x,y)

ok=pd.DataFrame(history)
ok=ok[ok.status=="OK"].sort_values(["objective","second_joint_est","pd_fair_est"])
best=ok.head(10).to_dict("records")
summary={"frozen_structure":{"AF":AF,"ZC":ZC,"width":WIDTH,"D_floor":DF,"Q":Q},
         "n_total":len(history),"n_ok":int(len(ok)),
         "best":best[0] if best else None,"top10":best,
         "criterion":"objective<0 => P+D fair<0 and at least 2 of 3 SN-added fair totals<0"}
(RES/"bo_summary.json").write_text(json.dumps(summary,indent=2))
print("BO_RAWDESI_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
