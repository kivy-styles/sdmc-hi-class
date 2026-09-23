#!/usr/bin/env python3
"""
Bayesian structural recovery around the coupled sobol036 SN-improved background.

Low-z geometry is frozen at the SN-friendly point:
  A_late=0.0013601897284388
  B_late=0.0147392870020121
  H0=69.83619736097754

Primordial amplitude uses the tiny Q improvement found in the preceding
ordinary pass:
  Q=2.9434975776090635

Reopen only perturbation/kinetic structure:
  A_F, z_c, width, D_floor

Because the explicit sdmc_full background is fixed, these directions are
intended to recover CMB likelihood without undoing the SN distance gain.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc, norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure107_seed015_sobol127_structural_bo"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Frozen SN-friendly background + ordinary sector
H0=69.83619736097754
OB=0.022092359690331246; OC=0.12294585637584288
NS=0.9626800328202773; TAU=0.055202901571989066
Q=2.9434975776090635; AS=math.exp(Q+2*TAU)/1e10
AL=0.01264600746287033; BL=0.012972991137765347
LAM=18.40625; ZT=16.173189924377947; DN=.5; TAUA=.25; TAUB=1.5
D0=0.34231919445927034

# edge024 structural center
AF0=0.0235834146537818; ZC0=3.4274108000472188
W0=0.3689727572351694; DF0=0.050407338682562114

# Fair-ledger constants used only for promotion ranking.
EDGE_P_LITE=1021.0852360745247
EDGE_PD_FAIR=-0.21418646229219362
EDGE_BAO=13.215120771731913
NEW_BAO=13.215120771731913
SN_PP=0.5864693196490407
SN_U3=0.5641153095639311
SN_D5=1.278687983751297

LOW=np.array([0.0226,3.36,0.345,0.048])
HIGH=np.array([0.0245,3.48,0.395,0.056])

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
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
    a=np.loadtxt(path); e=a[:,0].astype(int); n=int(e.max())+1; cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[e]=a[:,1]*cv; ee[e]=a[:,2]*cv; te[e]=a[:,3]*cv
    L=e.astype(float); pp[e]=a[:,5]*L*(L+1)
    return np.column_stack([e,tt[e],te[e],ee[e]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    ha,d=load_cls(path)
    def f(A):
        x=float(high.chi_squared(ha,A_planck=A))
        x+=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        x+=float(-2*lens.log_likelihood(d,**lp))
        x+=((A-1)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def ini(root,AF,ZC,W,DF):
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
parameters_smg={AF:.17g},{ZC:.17g},{W:.17g},{D0:.17g},1.0,{DF:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},{DN},{AL},{TAUA},{BL},{TAUB}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
write background=yes
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

def cand(tag,AF,ZC,W,DF):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,AF,ZC,W,DF))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    r={"id":tag,"AF":float(AF),"ZC":float(ZC),"width":float(W),"D_floor":float(DF),"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        r.update(min_D=float(bg["kin (D)"].min()),
                 min_cs2=float(bg["c_s^2"].min()),
                 max_cs2=float(bg["c_s^2"].max()),
                 max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1
        r["stable_subluminal"]=bool(stable)
        if stable:
            pc,Ap=pscore(clp); r.update(status="OK",chi2_planck=pc,A_planck=Ap)
            r["delta_planck_vs_edge024"]=pc-EDGE_P_LITE
            r["pd_proxy"]=EDGE_PD_FAIR+r["delta_planck_vs_edge024"]+.25*(NEW_BAO-EDGE_BAO)
            r["joint_pp_proxy"]=r["pd_proxy"]+SN_PP
            r["joint_u3_proxy"]=r["pd_proxy"]+SN_U3
            r["joint_d5_proxy"]=r["pd_proxy"]+SN_D5
            js=[r["joint_pp_proxy"],r["joint_u3_proxy"],r["joint_d5_proxy"]]
            r["second_joint_proxy"]=sorted(js)[1]
            r["n_sn_closed_proxy"]=sum(x<0 for x in js)
            r["goal_score"]=max(r["pd_proxy"],r["second_joint_proxy"])
    if r["status"]!="OK": r["error"]=cp.stdout[-500:].replace("\n"," | ")
    print("SN127BO_POINT",json.dumps(r,sort_keys=True),flush=True)
    return r

# Seed the GP with the center, known successful structures, coordinate anchors,
# and a small Sobol design.  The known points came from the previous structural
# searches and are valuable priors even though the late background has changed.
pts=[
 ("center",AF0,ZC0,W0,DF0),
 ("known_sobol009",0.024221895329654217,4.145109392143786,0.3370092310383916,0.04453643597662449),
 ("known_s110",0.026252536563202738,3.79610941009596,0.3229390993900597,0.046828035046346486),
 ("known_s013",0.02582530555780977,3.8840197445824742,0.3657897564768791,0.043821363696828486),
 ("known_s029",0.02663535825908184,3.7860057394951583,0.4024500951357186,0.05104022843576968),
 ("known_sobol006",0.02559563393890858,3.8428389857523144,0.388143997490406,0.04659694366157055),
 ("afm",AF0-.004,ZC0,W0,DF0),("afp",AF0+.004,ZC0,W0,DF0),
 ("zcm",AF0,ZC0-.20,W0,DF0),("zcp",AF0,ZC0+.20,W0,DF0),
 ("wm",AF0,ZC0,W0-.035,DF0),("wp",AF0,ZC0,W0+.035,DF0),
 ("dfm",AF0,ZC0,W0,DF0-.004),("dfp",AF0,ZC0,W0,DF0+.004)
]
sam=qmc.Sobol(d=4,scramble=True,seed=127411)
for i,p in enumerate(qmc.scale(sam.random_base2(m=4),LOW,HIGH)):
    pts.append((f"seed{i:03d}",*map(float,p)))

rows=[cand(*p) for p in pts]

def fit_gp(rows):
    ok=[r for r in rows if r.get("status")=="OK"]
    X=np.array([[r["AF"],r["ZC"],r["width"],r["D_floor"]] for r in ok],float)
    y=np.array([r["chi2_planck"] for r in ok],float)
    Xn=(X-LOW)/(HIGH-LOW)
    ker=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=np.ones(4)*0.25,length_scale_bounds=(0.03,3.0),nu=2.5)+WhiteKernel(1e-6,(1e-9,1e-2))
    gp=GaussianProcessRegressor(kernel=ker,normalize_y=True,n_restarts_optimizer=3,random_state=127412)
    gp.fit(Xn,y)
    return gp,float(y.min())

# Sequential expected-improvement Bayesian optimization.
rng=np.random.default_rng(127413)
for it in range(20):
    gp,ybest=fit_gp(rows)
    pool=rng.uniform(LOW,HIGH,size=(6000,4))
    # Add a local cloud around the current best to improve exploitation.
    ok=[r for r in rows if r.get("status")=="OK"]
    rb=min(ok,key=lambda r:r["chi2_planck"])
    xb=np.array([rb["AF"],rb["ZC"],rb["width"],rb["D_floor"]])
    local=xb+rng.normal(size=(2000,4))*0.06*(HIGH-LOW)
    local=np.clip(local,LOW,HIGH)
    pool=np.vstack([pool,local])
    xn=(pool-LOW)/(HIGH-LOW)
    mu,std=gp.predict(xn,return_std=True)
    imp=ybest-mu-0.002
    z=np.divide(imp,std,out=np.zeros_like(imp),where=std>1e-12)
    ei=imp*norm.cdf(z)+std*norm.pdf(z)
    # Avoid re-evaluating points too near existing samples.
    Xold=np.array([[r["AF"],r["ZC"],r["width"],r["D_floor"]] for r in rows],float)
    Xoldn=(Xold-LOW)/(HIGH-LOW)
    for j in np.argsort(ei)[::-1]:
        if np.min(np.sqrt(np.sum((Xoldn-xn[j])**2,axis=1)))>0.012:
            x=pool[j]; break
    rows.append(cand(f"bo{it:03d}",*map(float,x)))

df=pd.DataFrame(rows); df.to_csv(OUT/"snclosure107_seed015_sobol127_structural_bo.csv",index=False)
ok=df[df.status=="OK"].sort_values(["chi2_planck","goal_score"])
best=ok.head(15).to_dict("records")
summary={"background":{"H0":H0,"A":AL,"B":BL},"Q":Q,"n_ok":int(len(ok)),
         "best":best,"target":"sobol127 structural refinement; seek >=0.38 exact-equivalent Planck gain while preserving stability"}
(OUT/"snclosure107_seed015_sobol127_structural_bo_summary.json").write_text(json.dumps(summary,indent=2))
print("SN127BO_BEST",json.dumps(summary,sort_keys=True),flush=True)
