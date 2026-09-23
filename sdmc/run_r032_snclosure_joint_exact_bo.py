#!/usr/bin/env python3
"""
Exact-informed joint Bayesian search round 2 for the SDMC SN-closure problem.

Search manifold:
  s_bg    : interpolation/extrapolation bo002 -> sobol036 late background
  t_ord   : ordinary-sector displacement along the empirically tested O4 direction
  u_struct: interpolation bo002 -> seed015 perturbation structure

The Gaussian process is trained ONLY on measured exact fair Planck+raw-DESI
results.  Every proposal gets an exact covariance SN calculation and an exact
CLASS stability gate.  Acquisition winners are then promoted to the existing
full-Plik + official raw-DESI workflows.  This is deliberately multi-fidelity:
we do not spend a full Plik + raw DESI evaluation on hundreds of exploratory
points, but every point admitted to the exact training set came from those
likelihoods.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

OUT=Path("output/snclosure_joint_exact_bo_r2"); OUT.mkdir(parents=True,exist_ok=True)
TMP=OUT/"tmp"; TMP.mkdir(exist_ok=True)
SNROOT=Path("sn_data")
OR=4.17998772e-5

# Exact bo002 anchor.
A0=0.0013601897284388; B0=0.0147392870020121; H00=70.09917331933976
AF0=0.024290704212870225; ZC0=3.6096407580714724
W0=0.3631213944496205; DF0=0.04602317962943721

# Exact sobol036 + seed015 anchor.
A1=0.011249459300190211; B1=0.01280436261370778; H01=69.80878993044607
AF1=0.0235834146537818; ZC1=3.4274108000472188
W1=0.3689727572351694; DF1=0.050407338682562114

# Ordinary-sector interpolation: bo002 -> old O4 direction.
OB0=0.022083219194622913; OC0=0.12299536722293603
NS0=0.9625227132590487; AS0=2.119721074504234e-9
OB1=0.022215084896497428; OC1=0.12228109714277088
NS1=0.9647922897702084; AS1=2.120753981008663e-9
TAU=0.055202901571989066

LAM=18.40625; ZT=16.173189924377947; DNT=0.5; TAUA=0.25; TAUB=1.5
D0=0.34231919445927034

LOCAL_PD=-31.15223
LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}

# Exact full-Plik + raw-DESI training points already measured.
# y is fair P+D = dchi2_Planck + dchi2_DESI - LOCAL_PD.
TRAIN=[
 ("bo002",    0.0,  0.00, 0.0, -19.924920000000000, -11.902510000000000),
 ("tm10",     0.0, -0.10, 0.0, -19.843057215940917, -11.864786320198448),
 ("tm20",     0.0, -0.20, 0.0, -19.352524524435920, -11.826884658642768),
 ("t10",      0.0,  0.10, 0.0, -19.626250860538220, -11.937905571882538),
 ("line25",   0.0,  0.25, 0.0, -18.507734480853742, -11.990634046394064),
 ("sobol036", 1.0,  0.014789129979908472, 0.0, -19.581007560230773, -11.193284541589946),
 ("seed015",  1.0,  0.014789129979908472, 1.0, -19.930025292910614, -11.168825504830693),
]
LOW=np.array([0.30,-0.10,0.45],float)
HIGH=np.array([1.20, 0.10,1.18],float)

def interp(x0,x1,q): return x0+q*(x1-x0)

def params(s,t,u):
    A=interp(A0,A1,s); B=interp(B0,B1,s); H0=interp(H00,H01,s)
    ob=interp(OB0,OB1,t); oc=interp(OC0,OC1,t)
    ns=interp(NS0,NS1,t); As=interp(AS0,AS1,t)
    AF=interp(AF0,AF1,u); ZC=interp(ZC0,ZC1,u)
    W=interp(W0,W1,u); DF=interp(DF0,DF1,u)
    h=H0/100.; ox=1.-(ob+oc+OR)/(h*h)
    return dict(A=A,B=B,H0=H0,ob=ob,oc=oc,ns=ns,As=As,
                AF=AF,ZC=ZC,width=W,D_floor=DF,Ox=ox)

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
    q=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(q))

def chi_profile(pack,r):
    cf,one,den=pack
    y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

# Exact SN covariance products.
pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm))
m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]
pppack=setup_cov(Cpp)

p=SNROOT/"Union3/lcparam_full.txt"
cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float)
unpack=setup_cov(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv")
cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
depack=setup_cov(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

def mu_bg(bg,z,zh):
    zz=bg.z.to_numpy(); dm=bg["comov. dist."].to_numpy(); o=np.argsort(zz)
    DM=np.interp(z,zz[o],dm[o]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def sn_scores(bg):
    return {
      "pantheonplus":chi_profile(pppack,ppm-mu_bg(bg,ppz,ppzh)),
      "union3":chi_profile(unpack,unm-mu_bg(bg,unz,unz)),
      "desy5":chi_profile(depack,dem-mu_bg(bg,dez,dezh)),
    }

def ini_text(root,p):
    return textwrap.dedent(f"""\
H0={p['H0']:.17g}
omega_b={p['ob']:.17g}
omega_cdm={p['oc']:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={p['As']:.17e}
n_s={p['ns']:.17g}
tau_reio={TAU:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={p['AF']:.17g},{p['ZC']:.17g},{p['width']:.17g},{D0:.17g},1.0,{p['D_floor']:.17g}
expansion_model=sdmc_full
expansion_smg={p['Ox']:.17g},{LAM},{ZT},{DNT},{p['A']:.17g},{TAUA},{p['B']:.17g},{TAUB}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
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

def evaluate_background(tag,s,t,u):
    p=params(s,t,u); root=str(TMP/(tag+"_")); ip=TMP/(tag+".ini")
    ip.write_text(ini_text(root,p))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=120)
    r={"id":tag,"s_bg":float(s),"t_ord":float(t),"u_struct":float(u),**p,"status":"FAIL"}
    bgp=Path(root+"00_background.dat")
    if cp.returncode==0 and bgp.exists():
        bg=table(bgp)
        r["min_D"]=float(bg["kin (D)"].min())
        r["min_cs2"]=float(bg["c_s^2"].min())
        r["max_cs2"]=float(bg["c_s^2"].max())
        r["stable_subluminal"]=bool(r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1)
        if r["stable_subluminal"]:
            sn=sn_scores(bg)
            for k,v in sn.items():
                r[f"sn_{k}"]=v
                r[f"sn_delta_{k}"]=v-LOCAL_SN[k]
            r["status"]="OK"
    if r["status"]!="OK": r["error"]=cp.stdout[-500:].replace("\n"," | ")
    for q in TMP.glob(tag+"_*"):
        try:q.unlink()
        except:pass
    try:ip.unlink()
    except:pass
    return r

# GP on the measured exact P+D fair surface.
X=[]; y=[]; train_rows=[]
for name,s,t,u,pchi,dchi in TRAIN:
    pd_fair=float(pchi+dchi-LOCAL_PD)
    X.append([s,t,u]); y.append(pd_fair)
    train_rows.append(dict(id=name,s_bg=s,t_ord=t,u_struct=u,
                           delta_planck=pchi,delta_desi=dchi,pd_fair=pd_fair))
X=np.asarray(X,float); y=np.asarray(y,float)
Xn=(X-LOW)/(HIGH-LOW)
kernel=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=np.ones(3)*0.35,
        length_scale_bounds=(0.05,5.0),nu=2.5)+WhiteKernel(1e-5,(1e-9,1e-1))
gp=GaussianProcessRegressor(kernel=kernel,normalize_y=True,n_restarts_optimizer=8,random_state=1503611)
gp.fit(Xn,y)

sob=qmc.Sobol(d=3,scramble=True,seed=1503612)
pool=qmc.scale(sob.random_base2(m=8),LOW,HIGH)  # 256 exact-SN/stability proposals
rows=[]
for i,(s,t,u) in enumerate(pool):
    r=evaluate_background(f"joint{i:03d}",s,t,u)
    if r["status"]=="OK":
        xn=(np.array([[s,t,u]])-LOW)/(HIGH-LOW)
        mu,std=gp.predict(xn,return_std=True)
        mu=float(mu[0]); std=float(std[0])
        r["pd_exact_gp_mean"]=mu; r["pd_exact_gp_std"]=std
        js=[mu+r["sn_delta_pantheonplus"],mu+r["sn_delta_union3"],mu+r["sn_delta_desy5"]]
        r["joint_pp_mean"]=js[0]; r["joint_u3_mean"]=js[1]; r["joint_d5_mean"]=js[2]
        r["second_joint_mean"]=sorted(js)[1]
        # Lower confidence bound drives exploration toward plausible exact closure.
        beta=1.25
        r["pd_lcb"]=mu-beta*std
        r["second_joint_lcb"]=r["second_joint_mean"]-beta*std
        r["goal_mean"]=max(mu,r["second_joint_mean"])
        r["acquisition_lcb"]=max(r["pd_lcb"],r["second_joint_lcb"])
    rows.append(r)
    if i%32==0: print("JOINTBO_PROGRESS",i,json.dumps(r,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"joint_exact_bo_pool.csv",index=False)
ok=df[df.status.eq("OK")].sort_values(["acquisition_lcb","goal_mean"]).copy()

# Diversity-filter the promotion list.
chosen=[]; chosen_x=[]
for _,r in ok.iterrows():
    x=np.array([r.s_bg,r.t_ord,r.u_struct],float)
    xn=(x-LOW)/(HIGH-LOW)
    if all(np.linalg.norm(xn-z)>0.10 for z in chosen_x):
        chosen.append(r.to_dict()); chosen_x.append(xn)
    if len(chosen)>=8: break

summary={
 "training":train_rows,
 "kernel":str(gp.kernel_),
 "domain":{"low":LOW.tolist(),"high":HIGH.tolist()},
 "n_pool":int(len(df)),"n_stable":int(len(ok)),
 "promote":chosen,
 "stopping_rule":"Exact promotion succeeds when fair P+D < 0 and the second-smallest of [P+D+PP, P+D+U3, P+D+D5] is < 0.",
 "method":"GP trained on measured exact full-Plik+raw-DESI fair scores; exact SN covariance and CLASS stability evaluated for every acquisition candidate."
}
(OUT/"joint_exact_bo_summary.json").write_text(json.dumps(summary,indent=2))
print("JOINTBO_BEST",json.dumps(chosen,sort_keys=True),flush=True)
