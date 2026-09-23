#!/usr/bin/env python3
"""
Exact-joint Bayesian controller for the final SDMC SN-closure search.

Latent coordinates:
  u = 0 -> sobol107 late background; u = 1 -> coupled sobol036 background
  v = 0 -> bo002 structure;        v = 1 -> seed015 structure

Training targets are exact promoted results only:
  full Planck Plik + raw DESI full shape + exact SN covariance penalties.

Goal:
  P+D fair < 0 and at least two of
  (P+D+Pantheon+), (P+D+Union3), (P+D+DES-Y5) < 0.
"""
from pathlib import Path
import json, numpy as np, pandas as pd
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

OUT=Path("output/fulljoint_exact_bo_seed015"); OUT.mkdir(parents=True,exist_ok=True)
LOCAL_PD=-31.15223

# Exact promoted observations only.
obs=[
 {"id":"bo002","u":0.0,"v":0.0,"planck":-19.92492,"desi":-11.90251,
  "sn_pp":1.6312680542469025,"sn_u3":1.3954445793788182,"sn_d5":3.27561803907156},
 {"id":"sob036_bo002","u":1.0,"v":0.0,"planck":-19.581007560230773,"desi":-11.193284541589946,
  "sn_pp":0.5575426062569022,"sn_u3":0.5355814192298567,"sn_d5":1.2224958017468452},
 {"id":"sob036_seed015","u":1.0,"v":1.0,"planck":-19.930025292910614,"desi":-11.168825504830693,
  "sn_pp":0.5575426062569022,"sn_u3":0.5355814192298567,"sn_d5":1.2224958017468452},
 {"id":"exactbo000","u":0.9725,"v":0.97,"planck":-19.94159805526033,"desi":-11.190121031232081,
  "sn_pp":0.585484700743109,"sn_u3":0.5586358160071541,"sn_d5":1.2772395089268684},
]
for r in obs:
    r["pd_fair"]=r["planck"]+r["desi"]-LOCAL_PD
    js=[r["pd_fair"]+r["sn_pp"],r["pd_fair"]+r["sn_u3"],r["pd_fair"]+r["sn_d5"]]
    r["joint_pp"],r["joint_u3"],r["joint_d5"]=js
    r["second_sn_joint"]=sorted(js)[1]
    r["goal"]=max(r["pd_fair"],r["second_sn_joint"])

df=pd.DataFrame(obs)
X=df[["u","v"]].to_numpy(float); y=df["goal"].to_numpy(float)
kernel=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=[0.5,0.25],length_scale_bounds=(0.05,5.0),nu=2.5)+WhiteKernel(1e-8,(1e-10,1e-3))
gp=GaussianProcessRegressor(kernel=kernel,normalize_y=True,n_restarts_optimizer=12,random_state=15015)
gp.fit(X,y)

ug=np.linspace(0.50,1.00,201)
vg=np.linspace(0.75,1.00,151)
G=np.array([(u,v) for u in ug for v in vg],float)
mu,std=gp.predict(G,return_std=True)
best=float(y.min()); xi=0.005
imp=best-mu-xi
z=np.divide(imp,std,out=np.zeros_like(imp),where=std>1e-12)
ei=imp*norm.cdf(z)+std*norm.pdf(z)

order=np.argsort(ei)[::-1]
chosen=[]
for j in order:
    x=G[j]
    if np.min(np.sqrt(np.sum((X-x)**2,axis=1)))<0.04: continue
    if any(np.linalg.norm(x-c)<0.08 for c in chosen): continue
    chosen.append(x)
    if len(chosen)==4: break

A107=0.0013601897284388; B107=0.0147392870020121; H107=70.09917331933976
A036=0.011249459300190211; B036=0.01280436261370778; H036=69.80878993044607
T036=0.014789129979908472
BO=np.array([0.024290704212870225,3.6096407580714724,0.3631213944496205,0.04602317962943721])
SEED=np.array([0.0235834146537818,3.4274108000472188,0.3689727572351694,0.050407338682562114])

props=[]
for i,x in enumerate(chosen):
    u,v=map(float,x)
    bg=[A107+u*(A036-A107),B107+u*(B036-B107),H107+u*(H036-H107),u*T036]
    st=BO+v*(SEED-BO)
    j=np.argmin(np.sum((G-x)**2,axis=1))
    props.append({"id":f"exactbo{i:03d}","u":u,"v":v,
      "A":bg[0],"B":bg[1],"H0":bg[2],"t_ord":bg[3],
      "AF":float(st[0]),"ZC":float(st[1]),"width":float(st[2]),"D_floor":float(st[3]),
      "gp_mu_goal":float(mu[j]),"gp_sigma":float(std[j]),"expected_improvement":float(ei[j])})

df.to_csv(OUT/"training_exact.csv",index=False)
pd.DataFrame(props).to_csv(OUT/"proposals.csv",index=False)
summary={"kernel":str(gp.kernel_),"best_exact_goal":best,"training":obs,"proposals":props,
 "rule":"Every proposal must be promoted through exact full Plik, raw DESI full shape, and exact SN covariance likelihoods before it can update the GP."}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("FULLJOINT_EXACT_BO",json.dumps(summary,sort_keys=True),flush=True)
