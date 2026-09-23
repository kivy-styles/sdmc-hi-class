#!/usr/bin/env python3
"""
Hard-gate exact-informed Bayesian controller for SDMC.

Target:
  delta_Planck(local021) < 0
  delta_DESI(local021) < 0
  at least two exact SN covariance deltas(local021) < 0

Only fully exact promoted points train the GP.
Features are the four active coordinates in the current bridge:
  H0, late-lobe A, late-lobe B, structural AF.
"""
from pathlib import Path
import json, numpy as np, pandas as pd
from scipy.stats import norm, qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler

OUT=Path("output/hardgate_trilikelihood_bo_round1"); OUT.mkdir(parents=True,exist_ok=True)

P0=2781.4370287015954
D0=331.93726232610265
SN0={"pp":1406.214941392187,"u3":28.78345144550258,"d5":1650.6363433003426}

obs=[
 {"id":"exactbo000","H0":69.81677547364065,"A":0.010977504386967048,"B":0.012857573034386149,"AF":0.023604633340554453,
  "planck":2784.05650677952,"desi":329.34250909183334,"pp":1406.8001879565418,"u3":29.341775108259753,"d5":1651.9126336723566},
 {"id":"exactbo004","H0":69.82911676766862,"A":0.010557210430167612,"B":0.012939807320889084,"AF":0.023632917504307768,
  "planck":2784.0357790433236,"desi":329.30971730204203,"pp":1406.8435462792404,"u3":29.377469697079505,"d5":1651.997425943613},
 {"id":"exactbo2000","H0":69.80878993044607,"A":0.011249459300190211,"B":0.01280436261370778,"AF":0.023760237043553907,
  "planck":2784.0641657880965,"desi":329.36357787562383,"pp":1406.7722473018803,"u3":29.318722501484444,"d5":1651.857893191278},
 {"id":"late298","H0":69.58899342644959,"A":0.013419163235463202,"B":0.011442640743218362,"AF":0.023604633340554453,
  "planck":2784.4077241594355,"desi":330.1382662686478,"pp":1406.21248516161,"u3":28.839521281566704,"d5":1650.7453144714236},
 {"id":"late266","H0":69.45160505326464,"A":0.013036667098291217,"B":0.010544837576337158,"AF":0.023604633340554453,
  "planck":2784.6643534348414,"desi":330.69550645610434,"pp":1405.935751490295,"u3":28.592742397959228,"d5":1650.182109862566}
]

for r in obs:
    r["dP"]=r["planck"]-P0
    r["dD"]=r["desi"]-D0
    r["dPP"]=r["pp"]-SN0["pp"]
    r["dU3"]=r["u3"]-SN0["u3"]
    r["dD5"]=r["d5"]-SN0["d5"]
    r["dSN2"]=sorted([r["dPP"],r["dU3"],r["dD5"]])[1]
    r["hard_goal"]=max(r["dP"],r["dD"],r["dSN2"])
    r["closed"]=bool(r["dP"]<0 and r["dD"]<0 and r["dSN2"]<0)

df=pd.DataFrame(obs)
features=["H0","A","B","AF"]
X=df[features].to_numpy(float)
sc=StandardScaler().fit(X)
Xs=sc.transform(X)

def gpfit(y,seed):
    ker=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=np.ones(4),length_scale_bounds=(0.08,8.0),nu=2.5)+WhiteKernel(1e-7,(1e-10,1e-2))
    gp=GaussianProcessRegressor(kernel=ker,normalize_y=True,n_restarts_optimizer=20,random_state=seed)
    gp.fit(Xs,np.asarray(y,float))
    return gp

gpP=gpfit(df.dP,901)
gpD=gpfit(df.dD,902)
gpS=gpfit(df.dSN2,903)

# Search only the local bridge supported by exact points, expanded modestly.
mins=X.min(0); maxs=X.max(0); span=np.maximum(maxs-mins,[0.08,0.0008,0.0008,0.00008])
lo=mins-0.20*span; hi=maxs+0.20*span
# Keep physically sensible narrow bounds.
lo=np.maximum(lo,[69.30,0.0098,0.0098,0.02350])
hi=np.minimum(hi,[69.90,0.0142,0.0134,0.02386])

sob=qmc.Sobol(d=4,scramble=True,seed=23092301)
G=lo+(hi-lo)*sob.random_base2(m=15)
Gs=sc.transform(G)

mP,sP=gpP.predict(Gs,return_std=True)
mD,sD=gpD.predict(Gs,return_std=True)
mS,sS=gpS.predict(Gs,return_std=True)

# Probability each hard gate is satisfied.
pP=norm.cdf((0-mP)/np.maximum(sP,1e-9))
pD=norm.cdf((0-mD)/np.maximum(sD,1e-9))
pS=norm.cdf((0-mS)/np.maximum(sS,1e-9))
pfeas=pP*pD*pS

# Improvement in the worst hard gate, with feasibility-weighted exploration.
train_goal=np.maximum.reduce([df.dP.to_numpy(),df.dD.to_numpy(),df.dSN2.to_numpy()])
best=float(train_goal.min())
mG=np.maximum.reduce([mP,mD,mS])
# Conservative uncertainty for max-objective.
sG=np.maximum.reduce([sP,sD,sS])
imp=best-mG-0.02
z=imp/np.maximum(sG,1e-9)
ei=imp*norm.cdf(z)+sG*norm.pdf(z)
acq=ei*(0.15+0.85*pfeas)+0.08*pfeas+0.02*(sP+sD+sS)

order=np.argsort(acq)[::-1]
chosen=[]
for idx in order:
    x=G[idx]
    # exclude points too close to exact observations or one another
    dx=(x-X)/span
    if np.min(np.sqrt(np.sum(dx*dx,axis=1)))<0.10: continue
    if any(np.linalg.norm((x-c)/span)<0.14 for c in chosen): continue
    chosen.append(x)
    if len(chosen)==6: break

props=[]
for i,x in enumerate(chosen):
    j=int(np.argmin(np.sum((G-x)**2,axis=1)))
    props.append({
      "id":f"hardbo{i:03d}",
      **{k:float(v) for k,v in zip(features,x)},
      "gp_dP":float(mP[j]),"gp_dD":float(mD[j]),"gp_dSN2":float(mS[j]),
      "p_planck_closed":float(pP[j]),"p_desi_closed":float(pD[j]),"p_two_sn_closed":float(pS[j]),
      "p_all_hard_gates":float(pfeas[j]),"acquisition":float(acq[j])
    })

df.to_csv(OUT/"training_exact.csv",index=False)
pd.DataFrame(props).to_csv(OUT/"proposals.csv",index=False)
summary={
 "baselines":{"planck_local021":P0,"desi_local021":D0,"sn_local021":SN0},
 "features":features,
 "bounds":{"lo":lo.tolist(),"hi":hi.tolist()},
 "best_training_hard_goal":best,
 "training":obs,
 "kernels":{"P":str(gpP.kernel_),"D":str(gpD.kernel_),"SN2":str(gpS.kernel_)},
 "proposals":props,
 "rule":"A candidate counts only if exact full-Plik dP<0, raw DESI dD<0, and the second-smallest exact SN delta<0. Every proposal must be exact-promoted before updating the GP."
}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("HARDGATE_TRILIKELIHOOD_BO",json.dumps(summary,sort_keys=True),flush=True)
