#!/usr/bin/env python3
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
from scipy.stats import norm, qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

OUT=Path("output/exact_constrained_bo_h0ab"); OUT.mkdir(parents=True,exist_ok=True)

# Strict local021 baselines, regenerated with the same exact likelihood codes.
P0=2781.4370287015954
D0=331.93726232610265
SN0=np.array([1406.214941392187,28.78345144550258,1650.6363433003426])

# All observations below share the exactbo000 structure/ordinary sector and differ
# only in the late-background coordinates H0,A,B. Every target is an exact promoted
# full-Plik/raw-DESI/exact-SN result.
raw=[
 dict(id="exactbo000",H0=69.81677547364065,A=0.010977504386967048,B=0.012857573034386149,
      planck=2784.05650677952,desi=329.34250909183334,
      pp=1406.8001879565418,u3=29.341775108259753,d5=1651.9126336723566),
 dict(id="hybrid167",H0=69.51717766530813,A=0.008974749634042383,B=0.011274873235728592,
      planck=2784.9246844210033,desi=331.891151096234,
      pp=1406.2192904320545,u3=28.83534829504788,d5=1650.7610308378935),
 dict(id="hybrid043",H0=69.486625116691,A=0.007934502487070859,B=0.01100946151232347,
      planck=2784.92363121286,desi=332.2669697694018,
      pp=1406.1856789663434,u3=28.80288930531242,d5=1650.692088894546),
 dict(id="hybrid063",H0=69.44595765553416,A=0.010335846113041045,B=0.010714279787614942,
      planck=2784.790105139299,desi=332.68654218492713,
      pp=1406.0198865132406,u3=28.662211059985566,d5=1650.3558155372739),
 dict(id="late298",H0=69.58899342644959,A=0.013419163235463202,B=0.011442640743218362,
      planck=2784.4077241594355,desi=330.1382662686478,
      pp=1406.21248516161,u3=28.839521281566704,d5=1650.7453144714236),
 dict(id="late266",H0=69.45160505326464,A=0.013036667098291217,B=0.010544837576337158,
      planck=2784.6643534348414,desi=330.69550645610434,
      pp=1405.935751490295,u3=28.592742397959228,d5=1650.182109862566),
]
for r in raw:
    r["dp"]=r["planck"]-P0
    r["dd"]=r["desi"]-D0
    snd=np.array([r["pp"],r["u3"],r["d5"]])-SN0
    r["dpp"],r["du3"],r["dd5"]=map(float,snd)
    r["sn2"]=float(np.sort(snd)[1])
    r["feasible"]=bool(r["dd"]<0 and r["sn2"]<0)
    r["strict_goal"]=float(max(r["dp"],r["dd"],r["sn2"]))

df=pd.DataFrame(raw)
X=df[["H0","A","B"]].to_numpy(float)
lo=np.array([69.30,0.0065,0.0097]); hi=np.array([69.90,0.0142,0.0132])
Xs=(X-lo)/(hi-lo)

def gpfit(y,seed):
    ker=ConstantKernel(1.0,(1e-3,1e3))*Matern(length_scale=np.ones(3)*0.35,
        length_scale_bounds=(0.03,5.0),nu=2.5)+WhiteKernel(1e-7,(1e-10,1e-2))
    gp=GaussianProcessRegressor(kernel=ker,normalize_y=True,n_restarts_optimizer=20,
                                random_state=seed)
    gp.fit(Xs,np.asarray(y,float)); return gp

gpP=gpfit(df.dp,2301)
gpD=gpfit(df.dd,2302)
gpS=gpfit(df.sn2,2303)

sob=qmc.Sobol(d=3,scramble=True,seed=230923)
G=sob.random_base2(m=17)
# Map Sobol cube to explicit trust box.
phys=lo+G*(hi-lo)

mP,sP=gpP.predict(G,return_std=True)
mD,sD=gpD.predict(G,return_std=True)
mS,sS=gpS.predict(G,return_std=True)

# Current best exactly feasible point is late266; constrained EI tries to improve
# Planck while preserving DESI and the >=2-SN condition.
feas=df[df.feasible]
bestP=float(feas.dp.min()) if len(feas) else 0.0
imp=bestP-mP-0.01
z=np.divide(imp,sP,out=np.zeros_like(imp),where=sP>1e-12)
ei=imp*norm.cdf(z)+sP*norm.pdf(z)
pD=norm.cdf(np.divide(-mD,sD,out=np.where(mD<0,np.inf,-np.inf),where=sD>1e-12))
pS=norm.cdf(np.divide(-mS,sS,out=np.where(mS<0,np.inf,-np.inf),where=sS>1e-12))
acq=ei*pD*pS

order=np.argsort(acq)[::-1]
chosen=[]
for j in order:
    x=G[j]
    if np.min(np.linalg.norm(Xs-x,axis=1))<0.035: continue
    if any(np.linalg.norm(x-c)<0.075 for c in chosen): continue
    chosen.append(x)
    if len(chosen)>=8: break

props=[]
for i,x in enumerate(chosen):
    j=int(np.argmin(np.linalg.norm(G-x,axis=1)))
    H0,A,B=phys[j]
    h=H0/100.
    Ox=1.-(0.02208511574370519+0.12298509428428875+4.17998772e-5)/(h*h)
    props.append(dict(
      id=f"cbo{i:03d}",H0=float(H0),A=float(A),B=float(B),Ox=float(Ox),
      pred_dp=float(mP[j]),sigma_dp=float(sP[j]),
      pred_dd=float(mD[j]),sigma_dd=float(sD[j]),
      pred_sn2=float(mS[j]),sigma_sn2=float(sS[j]),
      p_desi_feasible=float(pD[j]),p_sn2_feasible=float(pS[j]),
      constrained_ei=float(acq[j]),
      pred_strict_goal=float(max(mP[j],mD[j],mS[j]))
    ))

df.to_csv(OUT/"training_exact.csv",index=False)
pd.DataFrame(props).to_csv(OUT/"proposals.csv",index=False)
summary={
 "objective":"constrained exact BO: minimize Planck delta while requiring DESI delta<0 and second-smallest SN delta<0",
 "success_definition":"Planck<local021, DESI<local021, and at least two of Pantheon+/Union3/DES-Y5<local021",
 "best_exact_feasible_planck_delta":bestP,
 "training":raw,
 "kernels":{"planck":str(gpP.kernel_),"desi":str(gpD.kernel_),"sn2":str(gpS.kernel_)},
 "proposals":props,
 "rule":"No proposal is accepted until exact covariant stability + full Plik + raw DESI full shape + exact SN covariance all complete."
}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("EXACT_CONSTRAINED_BO_H0AB",json.dumps(summary,sort_keys=True),flush=True)
