#!/usr/bin/env python3
"""
Expanded exact-informed hard-gate Bayesian controller.

Coordinates:
 H0, late-lobe A, late-lobe B, structural interpolation s
where s=0 is the exactbo000/seed-like kinetic point and s=1 is sobol006.

Hard success:
 Planck-local021 < 0
 DESI-local021 < 0
 second-smallest of 3 exact-SN deltas < 0
"""
from pathlib import Path
import json, numpy as np, pandas as pd
from scipy.stats import norm,qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel,Matern,WhiteKernel
from sklearn.preprocessing import StandardScaler

OUT=Path("output/hardgate_structural_bo_round2"); OUT.mkdir(parents=True,exist_ok=True)
P0=2781.4370287015954; D0=331.93726232610265
SN0={"pp":1406.214941392187,"u3":28.78345144550258,"d5":1650.6363433003426}

obs=[
{"id":"exactbo000","H0":69.81677547364065,"A":.010977504386967048,"B":.012857573034386149,"s":0.0,
 "P":2784.05650677952,"D":329.34250909183334,"pp":1406.8001879565418,"u3":29.341775108259753,"d5":1651.9126336723566},
{"id":"exactbo004","H0":69.82911676766862,"A":.010557210430167612,"B":.012939807320889084,"s":.04123711327583436,
 "P":2784.0357790433236,"D":329.30971730204203,"pp":1406.8435462792404,"u3":29.377469697079505,"d5":1651.997425943613},
{"id":"exactbo2000","H0":69.80878993044607,"A":.011249459300190211,"B":.01280436261370778,"s":.22680412371133876,
 "P":2784.0641657880965,"D":329.36357787562383,"pp":1406.7722473018803,"u3":29.318722501484444,"d5":1651.857893191278},
{"id":"late298","H0":69.58899342644959,"A":.013419163235463202,"B":.011442640743218362,"s":0.0,
 "P":2784.4077241594355,"D":330.1382662686478,"pp":1406.21248516161,"u3":28.839521281566704,"d5":1650.7453144714236},
{"id":"late266","H0":69.45160505326464,"A":.013036667098291217,"B":.010544837576337158,"s":0.0,
 "P":2784.6643534348414,"D":330.69550645610434,"pp":1405.935751490295,"u3":28.592742397959228,"d5":1650.182109862566},
{"id":"hybrid043","H0":69.486625116691,"A":.007934502487070859,"B":.01100946151232347,"s":0.0,
 "P":2784.92363121286,"D":332.2669697694018,"pp":1406.1856789663434,"u3":28.80288930531242,"d5":1650.692088894546},
{"id":"hybrid063","H0":69.44595765553416,"A":.010335846113041045,"B":.010714279787614942,"s":0.0,
 "P":2784.790105139299,"D":332.68654218492713,"pp":1406.0198865132406,"u3":28.662211059985566,"d5":1650.3558155372739},
{"id":"hybrid167","H0":69.51717766530813,"A":.008974749634042383,"B":.011274873235728592,"s":0.0,
 "P":2784.9246844210033,"D":331.891151096234,"pp":1406.2192904320545,"u3":28.83534829504788,"d5":1650.7610308378935},
{"id":"exactsn_e0231","H0":69.35292106793452,"A":.011241707768261,"B":.0102705562927805,"s":1.0,
 "P":2785.1734304764773,"D":331.1506431026279,"pp":1405.8037591041066,"u3":28.47129512768879,"d5":1649.913106791675}
]
for r in obs:
 r["dP"]=r["P"]-P0; r["dD"]=r["D"]-D0
 r["dPP"]=r["pp"]-SN0["pp"];r["dU3"]=r["u3"]-SN0["u3"];r["dD5"]=r["d5"]-SN0["d5"]
 r["dSN2"]=sorted([r["dPP"],r["dU3"],r["dD5"]])[1]
 r["goal"]=max(r["dP"],r["dD"],r["dSN2"])
 r["closed"]=bool(r["dP"]<0 and r["dD"]<0 and r["dSN2"]<0)

df=pd.DataFrame(obs)
features=["H0","A","B","s"]; X=df[features].to_numpy(float)
sc=StandardScaler().fit(X); Xs=sc.transform(X)

def fit(y,seed):
 ker=ConstantKernel(1.,(1e-3,1e3))*Matern(np.ones(4),(0.08,10.),nu=2.5)+WhiteKernel(1e-7,(1e-10,1e-2))
 gp=GaussianProcessRegressor(kernel=ker,normalize_y=True,n_restarts_optimizer=25,random_state=seed)
 gp.fit(Xs,np.asarray(y,float)); return gp
gpP=fit(df.dP,2201); gpD=fit(df.dD,2202); gpS=fit(df.dSN2,2203)

# Wider but still local/extrapolative region.
lo=np.array([69.20,.0070,.0096,-.20])
hi=np.array([69.92,.0144,.0135,1.25])
sob=qmc.Sobol(d=4,scramble=True,seed=23092302)
G=lo+(hi-lo)*sob.random_base2(m=16)
Gs=sc.transform(G)
mP,sP=gpP.predict(Gs,return_std=True); mD,sD=gpD.predict(Gs,return_std=True); mS,sS=gpS.predict(Gs,return_std=True)
pP=norm.cdf(-mP/np.maximum(sP,1e-9));pD=norm.cdf(-mD/np.maximum(sD,1e-9));pS=norm.cdf(-mS/np.maximum(sS,1e-9))
pfeas=pP*pD*pS
mG=np.maximum.reduce([mP,mD,mS]); sG=np.maximum.reduce([sP,sD,sS]); best=float(df.goal.min())
imp=best-mG-.01; z=imp/np.maximum(sG,1e-9); ei=imp*norm.cdf(z)+sG*norm.pdf(z)
# Stronger feasibility weight than round1, retain uncertainty exploration.
acq=ei*(.05+.95*pfeas)+.20*pfeas+.015*(sP+sD+sS)
order=np.argsort(acq)[::-1]; chosen=[]
span=hi-lo
for j in order:
 x=G[j]
 if np.min(np.linalg.norm((X-x)/span,axis=1))<.035: continue
 if any(np.linalg.norm((x-c)/span)<.07 for c in chosen): continue
 chosen.append(x)
 if len(chosen)==8: break

# Structural endpoint interpolation.
K0=np.array([.023604633340554453,3.4328776987879466,.36879721635160295,.050275813910968366])
K1=np.array([.024290704212870225,3.6096407580714724,.3631213944496205,.04602317962943721])
OB0=.02208511574370519; OB1=.022083219194622913
OC0=.12298509428428875; OC1=.12299536722293603
AS0=2.1197359302086e-9; AS1=2.119721074504234e-9
NS0=.962555355281866; NS1=.9625227132590487
props=[]
for i,x in enumerate(chosen):
 j=int(np.argmin(np.sum((G-x)**2,axis=1))); H,A,B,s=map(float,x); K=K0+s*(K1-K0)
 props.append({"id":f"structbo{i:03d}","H0":H,"A":A,"B":B,"s":s,
  "AF":float(K[0]),"ZC":float(K[1]),"width":float(K[2]),"D_floor":float(K[3]),
  "omega_b":float(OB0+s*(OB1-OB0)),"omega_cdm":float(OC0+s*(OC1-OC0)),
  "A_s":float(AS0+s*(AS1-AS0)),"n_s":float(NS0+s*(NS1-NS0)),"tau":.055202901571989066,
  "gp_dP":float(mP[j]),"gp_dD":float(mD[j]),"gp_dSN2":float(mS[j]),
  "pP":float(pP[j]),"pD":float(pD[j]),"pSN2":float(pS[j]),"p_all":float(pfeas[j]),"acq":float(acq[j])})
df.to_csv(OUT/"training_exact.csv",index=False);pd.DataFrame(props).to_csv(OUT/"proposals.csv",index=False)
summary={"best_exact_goal":best,"training":obs,"kernels":{"P":str(gpP.kernel_),"D":str(gpD.kernel_),"SN2":str(gpS.kernel_)},
 "proposals":props,"rule":"Exact full-Plik < local021 AND raw DESI < local021 AND at least 2/3 exact SN < local021. No surrogate-only promotion."}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("HARDGATE_STRUCTURAL_BO",json.dumps(summary,sort_keys=True),flush=True)
