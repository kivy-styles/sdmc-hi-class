#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,numpy as np
from scipy.optimize import minimize
S=importlib.util.spec_from_file_location("n","sdmc/run_late300_weyl_nuisance.py")
n=importlib.util.module_from_spec(S); S.loader.exec_module(n)
full=n.b.load_des(n.DATA,n.CUTS,False)
OUT=Path("output/ia_bound2_preview"); OUT.mkdir(parents=True,exist_ok=True)
starts={
"late300":np.array([-.01996655517480449,-.028269787372995955,.00434723235407561,-.008192635365377732,-.008551960903433619,-.021298921457268423,-.024790589474497034,-.038722417486105054,-1.0441760511772582,2.0]),
"local021":np.array([-.019966555142300856,-.017034100099934335,.006210056359096395,.0024667586291670257,-.004822310884252108,-.020930848744914042,-.023228585939621387,-.03677611564560728,-1.0329210331122793,2.0])
}
def fit(m):
 M=n.NativeWeylModel(m,full); q0=starts[m["name"]]
 bounds=[(-.08,.08)]*4+[(-.1,.1)]*4+[(-5,5),(-2,2)]
 r=minimize(M.objective,q0,method="L-BFGS-B",bounds=bounds,options=dict(maxiter=140,ftol=1e-8,maxls=20))
 p=r.x
 return dict(chi2=float(r.fun),success=bool(r.success),A1=float(p[8]),alpha1=float(p[9]),dz=p[:4].tolist(),m=p[4:8].tolist(),nfev=int(r.nfev))
rec={}
for m in [n.late,n.lcdm]:
 rec[m["name"]]=fit(m); print("IA_BOUND2",m["name"],json.dumps(rec[m["name"]],sort_keys=True),flush=True)
rec["delta_chi2"]=rec["late300"]["chi2"]-rec["local021"]["chi2"]
print("IA_BOUND2_RESULT",json.dumps(rec,sort_keys=True),flush=True)
(OUT/"result.json").write_text(json.dumps(rec,indent=2,sort_keys=True)+"\n")
