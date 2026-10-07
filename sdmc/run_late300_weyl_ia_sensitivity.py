#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,numpy as np
from scipy.optimize import minimize

S=importlib.util.spec_from_file_location("nuis","sdmc/run_late300_weyl_nuisance.py")
n=importlib.util.module_from_spec(S); S.loader.exec_module(n)

OUT=Path("output/weyl_ia_sensitivity"); OUT.mkdir(parents=True,exist_ok=True)
full=n.b.load_des(n.DATA,n.CUTS,False)
models=[n.late,n.lcdm]

starts={
"late300":np.array([-.01996655517480449,-.028269787372995955,.00434723235407561,-.008192635365377732,-.008551960903433619,-.021298921457268423,-.024790589474497034,-.038722417486105054,-1.0441760511772582,5.0]),
"local021":np.array([-.019966555142300856,-.017034100099934335,.006210056359096395,.0024667586291670257,-.004822310884252108,-.020930848744914042,-.023228585939621387,-.03677611564560728,-1.0329210331122793,5.0])
}

def fit(model,mode):
    M=n.NativeWeylModel(model,full); p0=starts[model["name"]].copy()
    base_bounds=[(-.08,.08)]*4+[(-.1,.1)]*4+[(-5,5)]
    if mode=="fixed0":
        q0=p0[:9]
        def obj(q): return M.objective(np.r_[q,0.0])
        bounds=base_bounds
    else:
        lim=float(mode.replace("bound",""))
        p0[9]=np.clip(p0[9],-lim,lim)
        q0=p0
        def obj(q): return M.objective(q)
        bounds=base_bounds+[(-lim,lim)]
    best=None
    for dA in [0.0,1.0,-1.0]:
        q=q0.copy(); q[8]=np.clip(q[8]+dA,-5,5)
        r=minimize(obj,q,method="L-BFGS-B",bounds=bounds,options=dict(maxiter=160,ftol=1e-9,maxls=30))
        if best is None or r.fun<best.fun: best=r
    if mode=="fixed0": p=np.r_[best.x,0.0]
    else: p=best.x
    return dict(chi2=float(best.fun),success=bool(best.success),nfev=int(best.nfev),nit=int(best.nit),
                dz=p[:4].tolist(),m=p[4:8].tolist(),A1=float(p[8]),alpha1=float(p[9]))

out={}
for mode in ["fixed0","bound2","bound5"]:
    rec={}
    for model in models:
        print("IA_SENS_START",mode,model["name"],flush=True)
        rec[model["name"]]=fit(model,mode)
        print("IA_SENS_DONE",mode,model["name"],json.dumps(rec[model["name"]],sort_keys=True),flush=True)
    rec["delta_chi2"]=rec["late300"]["chi2"]-rec["local021"]["chi2"]
    out[mode]=rec
    print("IA_SENS_RESULT",mode,json.dumps(rec,sort_keys=True),flush=True)
(OUT/"ia_sensitivity.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
