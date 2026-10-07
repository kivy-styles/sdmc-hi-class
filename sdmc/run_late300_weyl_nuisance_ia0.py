#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,numpy as np
from scipy.optimize import minimize

S=importlib.util.spec_from_file_location("n","sdmc/run_late300_weyl_nuisance.py")
n=importlib.util.module_from_spec(S); S.loader.exec_module(n)

def fit_fixed_alpha(model):
    # Optimize dz[4], m[4], A1 with alpha1 fixed to zero.
    x0=np.r_[np.zeros(4),n.M_MU,0.2]
    bounds=[(-0.08,0.08)]*4+[(-0.1,0.1)]*4+[(-5,5)]
    def obj(x):
        p=np.r_[x,0.0]
        return model.objective(p)
    best=None
    for A in [0.0,1.0,-1.0]:
        x=x0.copy(); x[8]=A
        r=minimize(obj,x,method="L-BFGS-B",bounds=bounds,options=dict(maxiter=220,ftol=1e-9,maxls=30))
        if best is None or r.fun<best.fun: best=r
    p=best.x
    return dict(chi2_total=float(best.fun),success=bool(best.success),
                dz=p[:4].tolist(),m=p[4:8].tolist(),A1=float(p[8]),alpha1=0.0,
                nfev=int(best.nfev),nit=int(best.nit),message=str(best.message))

def main():
    full=n.b.load_des(n.DATA,n.CUTS,False)
    cases={"official":n.subset(full,0,0),"c40_150":n.subset(full,40,150),"c100_250":n.subset(full,100,250)}
    out={}
    for label,des in cases.items():
        rec={"ndata":len(des["data"])}
        for m in (n.late,n.lcdm):
            print("IA0_START",label,m["name"],flush=True)
            fit=fit_fixed_alpha(n.NativeWeylModel(m,des))
            rec[m["name"]]=fit
            print("IA0_DONE",label,m["name"],json.dumps(fit,sort_keys=True),flush=True)
        rec["delta_chi2_late300_minus_local021"]=rec["late300"]["chi2_total"]-rec["local021"]["chi2_total"]
        out[label]=rec
        print("IA0_RESULT",label,json.dumps(rec,sort_keys=True),flush=True)
    R=Path("output/weyl_nuisance_ia0"); R.mkdir(parents=True,exist_ok=True)
    (R/"late300_weyl_nuisance_ia0.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
if __name__=="__main__": main()
