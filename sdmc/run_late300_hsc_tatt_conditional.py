#!/usr/bin/env python3
import sys,json,importlib.util
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

S=importlib.util.spec_from_file_location("tt","sdmc/run_late300_hsc_tatt_react.py")
tt=importlib.util.module_from_spec(S); S.loader.exec_module(tt)
h=tt.h

def load_best(model):
    root=Path("output/nla_best")/model
    hits=list(root.rglob(f"{model}_us_native.json"))
    if len(hits)!=1: raise RuntimeError((model,hits))
    d=json.loads(hits[0].read_text())["fit"]
    return d

def main():
    model=sys.argv[1]
    obj=tt.HSCReactTATT({"late300":h.late,"local021":h.lcdm}[model],"us_native")
    d=load_best(model)
    fixed=np.r_[d["dz"],d["m"],d["A1"],d["alpha1"],d["psf_u"]]
    # theory_tatt vector = dz4,m4,A1,alpha1,A2,alpha2,bias_ta,psf4
    def pack(x):
        return np.r_[fixed[:10],x[0],x[1],x[2],fixed[10:14]]
    def fun(x):
        return float(obj.objective_tatt(pack(x)))
    starts=[
      np.array([0.,0.,0.]),
      np.array([-.9,-4.,.05]),
      np.array([1.,0.,.5]),
      np.array([-1.,0.,1.]),
    ]
    bounds=[(-6,6),(-6,6),(0,2)]
    best=None
    for i,x in enumerate(starts):
        r=minimize(fun,x,method="L-BFGS-B",bounds=bounds,
                   options={"maxiter":240,"ftol":2e-10,"maxls":35})
        print("HSC_TATT_CONDITIONAL_START",model,i,float(r.fun),r.x.tolist(),bool(r.success),flush=True)
        if best is None or r.fun<best.fun: best=r
    nla=float(d["chi2_total"])
    out={
      "status":"conditional TATT extension at fixed completed NLA-z nuisance optimum",
      "model":model,"variant":"us_native","ndata":len(h.DATA),
      "nla_chi2":nla,
      "tatt_conditional_chi2":float(best.fun),
      "delta_chi2_tatt_minus_nla":float(best.fun-nla),
      "A2":float(best.x[0]),"alpha2":float(best.x[1]),"bias_ta":float(best.x[2]),
      "success":bool(best.success),
      "qualification":"Only A2, alpha2 and bias_ta are reoptimized. Photo-z, shear calibration, A1, alpha1 and PSF coordinates remain at each model's completed NLA-z optimum. This is a fast conditional diagnostic while the full joint TATT fit runs."
    }
    p=Path("output/hsc_tatt_conditional"); p.mkdir(parents=True,exist_ok=True)
    (p/f"{model}.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("HSC_TATT_CONDITIONAL_DONE",json.dumps(out,sort_keys=True),flush=True)
if __name__=="__main__": main()
