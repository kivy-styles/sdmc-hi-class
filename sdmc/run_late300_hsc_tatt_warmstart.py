#!/usr/bin/env python3
import sys,json,importlib.util
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

S=importlib.util.spec_from_file_location("tt","sdmc/run_late300_hsc_tatt_react.py")
tt=importlib.util.module_from_spec(S); S.loader.exec_module(tt)
h=tt.h

def load_one(root, needle):
    hits=[p for p in Path(root).rglob("*.json") if needle in p.name]
    if len(hits)!=1:
        raise RuntimeError((root,needle,[str(p) for p in hits]))
    return json.loads(hits[0].read_text())

def main():
    model=sys.argv[1]
    m={"late300":h.late,"local021":h.lcdm}[model]
    nla=load_one(f"warm_inputs/nla/{model}",f"{model}_us_native")
    con=load_one(f"warm_inputs/conditional/{model}",model)
    fit=nla["fit"]
    x0=np.r_[fit["dz"],fit["m"],fit["A1"],fit["alpha1"],
             con["A2"],con["alpha2"],con["bias_ta"],fit["psf_u"]]
    obj=tt.HSCReactTATT(m,"us_native")
    val=obj.validate_nla_limit()
    if val["max_rel_to_peak"]>5e-7 or val["objective_abs_diff"]>1e-5:
        raise RuntimeError("TATT NLA-limit validation failed")
    bounds=[(-1,1)]*4+[(-.1,.1)]*4+[(-6,6),(-6,6),(-6,6),(-6,6),(0,2)]+[(-5,5)]*4

    starts=[x0.copy()]
    x1=x0.copy()
    x1[8]*=0.9
    x1[10]*=1.1
    starts.append(x1)

    results=[]
    best=None
    outdir=Path("output/hsc_tatt_warm"); outdir.mkdir(parents=True,exist_ok=True)
    for i,x in enumerate(starts):
        r=minimize(obj.objective_tatt,x,method="L-BFGS-B",bounds=bounds,
                   options={"maxiter":600,"ftol":1e-10,"gtol":1e-7,"maxls":40})
        rec={"start":i,"chi2":float(r.fun),"success":bool(r.success),
             "message":str(r.message),"nfev":int(r.nfev),"nit":int(r.nit),
             "x":r.x.tolist()}
        results.append(rec)
        Path(outdir/f"{model}_start{i}.json").write_text(json.dumps(rec,indent=2)+"\n")
        print("HSC_TATT_WARM_START_DONE",model,json.dumps(rec,sort_keys=True),flush=True)
        if best is None or r.fun<best.fun: best=r

    p=best.x
    out={"status":"warm-start full joint HSC corrected-ReACT TATT refit",
         "model":model,"variant":"us_native","ndata":len(h.DATA),
         "validation":val,
         "conditional_start":{"A2":con["A2"],"alpha2":con["alpha2"],"bias_ta":con["bias_ta"],
                              "conditional_chi2":con["tatt_conditional_chi2"]},
         "fit":{"chi2_total":float(best.fun),"success":bool(best.success),"message":str(best.message),
                "dz":p[:4].tolist(),"m":p[4:8].tolist(),
                "A1":float(p[8]),"alpha1":float(p[9]),"A2":float(p[10]),
                "alpha2":float(p[11]),"bias_ta":float(p[12]),
                "psf_u":p[13:17].tolist(),"nfev":int(best.nfev),"nit":int(best.nit)},
         "all_starts":results}
    Path(outdir/f"{model}.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("HSC_TATT_WARM_RESULT",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
