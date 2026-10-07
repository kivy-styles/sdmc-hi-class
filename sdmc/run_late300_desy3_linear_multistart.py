#!/usr/bin/env python3
from pathlib import Path
import importlib.util, json, numpy as np
from scipy.optimize import minimize

S=importlib.util.spec_from_file_location("rn","sdmc/run_late300_desy3_react_nuisance.py")
rn=importlib.util.module_from_spec(S); S.loader.exec_module(rn)

# Two independently discovered official-cut minima from the earlier and ReACT-control runs.
KNOWN={
 "late300":[
  [-0.01996655517480449,-0.028269787372995955,0.00434723235407561,-0.008192635365377732,
   -0.008551960903433619,-0.021298921457268423,-0.024790589474497034,-0.038722417486105054,
   -1.0441760511772582,5.0],
  [-0.00998327736904212,-0.029242891733343948,0.004442062235434586,-0.0006975792533700345,
   -0.008851047587476835,-0.017153060206050848,-0.02835525327533222,-0.03754269234264069,
   -0.9804610389484189,5.0]
 ],
 "local021":[
  [-0.019966555142300856,-0.017034100099934335,0.006210056359096395,0.0024667586291670257,
   -0.004822310884252108,-0.020930848744914042,-0.023228585939621387,-0.03677611564560728,
   -1.0329210331122793,5.0],
  [-0.026134373258060492,-0.02994983277591963,0.0031570459428879394,0.000714718043450475,
   -0.005709400928765867,-0.020990728634128902,-0.023194463634909603,-0.04155436307281693,
   -1.1765475905161387,5.0]
 ]
}

def fit_model(m,des):
    obj=rn.ReactWeylModel(m,des,"linear")
    base=np.r_[np.zeros(4),rn.M_MU,0.2,0.0]
    starts=[base]
    starts += [np.asarray(x,float) for x in KNOWN[m["name"]]]
    # Add nearby IA alternatives to test the boundary structure.
    for A,alpha in [(-2,5),(0,5),(1,5),(-1,3),(-1,0),(0,0),(1,0)]:
        x=np.asarray(KNOWN[m["name"]][0],float).copy()
        x[8]=A; x[9]=alpha; starts.append(x)
    bounds=[(-0.08,0.08)]*4+[(-0.1,0.1)]*4+[(-5,5),(-5,5)]
    trials=[]
    best=None
    for i,x in enumerate(starts):
        before=float(obj.objective(x))
        r=minimize(obj.objective,x,method="L-BFGS-B",bounds=bounds,
                   options=dict(maxiter=500,ftol=1e-11,maxls=60,maxfun=12000))
        rec=dict(start=i,chi2_start=before,chi2=float(r.fun),success=bool(r.success),
                 message=str(r.message),nfev=int(r.nfev),nit=int(r.nit),p=r.x.tolist())
        trials.append(rec)
        print("DESY3_LINEAR_ROBUST_TRIAL",m["name"],json.dumps(rec,sort_keys=True),flush=True)
        if best is None or r.fun<best.fun: best=r
    return dict(chi2=float(best.fun),p=best.x.tolist(),success=bool(best.success),
                nfev=int(best.nfev),nit=int(best.nit),trials=trials)

def main():
    full=rn.b.load_des(rn.DATA,rn.CUTS,False)
    des=rn.subset(full,0,0)
    out={"status":"DES-Y3 official-cut linear nuisance multistart robustness",
         "ndata":len(des["data"]),"models":{}}
    for m in (rn.late,rn.lcdm):
        out["models"][m["name"]]=fit_model(m,des)
    out["delta_chi2_late300_minus_local021"]=out["models"]["late300"]["chi2"]-out["models"]["local021"]["chi2"]
    p=Path("output/desy3_linear_robust");p.mkdir(parents=True,exist_ok=True)
    (p/"summary.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("DESY3_LINEAR_ROBUST_SUMMARY",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
