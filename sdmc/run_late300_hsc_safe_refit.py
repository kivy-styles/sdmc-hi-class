#!/usr/bin/env python3
import sys, json, importlib.util, numpy as np
from pathlib import Path

S=importlib.util.spec_from_file_location("hsc","sdmc/run_late300_hscy3_weyl.py")
h=importlib.util.module_from_spec(S); S.loader.exec_module(h)

def configure_cut(ellmax):
    keep=np.array([i for i,p in enumerate(h.POINTS) if p["ell"]<=ellmax],dtype=int)
    raw=h.COV/float(h.HARTLAP)
    n=len(keep); r=1404.0
    fac=(r-1.0)/(r-n-2.0)
    h.DATA=h.DATA[keep]
    h.COV=raw[np.ix_(keep,keep)]*fac
    h.POINTS=[h.POINTS[i] for i in keep]
    h.HARTLAP=fac
    return n,fac

def main():
    model_name=sys.argv[1]
    ellmax=float(sys.argv[2])
    n,fac=configure_cut(ellmax)
    m={"late300":h.late,"local021":h.lcdm}[model_name]
    print("HSC_SAFE_FIT_START",model_name,ellmax,n,flush=True)
    fit=h.HSCNativeWeyl(m).fit()
    out={"status":"HSC-Y3 safe-scale refitted native-Weyl NLA-z pilot",
         "model":model_name,"ellmax":ellmax,"ndata":n,"hartlap_cov_factor":fac,
         "fit":fit,
         "missing":["full TATT A2/alpha2/bias_ta sector","validated nonlinear modified-gravity/baryonic prescription"]}
    p=Path("output/hsc_safe");p.mkdir(parents=True,exist_ok=True)
    (p/f"{model_name}_ell{int(ellmax)}.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("HSC_SAFE_FIT_DONE",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
