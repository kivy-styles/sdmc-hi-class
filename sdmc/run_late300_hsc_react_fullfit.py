#!/usr/bin/env python3
import sys,json,importlib.util
from pathlib import Path
import numpy as np
from scipy.stats import truncnorm
from dynesty import NestedSampler
S=importlib.util.spec_from_file_location("rx","sdmc/run_late300_hsc_react_cutscan.py")
rx=importlib.util.module_from_spec(S);S.loader.exec_module(rx);h=rx.h
name=sys.argv[1];m={"late300":h.late,"local021":h.lcdm}[name];obj=rx.HSCReact(m,"us_native")
def tn(u,s,lo,hi):
 return float(truncnorm.ppf(np.clip(u,1e-12,1-1e-12),lo/s,hi/s,loc=0,scale=s))
def pt(u):
 p=np.empty(14);p[0]=tn(u[0],h.DZ_SIG[0],-1,1);p[1]=tn(u[1],h.DZ_SIG[1],-1,1);p[2:4]=-1+2*u[2:4]
 for i in range(4,8):p[i]=tn(u[i],h.M_SIG,-.1,.1)
 p[8:10]=-6+12*u[8:10]
 for i in range(10,14):p[i]=tn(u[i],1.,-5,5)
 return p
def ll(p):
 t=obj.theory(p)
 if t is None or np.any(~np.isfinite(t)):return -1e300
 y=np.linalg.solve(obj.L,h.DATA-t);q=float(y@y)
 return -.5*q if np.isfinite(q) else -1e300
nlive=120;dlogz=.5;maxcall=180000
sam=NestedSampler(ll,pt,14,nlive=nlive,bound="multi",sample="rwalk",bootstrap=0)
sam.run_nested(dlogz=dlogz,maxcall=maxcall,print_progress=True);r=sam.results
ib=int(np.argmax(r.logl));nc=int(np.sum(r.ncall))
out={"status":"conditional fixed-cosmology HSC us_native nuisance evidence","model":name,"variant":"us_native",
"ndata":int(len(h.DATA)),"ndim_nuisance":14,"nlive":nlive,"requested_dlogz":dlogz,"maxcall":maxcall,"ncall":nc,
"terminated_by_maxcall":bool(nc>=maxcall),"logZ_shared_data_normalization_omitted":float(r.logz[-1]),
"logZerr":float(r.logzerr[-1]),"best_sample_data_chi2":float(-2*r.logl[ib]),
"qualification":"Fixed cosmology and corrected-ReACT us_native spectra; 14 HSC NLA-z nuisance coordinates marginalized under the same bounded Gaussian/top-hat priors as the completed fit. Shared data normalization cancels."}
Path("output/hsc_conditional_evidence").mkdir(parents=True,exist_ok=True)
Path(f"output/hsc_conditional_evidence/{name}.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("HSC_CONDITIONAL_EVIDENCE",json.dumps(out,sort_keys=True),flush=True)
