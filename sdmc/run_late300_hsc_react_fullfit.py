#!/usr/bin/env python3
import sys,json,math,importlib.util
from pathlib import Path
import numpy as np
from scipy.special import ndtr
S=importlib.util.spec_from_file_location("rx","sdmc/run_late300_hsc_react_cutscan.py")
rx=importlib.util.module_from_spec(S);S.loader.exec_module(rx);h=rx.h
name=sys.argv[1];m={"late300":h.late,"local021":h.lcdm}[name];obj=rx.HSCReact(m,"us_native")
root=Path("output/hsc_map")/name
hits=list(root.rglob(f"{name}_us_native.json"))
if len(hits)!=1: raise RuntimeError(("MAP artifact",hits))
fit=json.loads(hits[0].read_text())["fit"]
p=np.r_[fit["dz"],fit["m"],fit["A1"],fit["alpha1"],fit["psf_u"]].astype(float)
scale=np.r_[.024,.022,.10,.10,np.full(4,.01),1.,1.,np.ones(4)]
step=.05*scale
def nlp(x): return .5*float(obj.objective(np.asarray(x,float)))
f0=nlp(p); n=len(p); H=np.empty((n,n),float)
for i in range(n):
 ei=np.zeros(n);ei[i]=step[i]
 H[i,i]=(nlp(p+ei)-2*f0+nlp(p-ei))/(step[i]**2)
for i in range(n):
 for j in range(i+1,n):
  ei=np.zeros(n);ej=np.zeros(n);ei[i]=step[i];ej[j]=step[j]
  v=(nlp(p+ei+ej)-nlp(p+ei-ej)-nlp(p-ei+ej)+nlp(p-ei-ej))/(4*step[i]*step[j])
  H[i,j]=H[j,i]=v
eig=np.linalg.eigvalsh(H)
sgn,logdet=np.linalg.slogdet(H)
prior=float((p[0]/.024)**2+(p[1]/.022)**2+np.sum((p[4:8]/.01)**2)+np.sum(p[10:14]**2))
data=float(2*f0-prior)
def log_tnorm_const(sig,lo,hi):
 return math.log(sig)+.5*math.log(2*math.pi)+math.log(ndtr(hi/sig)-ndtr(lo/sig))
norm=sum([log_tnorm_const(.024,-1,1),log_tnorm_const(.022,-1,1)])
norm+=4*log_tnorm_const(.01,-.1,.1)+4*log_tnorm_const(1.,-5,5)
norm+=2*math.log(2.)+2*math.log(12.)
logZ=None
if sgn>0 and eig[0]>0:
 logZ=-f0-norm+.5*n*math.log(2*math.pi)-.5*logdet
out={"status":"conditional HSC us_native Laplace nuisance evidence","model":name,"ndim":n,
"map_total_chi2":2*f0,"map_data_chi2":data,"map_prior_chi2":prior,
"hessian_step_fraction":.05,"hessian_min_eigenvalue":float(eig[0]),"hessian_max_eigenvalue":float(eig[-1]),
"hessian_positive_definite":bool(sgn>0 and eig[0]>0),"logdet_hessian":float(logdet),"logZ_laplace":None if logZ is None else float(logZ),
"log_prior_normalization_cost":float(norm),
"qualification":"Laplace approximation around the completed corrected-ReACT us_native MAP using the exact 14D HSC objective and normalized bounded priors; shared data normalization omitted."}
Path("output/hsc_laplace").mkdir(parents=True,exist_ok=True)
Path(f"output/hsc_laplace/{name}.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("HSC_LAPLACE_EVIDENCE",json.dumps(out,sort_keys=True),flush=True)
