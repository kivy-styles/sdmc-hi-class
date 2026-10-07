#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,numpy as np,sacc
from scipy.optimize import minimize
S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)
SP=Path("cosmosis-standard-library/likelihood/hsc_cosmic_shear/hsc_y3_fourier_shear.sacc")
s=sacc.Sacc.load_fits(str(SP)); TR=[f"wl_{i}" for i in range(4)]
pairs=[]; inds=[]; dat=[]
for i in range(4):
 for j in range(i,4):
  ell,cl,ix=s.get_ell_cl("cl_ee",TR[i],TR[j],return_ind=True)
  ell=np.asarray(ell); ix=np.asarray(ix); m=(ell>=300)&(ell<=1800); pos=np.where(m)[0]
  w=s.get_bandpower_windows(ix)
  pairs.append((i,j,pos,np.asarray(w.values,float),np.asarray(w.weight,float)[:,pos]))
  inds.extend(ix[pos]); dat.extend(np.asarray(cl)[pos])
inds=np.asarray(inds,int); dat=np.asarray(dat,float)
cov=np.asarray(s.covariance.covmat)[np.ix_(inds,inds)]; L=np.linalg.cholesky(cov)
ztr=[np.asarray(s.tracers[t].z,float) for t in TR]; nztr=[np.asarray(s.tracers[t].nz,float) for t in TR]
hL=.6971482083084993; omL=(.022083219194622913+.12299536722293603)/hL**2
hC=.6856859; omC=(.02240637+.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,omL); lcdm=b.load_model("local021","output/weyl_spectra",hC,omC)
ELL=np.geomspace(10,2e4,360); C1=.0134
def tail(z,y):
 zr=z[::-1]; yr=y[::-1]; c=np.r_[0,np.cumsum((yr[1:]+yr[:-1])*.5*np.diff(zr))]; return (-c)[::-1]
def prep(m,lam):
 z=np.linspace(.01,min(3.5,m["zq"].max(),m["zp"].max()),220); chi=np.interp(z,m["bg"]["z"],m["bg"]["chi"]); H=np.interp(z,m["bg"]["z"],m["bg"]["H"])
 zm=np.broadcast_to(z,(len(ELL),len(z))); kk=(ELL[:,None]+.5)/chi
 Q=b.eval_cube(m["Iq"],zm,kk); P=b.eval_cube(m["Ip"],zm,kk); N=b.eval_cube(m["Ipnl"],zm,kk)
 ok=np.isfinite(Q)&np.isfinite(P); Q=np.where(ok,Q,0); P=np.where(ok,P,0); N=np.where(np.isfinite(N),N,P)
 boost=np.clip(np.divide(N,P,out=np.ones_like(P),where=P>0),1,50); fac=boost**lam; Q*=fac; P*=fac
 pr=np.exp(m["Ip"](np.column_stack([z,np.full_like(z,np.log(.05))]))); p0=float(np.exp(m["Ip"]([[0,np.log(.05)]])[0]))
 F=-C1*m["Om"]/np.maximum(np.sqrt(pr/p0),1e-4)
 return z,chi,H,Q,P,np.sqrt(np.maximum(Q*P,0)),F
def src(z,chi,H):
 dz=[0,0,.07,.15]; nz=[]
 for i in range(4):
  y=np.maximum(np.interp(z-dz[i],ztr[i],nztr[i],left=0,right=0),0); y/=np.trapezoid(y,z); nz.append(y)
 nz=np.asarray(nz); g=[]
 for y in nz:g.append(tail(z,y)-chi*tail(z,y/chi))
 return nz,np.asarray(g),nz*H
def theory(m,lam,x):
 A,eta=x[0],x[1]; mc=np.asarray(x[2:6]); z,chi,H,Q,P,X,F=prep(m,lam); nz,g,nchi=src(z,chi,H); F=F*((1+z)/1.62)**eta
 cls={}
 for i in range(4):
  for j in range(i,4):
   gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
   gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*F/chi)[None,:]*X,x=chi,axis=1)
   ii=np.trapezoid((nchi[i]*nchi[j]*F**2/chi**2)[None,:]*P,x=chi,axis=1)
   cls[(i,j)]=(gg+A*gi+A*A*ii)*(1+mc[i])*(1+mc[j])
 out=[]
 for i,j,pos,wv,ww in pairs:
  cl=cls[(i,j)]; cv=np.zeros_like(wv); q=wv>0; cv[q]=np.interp(np.log(wv[q]),np.log(ELL),cl,left=cl[0],right=cl[-1]); out.extend(cv@ww)
 return np.asarray(out)
def fit(m,lam):
 def f(x):
  y=np.linalg.solve(L,dat-theory(m,lam,x)); return float(y@y+np.sum((x[2:6]/.01)**2))
 best=None
 for A0 in [0,1,-1]:
  x=np.r_[A0,0,np.zeros(4)]; r=minimize(f,x,method="L-BFGS-B",bounds=[(-6,6),(-6,6)]+[(-.1,.1)]*4,options={"maxiter":80})
  if best is None or r.fun<best.fun:best=r
 return {"chi2":float(best.fun),"AIA":float(best.x[0]),"etaIA":float(best.x[1]),"m":best.x[2:6].tolist(),"success":bool(best.success)}
out={"status":"intermediate HSC-Y3 cl_ee native-Weyl pilot; NLA+m profile, photo-z fixed to HSC fiducials, no TATT/PSF","ndata":len(dat),"lambda":{}}
for lam in [0,.55,1]:
 a=fit(late,lam); c=fit(lcdm,lam); out["lambda"][str(lam)]={"late300":a,"local021":c,"delta_chi2":a["chi2"]-c["chi2"]}
 print("HSCY3_WEYL",lam,json.dumps(out["lambda"][str(lam)],sort_keys=True),flush=True)
O=Path("output/hscy3_weyl"); O.mkdir(parents=True,exist_ok=True); (O/"summary.json").write_text(json.dumps(out,indent=2)+"\n")
