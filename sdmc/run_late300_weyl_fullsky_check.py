#!/usr/bin/env python3
from pathlib import Path
import importlib.util,sys,json,numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from astropy.io import fits

S=importlib.util.spec_from_file_location("base","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)
sys.path.insert(0,"cosmosis-standard-library/shear/cl_to_xi_fullsky")
import legendre

DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
OUT=Path("output/weyl_fullsky_check"); OUT.mkdir(parents=True,exist_ok=True)
DZSIG=np.array([.018,.015,.011,.017]); MMU=np.array([-.0063,-.0198,-.0241,-.0369]); MSIG=np.array([.0091,.0078,.0076,.0076])
ZPIV=.62; C1=b.C1RHO
LATE_P=np.array([-.0016437758106129172,.002763865120787014,-.001113726521483655,.0003762068779416394,-.006399741534384516,-.019432625662701288,-.024724494015670357,-.036738989606425715,.17135058117372787,-.0010573037890486343])
LCDM_P=np.array([-.0016902354365068923,.0028367753655819472,-.001100532953537646,.000373416518216383,-.0064030639909885365,-.019422936812064837,-.02472174071973695,-.036740145007850084,.15929002527407424,-.0005098877193964763])

hL=.6971482083084993; OmL=(.022083219194622913+.12299536722293603)/hL**2
hC=.6856859; OmC=(.02240637+.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)
full=b.load_des(DATA,CUTS,False)

def sub(des):
 use=[i for i,r in enumerate(des["rows"]) if r["ang"] >= (100. if r["typ"]=="xip" else 250.)]
 ix=np.array(use,int)
 return dict(rows=[des["rows"][i] for i in use],data=des["data"][ix],cov=des["cov"][np.ix_(ix,ix)],zsrc=des["zsrc"],nz=des["nz"])
des=sub(full)

def edges():
 with fits.open(DATA) as f:
  d=f["xip"].data; names=[x.upper() for x in d.names]
  q=(d["BIN1"]==d["BIN1"][0])&(d["BIN2"]==d["BIN2"][0])
  dd=d[q]
  def col(cands):
   for c in cands:
    if c in names: return d.names[names.index(c)]
   return None
  cmin=col(["ANGMIN","ANGLEMIN","THETA_MIN","ANG_MIN"])
  cmax=col(["ANGMAX","ANGLEMAX","THETA_MAX","ANG_MAX"])
  if cmin and cmax:
   lo=np.asarray(dd[cmin],float); hi=np.asarray(dd[cmax],float)
   o=np.argsort(lo); return np.r_[lo[o][0],hi[o]]
  ang=np.sort(np.unique(np.asarray(dd["ANG"],float)))
  le=np.log(ang); mid=(le[:-1]+le[1:])/2
  e=np.empty(len(ang)+1); e[1:-1]=np.exp(mid); e[0]=np.exp(2*le[0]-mid[0]); e[-1]=np.exp(2*le[-1]-mid[-1])
  return e
theta_edges_arcmin=edges(); theta_edges=np.deg2rad(theta_edges_arcmin/60.)

def tail(z,y): return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def theory(m,p,ellmax):
 zmax=min(m["zq"].max(),m["zp"].max(),float(np.max(des["zsrc"])))
 z=np.linspace(.01,zmax,300); H=np.interp(z,m["bg"]["z"],m["bg"]["H"]); chi=np.interp(z,m["bg"]["z"],m["bg"]["chi"])
 dz=p[:4]; mm=p[4:8]; A=p[8]; alpha=p[9]
 nz=[]
 for i in range(4):
  f=interp1d(des["zsrc"],des["nz"][i],bounds_error=False,fill_value=0.0)
  y=np.maximum(f(z-dz[i]),0); y/=np.trapezoid(y,z); nz.append(y)
 nz=np.asarray(nz); g=[]
 for y in nz:
  g.append(tail(z,y)-chi*tail(z,y/np.maximum(chi,1e-12)))
 g=np.asarray(g); nchi=nz*H[None,:]
 pref=np.exp(m["Ip"](np.column_stack([z,np.full_like(z,np.log(.05))]))); p0=float(np.exp(m["Ip"]([[0.,np.log(.05)]])[0])); D=np.maximum(np.sqrt(pref/p0),1e-4)
 F=-A*C1*m["Om"]/D*((1+z)/(1+ZPIV))**alpha
 ells=np.arange(2,ellmax+1); zm=np.broadcast_to(z[None,:],(len(ells),len(z))); km=(ells[:,None]+.5)/np.maximum(chi[None,:],1e-12)
 Q=b.eval_cube(m["Iq"],zm,km); P=b.eval_cube(m["Ip"],zm,km); ok=np.isfinite(Q)&np.isfinite(P); Q=np.where(ok,Q,0); P=np.where(ok,P,0); X=np.sqrt(np.maximum(Q*P,0))
 cls={}
 for i in range(4):
  for j in range(i,4):
   gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
   gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*F/np.maximum(chi,1e-12))[None,:]*X,x=chi,axis=1)
   ii=np.trapezoid((nchi[i]*nchi[j]*F**2/np.maximum(chi,1e-12)**2)[None,:]*P,x=chi,axis=1)
   cl=np.zeros(ellmax+1); cl[2:]=(gg+gi+ii)*(1+mm[i])*(1+mm[j]); cls[(i+1,j+1)]=cl
 lp,lm=legendre.get_legfactors_22_binav(np.arange(ellmax+1),theta_edges)
 lp=legendre.apply_filter(ellmax,.75,lp); lm=legendre.apply_filter(ellmax,.75,lm)
 # map each data angle to the common angular bin
 centers=(2/3)*(theta_edges[1:]**3-theta_edges[:-1]**3)/(theta_edges[1:]**2-theta_edges[:-1]**2)
 tv=[]
 for r in des["rows"]:
  a=np.deg2rad(r["ang"]/60.); ib=int(np.argmin(np.abs(centers-a)))
  fac=lp[ib] if r["typ"]=="xip" else lm[ib]
  tv.append(np.sum(fac*cls[(min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))]))
 return np.asarray(tv)

def score(m,p,ellmax):
 t=theory(m,p,ellmax); L=np.linalg.cholesky(des["cov"]); y=np.linalg.solve(L,des["data"]-t)
 datachi=float(y@y); prior=float(np.sum((p[:4]/DZSIG)**2)+np.sum((p[4:8]-MMU)**2/MSIG**2))
 return datachi+prior,datachi,prior

out={"ndata":len(des["data"]),"theta_edges_arcmin":theta_edges_arcmin.tolist(),"results":{}}
for ellmax in [3000,6000,10000]:
 a=score(late,LATE_P,ellmax); c=score(lcdm,LCDM_P,ellmax)
 out["results"][str(ellmax)]={"late300":{"chi2_total":a[0],"chi2_data":a[1],"chi2_prior":a[2]},"local021":{"chi2_total":c[0],"chi2_data":c[1],"chi2_prior":c[2]},"delta_chi2":a[0]-c[0]}
 print("FULLSKY_CHECK",ellmax,json.dumps(out["results"][str(ellmax)],sort_keys=True),flush=True)
(OUT/"fullsky_check.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
