#!/usr/bin/env python3
from pathlib import Path
import importlib.util,numpy as np,json
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)
DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
OUT=Path("output/kids_screening"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2,1e5,1800)
C=299792.458

def tail(z,y): return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def data():
  with fits.open(DATA) as f:
    nz=f["NZ_SOURCE"].data
    return np.asarray(nz["Z_MID"],float),np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])

def calc(m,zsrc,nzbins,mode,ks_h=None,n=2):
  zmax=min(m["zq"].max(),float(np.max(zsrc))); z=np.linspace(.01,zmax,420)
  chi=np.interp(z,m["bg"]["z"],m["bg"]["chi"])
  nz=[]
  for a in nzbins:
    y=np.maximum(np.interp(z,zsrc,a,left=0,right=0),0); y/=np.trapezoid(y,z); nz.append(y)
  nz=np.asarray(nz); g=np.asarray([tail(z,y)-chi*tail(z,y/np.maximum(chi,1e-12)) for y in nz])
  zm=np.broadcast_to(z[None,:],(len(ELL),len(z))); km=(ELL[:,None]+.5)/np.maximum(chi[None,:],1e-12)
  Q=b.eval_cube(m["Iq"],zm,km); Pl=b.eval_cube(m["Ip"],zm,km); Pn=b.eval_cube(m["Ipnl"],zm,km)
  ok=np.isfinite(Q)&np.isfinite(Pl)
  Q=np.where(ok,Q,0); Pl=np.where(ok,Pl,0); Pn=np.where(np.isfinite(Pn),Pn,Pl)
  H0=100*m["h"]/C
  Agr=(1.5*m["Om"]*H0**2*(1+z))[None,:]
  sig2=np.divide(Q,Agr**2*Pl,out=np.ones_like(Q),where=(Agr**2*Pl)>0)
  sig=np.sqrt(np.clip(sig2,0.2,2.0))
  if mode=="unscreened":
    sigeff=sig
  elif mode=="full_screen":
    sigeff=np.ones_like(sig)
  elif mode=="transition":
    x=km/(float(ks_h)*m["h"])
    sigeff=1.0+(sig-1.0)/(1.0+x**n)
  else: raise ValueError(mode)
  Qnl=Agr**2*Pn*sigeff**2
  Qnl=np.where(ok,Qnl,0)
  out={"ell":ELL}
  for i in range(6):
    for j in range(i,6):
      out[f"bin_{j+1}_{i+1}"]=np.trapezoid(g[i][None,:]*g[j][None,:]*Qnl,x=chi,axis=1)
  meta={"mode":mode,"ks_h_Mpc":ks_h,"Sigma_linear_median":float(np.median(sig[ok]))}
  return out,meta

def main():
  zsrc,nzbins=data()
  hL=.6971482083084993; OmL=(.022083219194622913+.12299536722293603)/hL**2
  hC=.6856859; OmC=(.02240637+.11824151)/hC**2
  ms=[b.load_model("late300","output/weyl_spectra",hL,OmL),b.load_model("local021","output/weyl_spectra",hC,OmC)]
  specs=[("unscreened","unscreened",None),("screen_k01","transition",.1),
         ("screen_k03","transition",.3),("screen_k10","transition",1.0),
         ("full_screen","full_screen",None)]
  meta={}
  for tag,mode,ks in specs:
    meta[tag]={}
    for m in ms:
      x,mm=calc(m,zsrc,nzbins,mode,ks); np.savez(OUT/f"{m['name']}_{tag}_kids_shear_cl.npz",**x)
      meta[tag][m["name"]]=mm
      print("KIDS_SCREEN_CL",tag,m["name"],json.dumps(mm,sort_keys=True),flush=True)
  (OUT/"screening_metadata.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
