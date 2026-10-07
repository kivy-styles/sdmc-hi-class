#!/usr/bin/env python3
from pathlib import Path
import importlib.util,numpy as np,json
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)
DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
OUT=Path("output/kids_native_weyl"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2,1e5,1800)
def tail(z,y): return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]
def load():
  with fits.open(DATA) as f:
    nz=f["NZ_SOURCE"].data
    return np.asarray(nz["Z_MID"],float),np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])
def calc(m,zsrc,nzbins,variant):
  zmax=min(m["zq"].max(),float(np.max(zsrc))); z=np.linspace(0.01,zmax,420)
  chi=np.interp(z,m["bg"]["z"],m["bg"]["chi"])
  nz=[]
  for a in nzbins:
    y=np.maximum(np.interp(z,zsrc,a,left=0,right=0),0); y/=np.trapezoid(y,z); nz.append(y)
  nz=np.asarray(nz); g=np.asarray([tail(z,y)-chi*tail(z,y/np.maximum(chi,1e-12)) for y in nz])
  zm=np.broadcast_to(z[None,:],(len(ELL),len(z))); km=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
  Q=b.eval_cube(m["Iq"],zm,km); P=b.eval_cube(m["Ip"],zm,km)
  ok=np.isfinite(Q)&np.isfinite(P); Q=np.where(ok,Q,0); P=np.where(ok,P,0)
  if variant=="halofit_envelope":
    Pn=b.eval_cube(m["Ipnl"],zm,km); Pn=np.where(np.isfinite(Pn),Pn,P)
    boost=np.divide(Pn,P,out=np.ones_like(Pn),where=P>0); boost=np.clip(boost,0.2,50)
    Q=Q*boost
  out={"ell":ELL}
  for i in range(6):
    for j in range(i,6):
      out[f"bin_{j+1}_{i+1}"]=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
  return out
def main():
  zsrc,nzbins=load()
  hL=.6971482083084993; OmL=(.022083219194622913+.12299536722293603)/hL**2
  hC=.6856859; OmC=(.02240637+.11824151)/hC**2
  ms=[b.load_model("late300","output/weyl_spectra",hL,OmL),b.load_model("local021","output/weyl_spectra",hC,OmC)]
  for v in ["linear","halofit_envelope"]:
    for m in ms:
      x=calc(m,zsrc,nzbins,v); np.savez(OUT/f"{m['name']}_{v}_kids_shear_cl.npz",**x)
      print("KIDS_CL_VARIANT",m["name"],v,len(x["ell"]))
if __name__=="__main__":main()
