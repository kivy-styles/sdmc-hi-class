#!/usr/bin/env python3
from pathlib import Path
import importlib.util,numpy as np,json
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)
DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
OUT=Path("output/kids_ia_components"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2,1e5,1800)
C=299792.458
C1RHO=0.0134

def tail(z,y): return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]
def source():
  with fits.open(DATA) as f:
    nz=f["NZ_SOURCE"].data
    return np.asarray(nz["Z_MID"],float),np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])

def components(m,zsrc,nzbins,variant):
  zmax=min(m["zq"].max(),float(np.max(zsrc))); z=np.linspace(.01,zmax,420)
  H=np.interp(z,m["bg"]["z"],m["bg"]["H"]); chi=np.interp(z,m["bg"]["z"],m["bg"]["chi"])
  nz=[]
  for a in nzbins:
    y=np.maximum(np.interp(z,zsrc,a,left=0,right=0),0); y/=np.trapezoid(y,z); nz.append(y)
  nz=np.asarray(nz); nchi=nz*H[None,:]
  g=np.asarray([tail(z,y)-chi*tail(z,y/np.maximum(chi,1e-12)) for y in nz])
  zm=np.broadcast_to(z[None,:],(len(ELL),len(z))); km=(ELL[:,None]+.5)/np.maximum(chi[None,:],1e-12)
  Q=b.eval_cube(m["Iq"],zm,km); Pl=b.eval_cube(m["Ip"],zm,km)
  ok=np.isfinite(Q)&np.isfinite(Pl); Q=np.where(ok,Q,0); Pl=np.where(ok,Pl,0)
  P=Pl
  if variant=="unscreened_nl":
    Pn=b.eval_cube(m["Ipnl"],zm,km); Pn=np.where(np.isfinite(Pn),Pn,Pl)
    Q=Q*np.divide(Pn,Pl,out=np.ones_like(Pn),where=Pl>0)
    P=Pn
  pref=np.exp(m["Ip"](np.column_stack([z,np.full_like(z,np.log(.05))])))
  p0=float(np.exp(m["Ip"]([[0.0,np.log(.05)]])[0]))
  D=np.maximum(np.sqrt(pref/p0),1e-4)
  F=-C1RHO*m["Om"]/D
  X=np.sqrt(np.maximum(Q*P,0))
  out={q:{"ell":ELL} for q in ["gg","gi","ii"]}
  for i in range(6):
    for j in range(i,6):
      key=f"bin_{j+1}_{i+1}"
      out["gg"][key]=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
      out["gi"][key]=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*F/np.maximum(chi,1e-12))[None,:]*X,x=chi,axis=1)
      out["ii"][key]=np.trapezoid((nchi[i]*nchi[j]*F**2/np.maximum(chi,1e-12)**2)[None,:]*P,x=chi,axis=1)
  return out

def main():
  zsrc,nzbins=source()
  hL=.6971482083084993; OmL=(.022083219194622913+.12299536722293603)/hL**2
  hC=.6856859; OmC=(.02240637+.11824151)/hC**2
  ms=[b.load_model("late300","output/weyl_spectra",hL,OmL),b.load_model("local021","output/weyl_spectra",hC,OmC)]
  for variant in ["linear","unscreened_nl"]:
    for m in ms:
      x=components(m,zsrc,nzbins,variant)
      for comp,d in x.items():
        np.savez(OUT/f"{m['name']}_{variant}_{comp}.npz",**d)
        print("KIDS_IA_COMPONENT",m["name"],variant,comp,flush=True)
if __name__=="__main__":main()
