#!/usr/bin/env python3
from pathlib import Path
import importlib.util, numpy as np, json
from astropy.io import fits

S=importlib.util.spec_from_file_location("base","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
OUT=Path("output/kids_native_weyl"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2,1e5,1800)

def load_kids():
    with fits.open(DATA) as f:
        nz=f["NZ_SOURCE"].data
        z=np.asarray(nz["Z_MID"],float)
        bins=np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])
        data=np.asarray(f["En"].data["VALUE"],float)
        cov=np.asarray(f["COVMAT"].data,float)
    return z,bins,data,cov

def tail(z,y):
    from scipy.integrate import cumulative_trapezoid
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def build(model,zsrc,nzbins):
    zmax=min(model["zq"].max(),float(np.max(zsrc)))
    z=np.linspace(0.01,zmax,420)
    H=np.interp(z,model["bg"]["z"],model["bg"]["H"])
    chi=np.interp(z,model["bg"]["z"],model["bg"]["chi"])
    nz=[]
    for a in nzbins:
        y=np.maximum(np.interp(z,zsrc,a,left=0,right=0),0)
        y/=np.trapezoid(y,z)
        nz.append(y)
    nz=np.asarray(nz)
    g=[]
    for y in nz:
        g.append(tail(z,y)-chi*tail(z,y/np.maximum(chi,1e-12)))
    g=np.asarray(g)
    zm=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    km=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
    Q=b.eval_cube(model["Iq"],zm,km)
    Q=np.where(np.isfinite(Q),Q,0.0)
    cls={}
    for i in range(6):
      for j in range(i,6):
        cl=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
        cls[f"bin_{j+1}_{i+1}"]=cl
    return dict(ell=ELL,**cls)

def main():
    zsrc,nzbins,data,cov=load_kids()
    hL=0.6971482083084993; OmL=(0.022083219194622913+0.12299536722293603)/hL**2
    hC=0.6856859; OmC=(0.02240637+0.11824151)/hC**2
    late=b.load_model("late300","output/weyl_spectra",hL,OmL)
    lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)
    for m in (late,lcdm):
        x=build(m,zsrc,nzbins)
        np.savez(OUT/f"{m['name']}_kids_shear_cl.npz",**x)
        print("KIDS_CL_READY",m["name"],len(x["ell"]),len(x)-1)
    meta={"ndata":int(len(data)),"cov_shape":list(cov.shape),"nbin":6,"nmode":6,
          "status":"linear native-Weyl KiDS COSEBI pipeline validation; no IA/nonlinear-MG"}
    (OUT/"metadata.json").write_text(json.dumps(meta,indent=2)+"\n")
if __name__=="__main__": main()
