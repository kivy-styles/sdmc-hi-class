#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from astropy.io import fits

S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
OUT=Path("output/kids_native_weyl"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2.0,2.0e4,800)

hL=0.6971482083084993; OmL=(0.022083219194622913+0.12299536722293603)/hL**2
hC=0.6856859; OmC=(0.02240637+0.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)

with fits.open(DATA) as f:
    nz=f["NZ_SOURCE"].data
    zsrc=np.asarray(nz["Z_MID"],float)
    bins=np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def make(model):
    zmax=min(model["zq"].max(),float(zsrc.max()))
    z=np.linspace(max(0.01,float(zsrc[zsrc>0].min())),zmax,420)
    chi=np.interp(z,model["bg"]["z"],model["bg"]["chi"])
    nzs=[]
    for arr in bins:
        y=np.maximum(interp1d(zsrc,arr,bounds_error=False,fill_value=0.0)(z),0.0)
        y/=np.trapezoid(y,z)
        nzs.append(y)
    nzs=np.asarray(nzs)
    g=[]
    for y in nzs:
        q0=tail(z,y); q1=tail(z,y/np.maximum(chi,1e-12))
        g.append(q0-chi*q1)
    g=np.asarray(g)
    zm=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    km=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
    Q=b.eval_cube(model["Iq"],zm,km)
    valid=np.isfinite(Q)
    Q=np.where(valid,Q,0.0)
    out={"ell":ELL}
    for i in range(6):
        for j in range(i,6):
            cl=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
            out[f"bin_{j+1}_{i+1}"]=cl
            out[f"bin_{i+1}_{j+1}"]=cl
    meta=dict(model=model["name"],zmin=float(z.min()),zmax=float(z.max()),
              valid_fraction=float(valid.mean()),nbin=6,ell_min=float(ELL.min()),ell_max=float(ELL.max()),
              note="Native linear Weyl GG only; fiducial KiDS n(z); no IA/photo-z nuisance in baseline.")
    return out,meta

summary={}
for model in (late,lcdm):
    arrays,meta=make(model)
    np.savez(OUT/f"{model['name']}_kids_shear_cl.npz",**arrays)
    summary[model["name"]]=meta
    print("KIDS_NATIVE_WEYL_CL",json.dumps(meta,sort_keys=True),flush=True)
(OUT/"kids_native_weyl_cl_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
