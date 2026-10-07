#!/usr/bin/env python3
from pathlib import Path
import json,numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits

DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
WDIR=Path("cosmosis-standard-library/shear/cosebis/WnLog")
NPZDIR=Path("output/kids_native_weyl")
OUT=Path("output/kids_fast_validate"); OUT.mkdir(parents=True,exist_ok=True)

with fits.open(DATA) as f:
    d=np.asarray(f["En"].data["VALUE"],float)
    rows=f["En"].data
    cov=np.asarray(f["COVMAT"].data,float)
L=np.linalg.cholesky(cov)

def windows(ell):
    out={}
    x=np.log(ell)
    for n in range(1,7):
        a=np.loadtxt(WDIR/f"WnLog{n}-2.00-300.00.table")
        # Stored x coordinate is log(ell) because logTable=1.
        w=np.interp(x,a[:,0],a[:,1],left=0.0,right=0.0)
        out[n]=w
    return out

def theory(npz):
    a=np.load(npz)
    ell=np.asarray(a["ell"],float)
    W=windows(ell)
    En={}
    for i in range(1,7):
        for j in range(i,7):
            cl=np.asarray(a[f"bin_{j}_{i}"],float)
            En[(i,j)]={}
            for n in range(1,7):
                En[(i,j)][n]=np.trapezoid(ell*W[n]*cl,x=ell)/(2*np.pi)
    t=[]
    for r in rows:
        i=int(r["BIN1"]); j=int(r["BIN2"]); n=int(r["ANGBIN"])
        t.append(En[(min(i,j),max(i,j))][n])
    return np.asarray(t)

out={}
official={"late300":392.7860933581886,"local021":404.42821237195534}
for m in ["late300","local021"]:
    t=theory(NPZDIR/f"{m}_kids_shear_cl.npz")
    y=np.linalg.solve(L,d-t); chi=float(y@y)
    out[m]={"chi2_fast":chi,"chi2_official_cpp":official[m],
            "difference":chi-official[m],"max_abs_theory":float(np.max(np.abs(t)))}
out["delta_chi2_fast"]=out["late300"]["chi2_fast"]-out["local021"]["chi2_fast"]
out["delta_chi2_official"]=official["late300"]-official["local021"]
print("KIDS_FAST_COSEBI_VALIDATE",json.dumps(out,sort_keys=True))
(OUT/"kids_fast_cosebi_validate.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
