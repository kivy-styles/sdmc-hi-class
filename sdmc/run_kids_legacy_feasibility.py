#!/usr/bin/env python3
from pathlib import Path
import json, numpy as np
from astropy.io import fits

p=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
if not p.exists():
    raise SystemExit(f"missing {p}")

out={"file":str(p),"size_bytes":p.stat().st_size,"hdus":[]}
with fits.open(p) as f:
    for h in f:
        rec={"name":h.name,"type":type(h).__name__}
        if getattr(h,"data",None) is not None:
            try: rec["shape"]=list(h.data.shape)
            except: pass
            try: rec["columns"]=list(h.columns.names)
            except: pass
        out["hdus"].append(rec)
    # generic discovery of common two-point content
    for name in ["En","NZ_SOURCE","nz_source","COVMAT","covmat"]:
        if name in f:
            h=f[name]
            out[name]={"shape":list(h.data.shape) if getattr(h,"data",None) is not None else None}
            try: out[name]["columns"]=list(h.columns.names)
            except: pass

# auxiliary priors/covariances
for q in ["Nz_covariance.txt","massdep_cov.txt"]:
    a=np.loadtxt(Path("cosmosis-standard-library/likelihood/KiDS-Legacy")/q)
    out[q]={"shape":list(a.shape),"condition":float(np.linalg.cond(a)) if a.ndim==2 else None}

Path("output/kids_feasibility").mkdir(parents=True,exist_ok=True)
Path("output/kids_feasibility/kids_legacy_data_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("KIDS_LEGACY_DATA_AUDIT",json.dumps(out,sort_keys=True))
