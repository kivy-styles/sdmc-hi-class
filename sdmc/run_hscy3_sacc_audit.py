#!/usr/bin/env python3
from pathlib import Path
import json, numpy as np
import sacc

P=Path("cosmosis-standard-library/likelihood/hsc_cosmic_shear/hsc_y3_fourier_shear.sacc")
OUT=Path("output/hscy3_sacc_audit"); OUT.mkdir(parents=True,exist_ok=True)
s=sacc.Sacc.load_fits(str(P))

summary={}
summary["ndata"]=int(len(s.mean))
summary["data_types"]=list(s.get_data_types())
summary["tracers"]=list(s.tracers.keys())
summary["cov_shape"]=list(s.covariance.covmat.shape) if s.has_covariance() else None
summary["cov_condition"]=float(np.linalg.cond(s.covariance.covmat)) if s.has_covariance() else None
summary["tracer_info"]={}
for name,tr in s.tracers.items():
    rec={"type":type(tr).__name__}
    if hasattr(tr,"z"):
        z=np.asarray(tr.z,float); rec["zmin"]=float(z.min()); rec["zmax"]=float(z.max()); rec["nz"]=int(len(z))
    if hasattr(tr,"nz"):
        nz=np.asarray(tr.nz,float); rec["nz_integral_trapz"]=float(np.trapezoid(nz,np.asarray(tr.z,float)))
    summary["tracer_info"][name]=rec

rows=[]
for dt in s.get_data_types():
    inds=s.indices(data_type=dt)
    for ii in inds:
        d=s.data[ii]
        rec={"index":int(ii),"data_type":dt,"tracers":list(d.tracers),"value":float(d.value)}
        tags={}
        for key in ["ell","theta"]:
            try:
                val=d.get_tag(key)
                if val is not None: tags[key]=float(val)
            except Exception: pass
        rec["tags"]=tags
        w=None
        try: w=d.get_tag("window")
        except Exception: pass
        if w is not None:
            rec["window_type"]=type(w).__name__
            for attr in ["values","weight"]:
                if hasattr(w,attr):
                    a=np.asarray(getattr(w,attr))
                    rec[f"window_{attr}_shape"]=list(a.shape)
            if hasattr(w,"values"):
                a=np.asarray(w.values,float)
                if a.size:
                    rec["window_ell_min"]=float(np.nanmin(a)); rec["window_ell_max"]=float(np.nanmax(a))
        rows.append(rec)
summary["rows_preview"]=rows[:30]
summary["ell_stats"]={}
for dt in summary["data_types"]:
    try:
        ell=np.asarray(s.get_tag("ell",data_type=dt),float)
        summary["ell_stats"][dt]={"n":int(len(ell)),"min":float(np.min(ell)),"max":float(np.max(ell)),"unique":int(len(np.unique(ell)))}
    except Exception as e:
        summary["ell_stats"][dt]={"error":repr(e)}

# Explicit HSC range counts 300..1800 for cl_ee
try:
    inds=s.indices(data_type="cl_ee")
    kept=[]
    for ii in inds:
        d=s.data[ii]
        ell=float(d.get_tag("ell"))
        if 300.0 <= ell <= 1800.0: kept.append(int(ii))
    summary["cl_ee_300_1800_count"]=len(kept)
except Exception as e:
    summary["cl_ee_300_1800_error"]=repr(e)

(OUT/"hscy3_sacc_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
print("HSCY3_SACC_AUDIT",json.dumps(summary,sort_keys=True),flush=True)
