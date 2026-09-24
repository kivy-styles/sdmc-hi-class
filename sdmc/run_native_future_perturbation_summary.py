#!/usr/bin/env python3
"""
Summarize isolated native future hi_class scalar perturbation propagation.

The workflow writes perturbation histories at three representative k values
for both mature structural endpoints.  This script checks that the native
histories actually extend beyond a=1, remain finite, and reports the future
amplitude behavior of the Horndeski scalar, metric potentials and matter.

This is deliberately diagnostic: it does not impose a stability threshold on
gauge-dependent raw amplitudes.  The hard gates are successful CLASS return
code, a_max>1, finite data, and the absence of numerical blow-up/non-finite
entries.  Gauge-invariant interpretation remains tied to the action-level
Q_s/c_s^2 audit and the curvature-envelope audit.
"""
from pathlib import Path
import json,re,math
import numpy as np

OUT=Path("output/native_future_perturbation_summary.json")
RUNS={
  "matched_chi_F":{
    "prefix":"native_future_perturb_matched_",
    "rc":Path("output/native_future_perturb_matched.rc"),
    "log":Path("output/native_future_perturb_matched.log"),
  },
  "bounded_chi":{
    "prefix":"native_future_perturb_bounded_",
    "rc":Path("output/native_future_perturb_bounded.rc"),
    "log":Path("output/native_future_perturb_bounded.log"),
  },
}

def read_titles(path):
    lines=path.read_text(errors="replace").splitlines()
    k=None
    for line in lines[:5]:
        m=re.search(r"mode k\s*=\s*([0-9eE+\-.]+)",line)
        if m: k=float(m.group(1))
    hdr=None
    for line in lines:
        if line.startswith("#") and re.search(r"(?:^|\s)1\s*:",line):
            hdr=line.lstrip("#").strip()
    if hdr is None:
        raise RuntimeError(f"no numbered header in {path}")
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    arr=np.loadtxt(path)
    if arr.ndim==1: arr=arr[None,:]
    if arr.shape[1]!=len(names):
        raise RuntimeError(f"{path}: {arr.shape[1]} columns vs {len(names)} titles")
    return k,names,arr

def col(names,arr,*candidates):
    for c in candidates:
        if c in names:
            return arr[:,names.index(c)],c
    return None,None

def amp_summary(a,y):
    fut=a>=1.
    if not np.any(fut):
        return {"status":"no_future_rows"}
    aa=a[fut]; yy=np.asarray(y)[fut]
    j0=0
    scale=max(abs(float(yy[j0])),1e-300)
    finite=np.isfinite(yy)
    return {
      "finite":bool(np.all(finite)),
      "present_or_first_future":float(yy[j0]),
      "max_abs_future":float(np.nanmax(np.abs(yy))),
      "max_abs_over_first_future":float(np.nanmax(np.abs(yy))/scale),
      "final":float(yy[-1]),
      "final_abs_over_first_future":float(abs(yy[-1])/scale),
      "a_at_max_abs":float(aa[int(np.nanargmax(np.abs(yy)))]),
    }

out={"status":"native future perturbation summary","runs":{}}
for name,cfg in RUNS.items():
    rc=None
    if cfg["rc"].exists():
        try: rc=int(cfg["rc"].read_text().strip())
        except Exception: pass
    log_tail=cfg["log"].read_text(errors="replace").splitlines()[-30:] if cfg["log"].exists() else []
    files=sorted(Path("output").glob(cfg["prefix"]+"*perturbations_k*_s.dat"))
    run={"returncode":rc,"log_tail":log_tail,"files":[]}
    if not files:
        run["status"]="missing_perturbation_files"
        out["runs"][name]=run
        continue
    all_ok=(rc==0)
    for path in files:
        try:
            k,names,arr=read_titles(path)
            a,an=col(names,arr,"a")
            if a is None: raise RuntimeError("missing a column")
            future=a>=1.
            finite_all=bool(np.all(np.isfinite(arr)))
            amax=float(np.nanmax(a))
            mode={
              "path":str(path),"k_Mpc_inv":k,"rows":int(arr.shape[0]),
              "a_min":float(np.nanmin(a)),"a_max":amax,
              "ln_a_max":float(math.log(amax)) if amax>0 else None,
              "all_columns_finite":finite_all,
              "future_rows":int(np.count_nonzero(future)),
              "columns":names,
            }
            for key,cands in {
              "scalar_field":("delta_phi_smg (sync)","V_x_smg (sync)"),
              "scalar_field_prime":("delta_phi_prime_smg (sync)","V_x_prime_smg (sync)"),
              "scalar_field_prime_prime":("delta_phi_prime_prime_smg (sync)","V_x_prime_prime_smg (sync)"),
              "metric_psi":("psi",),
              "metric_phi":("phi",),
              "delta_cdm":("delta_cdm",),
              "delta_b":("delta_b",),
            }.items():
                y,used=col(names,arr,*cands)
                if y is not None:
                    mode[key]={"column":used,**amp_summary(a,y)}
            mode["endpoint_gate_pass"]=bool(rc==0 and finite_all and amax>1.01)
            all_ok=all_ok and mode["endpoint_gate_pass"]
            run["files"].append(mode)
        except Exception as exc:
            run["files"].append({"path":str(path),"status":"parse_error","error":str(exc)})
            all_ok=False
    run["status"]="ok" if all_ok else "incomplete_or_failed"
    run["all_endpoint_gates_pass"]=bool(all_ok)
    out["runs"][name]=run

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NATIVE_FUTURE_PERTURBATION_SUMMARY")
for name,r in out["runs"].items():
    compact={"status":r.get("status"),"returncode":r.get("returncode"),
             "all_endpoint_gates_pass":r.get("all_endpoint_gates_pass")}
    if r.get("files"):
        compact["modes"]=[{
          "k":m.get("k_Mpc_inv"),"a_max":m.get("a_max"),
          "finite":m.get("all_columns_finite"),
          "scalar_max_ratio":m.get("scalar_field",{}).get("max_abs_over_first_future"),
          "scalar_final_ratio":m.get("scalar_field",{}).get("final_abs_over_first_future"),
          "gate":m.get("endpoint_gate_pass")
        } for m in r["files"]]
    print(name,json.dumps(compact,sort_keys=True))
