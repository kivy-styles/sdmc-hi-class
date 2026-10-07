#!/usr/bin/env python3
from pathlib import Path
import sys,re,json,numpy as np

EPS=[0.0,1.90e-4,1.92e-4,2.00e-4,3.00e-4,5.00e-4,1.00e-3]
DPHI=[-0.75,-0.50,-0.25,0.0]
DVEL=[0.0,0.5,1.0,1.5,2.0]

def read_bg(path):
    ls=Path(path).read_text().splitlines()
    h=[x for x in ls if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",h)); ns=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(h)
        ns.append(h[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(ns)}

def prepare():
    base=Path("output/linear_cov_free.ini").read_text()
    meta=[]
    for ie,eps in enumerate(EPS):
      for ip,dp in enumerate(DPHI):
       for iv,dv in enumerate(DVEL):
        key=f"e{ie:02d}_p{ip:02d}_v{iv:02d}"
        s=re.sub(r"(?m)^\s*parameters_smg\s*=.*$",
                 f"parameters_smg = 0.0, {dp:.17g}, {dv:.17g}, {eps:.17g}",base)
        s=re.sub(r"(?m)^\s*root\s*=.*$",f"root = output/barea_{key}_",s)
        s="\n".join(line for line in s.splitlines() if not re.match(r"^\s*output\s*=",line))
        s=re.sub(r"(?m)^\s*lensing\s*=.*$","lensing = no",s)
        s=re.sub(r"(?m)^\s*write_thermodynamics\s*=.*$","write_thermodynamics = no",s)
        s += "\nskip_stability_tests_smg = yes\nreio_parametrization = reio_none\n"
        Path(f"output/barea_{key}.ini").write_text(s+"\n")
        meta.append({"key":key,"epsilon":eps,"delta_phi":dp,"delta_velocity_fraction":dv})
    Path("output/null_g2_basin_area_points.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")
    print("NULL_G2_BASIN_AREA_PREPARED",len(meta),flush=True)

def analyze():
    meta=json.loads(Path("output/null_g2_basin_area_points.json").read_text())
    ref=read_bg("output/barea_e00_p03_v00_00_background.dat")
    Nr=-np.log1p(ref["z"]); oo=np.argsort(Nr)
    rows=[]
    for q in meta:
        p=Path("output")/f"barea_{q['key']}_00_background.dat"
        rcpath=Path("output")/f"barea_{q['key']}.rc"
        rc=int(rcpath.read_text()) if rcpath.exists() else 999
        r=dict(q,returncode=rc,success=(rc==0 and p.exists()))
        if r["success"]:
            d=read_bg(p); N=-np.log1p(d["z"]); o=np.argsort(N); z=d["z"][o]; m=z<=100
            Href=np.interp(N[o],Nr[oo],ref["H [1/Mpc]"][oo])
            Fref=np.interp(N[o],Nr[oo],ref["M*^2_smg"][oo])
            r.update({
              "min_D_z100":float(np.min(d["kin (D)"][o][m])),
              "min_cs2_z100":float(np.min(d["c_s^2"][o][m])),
              "max_cs2_z100":float(np.max(d["c_s^2"][o][m])),
              "max_abs_dH_H_vs_base_z100":float(np.max(np.abs(d["H [1/Mpc]"][o][m]/Href[m]-1))),
              "max_abs_dF_F_vs_base_z100":float(np.max(np.abs(d["M*^2_smg"][o][m]/Fref[m]-1))),
            })
            r["stable"]=bool(r["min_D_z100"]>0 and r["min_cs2_z100"]>0 and r["max_cs2_z100"]<=1.0+1e-8)
        else:
            r["stable"]=False
        rows.append(r)
    summary=[]
    for eps in EPS:
        rr=[r for r in rows if r["epsilon"]==eps]
        success=[r for r in rr if r["success"]]
        stable=[r for r in rr if r["stable"]]
        grad=[r for r in success if r.get("min_cs2_z100",1)>0]
        ghost=[r for r in success if r.get("min_D_z100",1)>0]
        summary.append({
          "epsilon":eps,
          "n_points":len(rr),
          "n_success":len(success),
          "n_stable":len(stable),
          "stable_fraction":len(stable)/len(rr),
          "n_positive_cs2":len(grad),
          "n_positive_D":len(ghost),
          "worst_min_cs2":min((r["min_cs2_z100"] for r in success),default=None),
          "worst_min_D":min((r["min_D_z100"] for r in success),default=None),
          "stable_points":[[r["delta_phi"],r["delta_velocity_fraction"]] for r in stable]
        })
    out={
      "status":"reionization-free finite-amplitude null-G2 basin-area audit",
      "deformation":"Delta G2 = epsilon*(X-Xstar(phi))^4/Xstar(phi)^3",
      "grid":{"delta_phi":DPHI,"delta_velocity_fraction":DVEL,"n_points_per_epsilon":len(DPHI)*len(DVEL)},
      "epsilon_values":EPS,
      "summary":summary,
      "rows":rows,
      "qualification":[
        "This sampled rectangle is a diagnostic basin measure, not a fundamental initial-condition prior.",
        "Thermodynamics/reionization is disabled to isolate covariant background scalar stability.",
        "The accepted late300 trajectory is not retuned."
      ]
    }
    Path("output/null_g2_basin_area.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("NULL_G2_BASIN_AREA",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__":
    if len(sys.argv)!=2 or sys.argv[1] not in ("prepare","analyze"):
        raise SystemExit("usage: run_null_g2_basin_area.py prepare|analyze")
    {"prepare":prepare,"analyze":analyze}[sys.argv[1]]()
