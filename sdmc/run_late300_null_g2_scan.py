#!/usr/bin/env python3
from pathlib import Path
import sys,re,json,numpy as np

EPS=[0.0,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1.0]
POINTS=[
 ("base",0.0,0.0),
 ("A_bad",-0.75,0.50),("A_good",-0.75,0.60),
 ("B1_bad",-0.50,0.60),("B2_bad",-0.50,1.00),("B_good",-0.50,1.20),
 ("C_good",0.0,1.50),("C1_bad",0.0,1.80),("C2_bad",0.0,2.00),
]
BAD={"A_bad","B1_bad","B2_bad","C1_bad","C2_bad"}

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
        for tag,dp,dv in POINTS:
            key=f"e{ie:02d}_{tag}"
            s=re.sub(r"(?m)^\s*parameters_smg\s*=.*$",
                     f"parameters_smg = 0.0, {dp:.17g}, {dv:.17g}, {eps:.17g}",base)
            s=re.sub(r"(?m)^\s*root\s*=.*$",f"root = output/null_{key}_",s)
            s="\n".join(line for line in s.splitlines() if not re.match(r"^\s*output\s*=",line))
            s=re.sub(r"(?m)^\s*lensing\s*=.*$","lensing = no",s)
            s=re.sub(r"(?m)^\s*write_thermodynamics\s*=.*$","write_thermodynamics = no",s)
            s += "\nskip_stability_tests_smg = yes\n"
            Path(f"output/null_{key}.ini").write_text(s+"\n")
            meta.append({"key":key,"epsilon":eps,"point":tag,"delta_phi":dp,"delta_velocity_fraction":dv})
    Path("output/null_g2_points.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")
    print("NULL_G2_PREPARED",len(meta))

def analyze():
    meta=json.loads(Path("output/null_g2_points.json").read_text())
    ref=read_bg("output/null_e00_base_00_background.dat")
    Nr=-np.log1p(ref["z"]); oo=np.argsort(Nr)
    rows=[]
    for q in meta:
        p=Path("output")/f"null_{q['key']}_00_background.dat"
        rcpath=Path("output")/f"null_{q['key']}.rc"
        rc=int(rcpath.read_text()) if rcpath.exists() else 999
        r=dict(q,returncode=rc,success=(rc==0 and p.exists()))
        if not r["success"]:
            rows.append(r); continue
        d=read_bg(p); N=-np.log1p(d["z"]); o=np.argsort(N); z=d["z"][o]; m=z<=100
        Href=np.interp(N[o],Nr[oo],ref["H [1/Mpc]"][oo])
        Fref=np.interp(N[o],Nr[oo],ref["M*^2_smg"][oo])
        r.update({
          "min_D_z100":float(np.min(d["kin (D)"][o][m])),
          "min_cs2_z100":float(np.min(d["c_s^2"][o][m])),
          "max_cs2_z100":float(np.max(d["c_s^2"][o][m])),
          "max_abs_dH_H_vs_base_z100":float(np.max(np.abs(d["H [1/Mpc]"][o][m]/Href[m]-1))),
          "max_abs_dF_F_vs_base_z100":float(np.max(np.abs(d["M*^2_smg"][o][m]/Fref[m]-1))),
          "max_abs_noslip_z100":float(np.max(np.abs(d["braiding_smg"][o][m]+2*d["M2_running_smg"][o][m]))),
        })
        r["stable"]=bool(r["min_D_z100"]>0 and r["min_cs2_z100"]>0 and r["max_cs2_z100"]<=1.0+1e-8)
        rows.append(r)

    summary=[]
    for ie,eps in enumerate(EPS):
        rr=[r for r in rows if r["epsilon"]==eps]
        base=next(r for r in rr if r["point"]=="base")
        bad=[r for r in rr if r["point"] in BAD]
        good=[r for r in rr if r["point"] not in BAD and r["point"]!="base"]
        summary.append({
          "epsilon":eps,
          "base_stable":base.get("stable",False),
          "base_max_abs_dH_H_z100":base.get("max_abs_dH_H_vs_base_z100"),
          "base_max_abs_dF_F_z100":base.get("max_abs_dF_F_vs_base_z100"),
          "base_min_cs2_z100":base.get("min_cs2_z100"),
          "bad_points_stabilized":sum(bool(r.get("stable")) for r in bad),
          "bad_points_total":len(bad),
          "good_neighbors_still_stable":sum(bool(r.get("stable")) for r in good),
          "good_neighbors_total":len(good),
          "worst_min_cs2_bad":min((r.get("min_cs2_z100",float("nan")) for r in bad),default=float("nan")),
          "all_original_bad_stable":all(bool(r.get("stable")) for r in bad),
        })
    passing=[s for s in summary if s["all_original_bad_stable"] and s["base_stable"]
             and (s["base_max_abs_dH_H_z100"] or 0)<1e-8 and (s["base_max_abs_dF_F_z100"] or 0)<1e-8]
    out={
      "status":"exploratory null-G2 off-trajectory stabilizer scan",
      "deformation":"Delta G2 = epsilon*(X-Xstar(phi))^4/Xstar(phi)^3",
      "on_trajectory_property":"Delta G2 and all derivatives used through total order 3 vanish at X=Xstar; target late300 should be unchanged",
      "epsilon_values":EPS,"rows":rows,"summary":summary,
      "smallest_strict_passing_epsilon":passing[0]["epsilon"] if passing else None,
      "qualification":"diagnostic of covariant-completion non-uniqueness; not adopted into accepted SDMC action"
    }
    Path("output/null_g2_stabilizer_scan.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("NULL_G2_SCAN_SUMMARY",json.dumps(out,sort_keys=True))

if __name__=="__main__":
    if len(sys.argv)!=2 or sys.argv[1] not in ("prepare","analyze"):
        raise SystemExit("usage: run_late300_null_g2_scan.py prepare|analyze")
    {"prepare":prepare,"analyze":analyze}[sys.argv[1]]()
