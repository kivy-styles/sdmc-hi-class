#!/usr/bin/env python3
from pathlib import Path
import re, subprocess, json
import numpy as np

BASE=Path("input_artifact/unpacked/candidate.ini")
if not BASE.exists():
    raise SystemExit("candidate.ini missing")
base=BASE.read_text()

cases=[
    ("a1e7","1.e-7","1.e-7","1.e6","1.e-7"),
    ("a3e8","3.e-8","3.e-8","3.333333333e7","3.e-8"),
    ("a1e8","1.e-8","1.e-8","1.e7","1.e-8"),
    ("a3e9","3.e-9","3.e-9","3.333333333e8","3.e-9"),
    ("a1e9","1.e-9","1.e-9","1.e9","1.e-9"),
    ("a1e10","1.e-10","1.e-10","1.e10","1.e-10"),
]

def set_value(text,key,val):
    pat=rf"(?m)^\s*{re.escape(key)}\s*=.*$"
    out,n=re.subn(pat,f"{key} = {val}",text,count=1)
    if n!=1:
        raise RuntimeError(f"{key}: replacement count {n}")
    return out

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return names,a

rows=[]
for tag,aini,aqs,zref,amin in cases:
    root=f"output/stabrob_{tag}_"
    s=base
    for key,val in [
        ("a_ini_over_a_today_default",aini),
        ("a_ini_test_qs_smg",aqs),
        ("pert_ic_ini_z_ref_smg",zref),
        ("a_min_stability_test_smg",amin),
        ("root",root),
    ]:
        s=set_value(s,key,val)
    ip=Path(f"output/stabrob_{tag}.ini")
    ip.write_text(s)
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    Path(f"output/stabrob_{tag}.log").write_text(cp.stdout)
    rec={"case":tag,"returncode":int(cp.returncode)}
    p=Path(root+"00_background.dat")
    if cp.returncode==0 and p.exists():
        names,a=table(p); d={n:a[:,i] for i,n in enumerate(names)}
        cs=d["c_s^2"]; D=d["kin (D)"]; z=d["z"]
        ii=int(np.argmin(cs)); lo=max(0,ii-4); hi=min(len(cs),ii+5)
        rec.update(
            min_cs2=float(cs[ii]),
            z_at_min_cs2=float(z[ii]),
            D_at_min_cs2=float(D[ii]),
            min_D=float(np.min(D)),
            max_cs2=float(np.max(cs)),
            max_abs_noslip=float(np.max(np.abs(d["braiding_smg"]+2*d["M2_running_smg"]))),
            neighborhood=[
                {"z":float(z[j]),"cs2":float(cs[j]),"D":float(D[j])}
                for j in range(lo,hi)
            ],
        )
    else:
        rec["tail"]=cp.stdout[-800:].replace("\n"," | ")
    print("EDGE024_QBEST_STABROB",json.dumps(rec,sort_keys=True),flush=True)
    rows.append(rec)

Path("output/edge024_qbest_stability_robustness.json").write_text(json.dumps(rows,indent=2))
good=[r for r in rows if r.get("returncode")==0 and "min_cs2" in r]
if len(good)!=len(cases):
    raise SystemExit("one or more robustness cases failed")
if any(r["min_cs2"]<=0 or r["min_D"]<=0 or r["max_cs2"]>1 for r in good):
    raise SystemExit("stability robustness failed")
vals=np.array([r["min_cs2"] for r in good],float)
zs=np.array([r["z_at_min_cs2"] for r in good],float)
summary={
    "n_good":len(good),
    "min_of_min_cs2":float(vals.min()),
    "max_of_min_cs2":float(vals.max()),
    "relative_spread":float((vals.max()-vals.min())/vals.mean()),
    "z_min_range":[float(zs.min()),float(zs.max())],
}
print("EDGE024_QBEST_STABROB_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
Path("output/edge024_qbest_stability_robustness_summary.json").write_text(json.dumps(summary,indent=2))
