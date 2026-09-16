#!/usr/bin/env python3
"""Compare fully-dynamic and quasi-static NKp-v2 No-Slip perturbations.

The background and EFT functions are identical in all rows. This is a solver
regime diagnostic, not a physical likelihood. A forced-QS failure is retained
as a diagnostic rather than silently interpreted as a model rejection.
"""
from __future__ import annotations
import csv, glob, math, os, subprocess
from pathlib import Path

TEMPLATE=Path("sdmc/config/nkp_v2_full_clustered_noslip.ini")
REF=Path("output/sdmc_nkp_v2_tracker_clustered_cs0003_00_cl_lensed.dat")
OUT=Path("output/sdmc_nkp_v2_noslip_qs_scan.csv")
METHODS=["fully_dynamic","automatic","quasi_static","quasi_static_debug"]


def replace_line(text,key,value):
    lines=text.splitlines(); out=[]; n=0
    for line in lines:
        if line.strip().startswith(key+" ="):
            out.append(f"{key} = {value}"); n+=1
        else: out.append(line)
    if n!=1: raise RuntimeError(f"expected one {key}, found {n}")
    return "\n".join(out)+"\n"


def table(path):
    rows=[]
    with open(path,encoding="utf-8") as f:
        for line in f:
            s=line.strip()
            if s and not s.startswith("#"): rows.append([float(x) for x in s.split()])
    return rows


def find1(pattern):
    x=sorted(glob.glob(pattern))
    if len(x)!=1: raise RuntimeError(f"{pattern}: found {len(x)}")
    return x[0]


def sigma8(rows):
    v=[]
    for r in rows:
        k,p=r[0],r[1]; x=8*k
        w=1. if abs(x)<1e-8 else 3*(math.sin(x)-x*math.cos(x))/x**3
        v.append((math.log(k),k**3*p*w*w/(2*math.pi**2)))
    q=sum(.5*(b[1]+a[1])*(b[0]-a[0]) for a,b in zip(v[:-1],v[1:]))
    return math.sqrt(max(q,0.))


def nearest(rows,ell,col): return min(rows,key=lambda r:abs(r[0]-ell))[col]

def clean(prefix):
    for p in glob.glob(prefix+"*"):
        try: os.remove(p)
        except OSError: pass


def main():
    base=TEMPLATE.read_text(encoding="utf-8")
    ref=table(REF)
    fields=["method","status","sigma8_z0","sigma8_z10","growth_0_over_10",
            "phiphi_ratio_l100","phiphi_ratio_l500","phiphi_ratio_l1000","phiphi_ratio_l1500",
            "TT_ratio_l200","TT_ratio_l1000","TT_ratio_l2000","error"]
    rows=[]
    for method in METHODS:
        prefix=f"output/sdmc_nkp_v2_noslip_qs_{method}_"
        clean(prefix)
        ini=Path(f"output/sdmc_nkp_v2_noslip_qs_{method}.ini")
        log=Path(f"output/sdmc_nkp_v2_noslip_qs_{method}.log")
        text=replace_line(base,"method_qs_smg",method)
        text=replace_line(text,"z_pk","0,10")
        text=replace_line(text,"root",prefix)
        for k in ("input_verbose","background_verbose","thermodynamics_verbose","output_verbose"):
            text=replace_line(text,k,"0")
        text=replace_line(text,"perturbations_verbose","2" if "quasi_static" in method else "0")
        ini.write_text(text,encoding="utf-8")
        row={k:"" for k in fields}; row["method"]=method
        try:
            cp=subprocess.run(["./class",str(ini)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=240)
            log.write_text(cp.stdout,encoding="utf-8")
            if cp.returncode:
                row["status"]="FAIL"
                tail=[x for x in cp.stdout.strip().splitlines() if x.strip()][-18:]
                row["error"]=" | ".join(tail)[:3000]
                print(method,"FAIL",row["error"]); rows.append(row); continue
            cl=table(find1(prefix+"*_cl_lensed.dat")); p0=table(find1(prefix+"*_z1_pk.dat")); p10=table(find1(prefix+"*_z2_pk.dat"))
            s0,s10=sigma8(p0),sigma8(p10)
            row.update(status="OK",sigma8_z0=f"{s0:.12g}",sigma8_z10=f"{s10:.12g}",growth_0_over_10=f"{s0/s10:.12g}")
            for e in (100,500,1000,1500): row[f"phiphi_ratio_l{e}"]=f"{nearest(cl,e,5)/nearest(ref,e,5):.12g}"
            for e in (200,1000,2000): row[f"TT_ratio_l{e}"]=f"{nearest(cl,e,1)/nearest(ref,e,1):.12g}"
            print(method,"OK","sigma8",s0,"growth",s0/s10,"phi1000/ref",row["phiphi_ratio_l1000"])
            rows.append(row)
        except subprocess.TimeoutExpired as exc:
            row["status"]="TIMEOUT"; row["error"]=f"timeout after {exc.timeout} s"
            print(method,"TIMEOUT"); rows.append(row)
        finally:
            clean(prefix)
            try: ini.unlink()
            except OSError: pass
    with OUT.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    print("Wrote",OUT)

if __name__=="__main__": main()
