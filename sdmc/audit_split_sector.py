#!/usr/bin/env python3
"""Audit the SDMC split-sector Boltzmann runs.

Compares each split-sector completion to its clustered-tracker-only reference.
No external Python packages are required.
"""
from __future__ import annotations
import csv, glob, math, os, re
from pathlib import Path

CASES = [
    ("Kp_core", "output/sdmc_kp_core_split_sector_", "output/sdmc_kp_tracker_clustered_cs0003_"),
    ("NKp_v2", "output/sdmc_nkp_v2_split_sector_", "output/sdmc_nkp_v2_tracker_clustered_cs0003_"),
]


def data(path):
    rows=[]
    with open(path,encoding="utf-8") as f:
        for line in f:
            s=line.strip()
            if s and not s.startswith("#"):
                rows.append([float(x) for x in s.split()])
    return rows


def first(pattern):
    xs=sorted(glob.glob(pattern))
    if not xs:
        raise FileNotFoundError(pattern)
    return xs[0]


def sigma8(path):
    rows=data(path)
    vals=[]
    for r in rows:
        k,p=r[0],r[1]
        x=8.*k
        W=1. if abs(x)<1e-10 else 3.*(math.sin(x)-x*math.cos(x))/x**3
        vals.append((math.log(k), k**3*p*W*W/(2.*math.pi**2)))
    s=0.
    for a,b in zip(vals[:-1],vals[1:]):
        s += .5*(a[1]+b[1])*(b[0]-a[0])
    return math.sqrt(max(s,0.))


def nearest(rows,ell,col):
    return min(rows,key=lambda r:abs(r[0]-ell))[col]


def header_cols(path):
    hdr=None
    with open(path,encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") and "1:" in line:
                hdr=line[1:].strip()
    if hdr is None:
        return []
    marks=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    out=[]
    for i,m in enumerate(marks):
        end=marks[i+1].start() if i+1<len(marks) else len(hdr)
        out.append(hdr[m.end():end].strip())
    return out


def bg_stats(path):
    cols=header_cols(path)
    rows=data(path)
    idx={c:i for i,c in enumerate(cols)}
    out={}
    aliases={
        "F":"M*^2_smg",
        "alphaM":"M2_running_smg",
        "D":"kin (D)",
        "cs2":"c_s^2",
        "rho_fld":"(.)rho_fld",
        "w_fld":"(.)w_fld",
    }
    for key,col in aliases.items():
        if col not in idx:
            continue
        j=idx[col]
        arr=[r[j] for r in rows]
        out[key+"_min"]=min(arr)
        out[key+"_max"]=max(arr)
        # background z column is first in CLASS output
        out[key+"_z0"]=min(rows,key=lambda r:abs(r[0]))[j]
    return out


def main():
    fields=[
        "case","status","sigma8_z0","sigma8_z10","growth_0_over_10",
        "phiphi_ratio_l100","phiphi_ratio_l500","phiphi_ratio_l1000","phiphi_ratio_l1500",
        "TT_ratio_l200","TT_ratio_l1000","TT_ratio_l2000",
        "EE_ratio_l200","EE_ratio_l1000","EE_ratio_l2000",
        "F_z0","F_min","F_max","alphaM_z0","alphaM_min","alphaM_max",
        "D_min","cs2_min","rho_fld_z0","rho_fld_max","w_fld_min","w_fld_max","error"
    ]
    out=[]
    for name,prefix,refprefix in CASES:
        row={k:"" for k in fields}; row["case"]=name
        try:
            cl=data(first(prefix+"*_cl_lensed.dat"))
            ref=data(first(refprefix+"*_cl_lensed.dat"))
            pk0=first(prefix+"*_z1_pk.dat")
            pk10=first(prefix+"*_z7_pk.dat")
            s0=sigma8(pk0); s10=sigma8(pk10)
            row["sigma8_z0"]=s0; row["sigma8_z10"]=s10; row["growth_0_over_10"]=s0/s10
            for ell in (100,500,1000,1500):
                row[f"phiphi_ratio_l{ell}"]=nearest(cl,ell,5)/nearest(ref,ell,5)
            for ell in (200,1000,2000):
                row[f"TT_ratio_l{ell}"]=nearest(cl,ell,1)/nearest(ref,ell,1)
                row[f"EE_ratio_l{ell}"]=nearest(cl,ell,2)/nearest(ref,ell,2)
            bg=bg_stats(first(prefix+"*_background.dat"))
            for k,v in bg.items():
                if k in row: row[k]=v
            row["status"]="OK"
        except Exception as exc:
            row["status"]="FAIL"; row["error"]=repr(exc)
        out.append(row)

    Path("output").mkdir(exist_ok=True)
    path=Path("output/sdmc_split_sector_audit.csv")
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(out)
    print(path)
    for r in out:
        print(r)

if __name__=="__main__":
    main()
