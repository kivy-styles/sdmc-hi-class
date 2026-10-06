#!/usr/bin/env python3
"""Compare coordinate, smooth-tracker, and clustered-tracker CMB spectra.

For each SDMC branch the script reports two differentials:
  (1) clustered minus coordinate: total early-tracker response;
  (2) clustered minus smooth: the isolated response to lowering the scalar
      propagation speed to c_phi^2=0.003 at fixed exact tracker background.
"""
from __future__ import annotations

import csv
import glob
import math
import os


def read_table(path):
    cols = None
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                named = []
                for tok in s.lstrip("#").split():
                    if ":" in tok and tok.split(":",1)[0].isdigit():
                        named.append(tok.split(":",1)[1])
                if len(named) >= 2:
                    cols = named
                continue
            rows.append([float(x) for x in s.split()])
    if not rows:
        raise RuntimeError(path)
    if cols is None or len(cols) != len(rows[0]):
        cols = ["l","TT","EE","TE","BB","phiphi","Tphi","Ephi"][:len(rows[0])]
    return cols, rows


def find_spectrum(prefix):
    cands = [prefix+"cl_lensed.dat", prefix+"cl.dat"]
    cands += sorted(glob.glob(prefix+"*_cl_lensed.dat"))
    cands += sorted(glob.glob(prefix+"*_cl.dat"))
    for p in cands:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(prefix)


def load(prefix):
    p = find_spectrum(prefix)
    c,r = read_table(p)
    i = {x:j for j,x in enumerate(c)}
    return p, i, {int(round(q[i["l"]])):q for q in r}


def differential(name, A, B, fields=("TT","EE","phiphi")):
    pa, ia, a = A
    pb, ib, b = B
    ells = sorted(set(a)&set(b))
    print(f"\n{name}")
    for lo,hi in [(2,29),(30,499),(500,1499),(1500,3000)]:
        ls=[l for l in ells if lo<=l<=hi]
        parts=[]
        for x in fields:
            if x not in ia or x not in ib:
                continue
            vals=[]
            for l in ls:
                va=a[l][ia[x]]; vb=b[l][ib[x]]
                if va != 0:
                    vals.append(abs((vb-va)/va))
            if vals:
                rms=math.sqrt(sum(v*v for v in vals)/len(vals))
                parts.append(f"{x}: max={100*max(vals):.3f}%, rms={100*rms:.3f}%")
        if all(x in ia and x in ib for x in ("TT","EE","TE")):
            vals=[]
            for l in ls:
                den=math.sqrt(abs(a[l][ia["TT"]]*a[l][ia["EE"]]))
                if den>0:
                    vals.append(abs((b[l][ib["TE"]]-a[l][ia["TE"]])/den))
            if vals:
                rms=math.sqrt(sum(v*v for v in vals)/len(vals))
                parts.append(f"TE/sqrt(TT EE): max={100*max(vals):.3f}%, rms={100*rms:.3f}%")
        print(f"  ell {lo:4d}-{hi:4d}: "+"; ".join(parts))

    if "TT" in ia and "TT" in ib:
        bands=[(100,350),(350,650),(650,950),(950,1250),(1250,1600),(1600,2000),(2000,2500)]
        print("  TT broad-band peaks A -> B:")
        for j,(lo,hi) in enumerate(bands,1):
            ls=[l for l in ells if lo<=l<=hi]
            if not ls: continue
            la=max(ls,key=lambda l:a[l][ia["TT"]])
            lb=max(ls,key=lambda l:b[l][ib["TT"]])
            print(f"    peak {j}: ell {la} -> {lb} (Delta ell={lb-la:+d})")


def rows_for_branch(label, coord, smooth, clustered):
    pc,ic,c=coord; ps,is_,s=smooth; pk,ik,k=clustered
    ells=sorted(set(c)&set(s)&set(k))
    rows=[]
    for l in ells:
        row={"branch":label,"ell":l}
        for x in ("TT","EE","TE","phiphi"):
            if x in ic and x in is_ and x in ik:
                vc=c[l][ic[x]]; vs=s[l][is_[x]]; vk=k[l][ik[x]]
                row[f"coord_{x}"]=vc
                row[f"smooth_{x}"]=vs
                row[f"clustered_{x}"]=vk
                row[f"frac_clustered_minus_coord_{x}"]=(vk-vc)/vc if vc!=0 else ""
                row[f"frac_clustered_minus_smooth_{x}"]=(vk-vs)/vs if vs!=0 else ""
        rows.append(row)
    return rows


def main():
    branches={
        "Kp": (
            load("output/sdmc_kp_coordinate_"),
            load("output/sdmc_kp_tracker_smooth_"),
            load("output/sdmc_kp_tracker_clustered_cs0003_"),
        ),
        "NKp-v2": (
            load("output/sdmc_nkp_v2_coordinate_"),
            load("output/sdmc_nkp_v2_tracker_smooth_"),
            load("output/sdmc_nkp_v2_tracker_clustered_cs0003_"),
        ),
    }
    allrows=[]
    for label,(coord,smooth,clustered) in branches.items():
        differential(f"{label}: clustered(c_s^2=0.003) minus coordinate", coord, clustered)
        differential(f"{label}: clustered(c_s^2=0.003) minus smooth tracker", smooth, clustered)
        allrows += rows_for_branch(label,coord,smooth,clustered)

    out="output/sdmc_tracker_clustered_cl_comparison.csv"
    keys=[]
    for r in allrows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(out,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=keys)
        w.writeheader(); w.writerows(allrows)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
