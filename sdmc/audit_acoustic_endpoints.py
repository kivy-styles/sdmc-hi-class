#!/usr/bin/env python3
"""Reconcile SDMC manuscript acoustic rulers with hi_class endpoint choices.

CLASS/HyRec reports the sound horizon at the maximum of the visibility function
and at its internally determined baryon-drag epoch.  The SDMC manuscript uses
working benchmark endpoints z_*=1090 and z_d~=1059.  Since the background
output tabulates the comoving sound horizon as a function of z, this script
compares the two conventions without changing lambda_e or any branch physics.
"""
from __future__ import annotations

import glob
import math

TARGETS = {
    "Kp": {"prefix":"output/sdmc_kp_full_background_", "rs":143.779000, "rd":146.518300},
    "NKp-v2": {"prefix":"output/sdmc_nkp_v2_full_background_", "rs":142.688139, "rd":145.395978},
}


def bgfile(prefix):
    c=sorted(glob.glob(prefix+"*_background.dat"))
    if not c: c=sorted(glob.glob(prefix+"background.dat"))
    if not c: raise FileNotFoundError(prefix)
    return c[0]


def load(path):
    rows=[]
    with open(path,"r",encoding="utf-8") as f:
        for line in f:
            s=line.strip()
            if not s or s.startswith("#"): continue
            q=[float(x) for x in s.split()]
            rows.append((q[0],q[7]))  # z, comov.snd.hrz.
    rows.sort()
    return rows


def interp(rows,z):
    if z<=rows[0][0]: return rows[0][1]
    if z>=rows[-1][0]: return rows[-1][1]
    lo=0; hi=len(rows)-1
    while hi-lo>1:
        m=(lo+hi)//2
        if rows[m][0] <= z: lo=m
        else: hi=m
    z0,r0=rows[lo]; z1,r1=rows[hi]
    t=(z-z0)/(z1-z0)
    return r0+t*(r1-r0)


def invert(rows,target):
    # Around recombination r_s decreases monotonically with z.
    rr=sorted((r,z) for z,r in rows)
    if target<=rr[0][0]: return rr[0][1]
    if target>=rr[-1][0]: return rr[-1][1]
    lo=0; hi=len(rr)-1
    while hi-lo>1:
        m=(lo+hi)//2
        if rr[m][0] <= target: lo=m
        else: hi=m
    r0,z0=rr[lo]; r1,z1=rr[hi]
    t=(target-r0)/(r1-r0)
    return z0+t*(z1-z0)


def main():
    for name,cfg in TARGETS.items():
        rows=load(bgfile(cfg["prefix"]))
        rs1090=interp(rows,1090.0)
        rd1059=interp(rows,1059.0)
        zrs=invert(rows,cfg["rs"])
        zrd=invert(rows,cfg["rd"])
        print(f"\n{name}")
        print(f"  r_s(z=1090) = {rs1090:.6f} Mpc; manuscript target = {cfg['rs']:.6f}; delta={rs1090-cfg['rs']:+.6f}")
        print(f"  r_s(z=1059) = {rd1059:.6f} Mpc; manuscript drag target = {cfg['rd']:.6f}; delta={rd1059-cfg['rd']:+.6f}")
        print(f"  exact tabulated z giving manuscript r_s target = {zrs:.6f}")
        print(f"  exact tabulated z giving manuscript r_d target = {zrd:.6f}")
    print("\nNo model parameter was retuned in this audit.")


if __name__ == "__main__":
    main()
