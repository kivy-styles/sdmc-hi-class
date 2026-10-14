#!/usr/bin/env python3
"""Compare exact-tracker smooth-EFT CMB controls to same-branch coordinates.

The comparison isolates the response produced by the frozen lambda_e=20
tracker background together with a deliberately smooth scalar perturbation
closure.  It is a control layer only; the final SDMC early perturbation target
has c_phi,e^2 ~= 0.003 and is handled separately.
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
                    if ":" in tok and tok.split(":", 1)[0].isdigit():
                        named.append(tok.split(":", 1)[1])
                if len(named) >= 2:
                    cols = named
                continue
            rows.append([float(x) for x in s.split()])
    if not rows:
        raise RuntimeError(f"No data in {path}")
    if cols is None or len(cols) != len(rows[0]):
        cols = ["l", "TT", "EE", "TE", "BB", "phiphi", "Tphi", "Ephi"][:len(rows[0])]
    return cols, rows


def find_spectrum(prefix):
    cands = [prefix + "cl_lensed.dat", prefix + "cl.dat"]
    cands += sorted(glob.glob(prefix + "*_cl_lensed.dat"))
    cands += sorted(glob.glob(prefix + "*_cl.dat"))
    for p in cands:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(prefix)


def compare(label, coord_prefix, smooth_prefix):
    cp = find_spectrum(coord_prefix)
    sp = find_spectrum(smooth_prefix)
    cc, cr = read_table(cp)
    sc, sr = read_table(sp)
    ci = {x: i for i, x in enumerate(cc)}
    si = {x: i for i, x in enumerate(sc)}
    common = [x for x in ("TT", "EE", "TE", "phiphi") if x in ci and x in si]
    cby = {int(round(r[ci["l"]])): r for r in cr}
    sby = {int(round(r[si["l"]])): r for r in sr}
    ells = sorted(set(cby) & set(sby))

    outrows = []
    for ell in ells:
        a, b = cby[ell], sby[ell]
        row = {"branch": label, "ell": ell}
        for x in common:
            va, vb = a[ci[x]], b[si[x]]
            row[f"coord_{x}"] = va
            row[f"tracker_smooth_{x}"] = vb
            row[f"frac_{x}"] = (vb-va)/va if va != 0 else ""
        if all(x in common for x in ("TT", "EE", "TE")):
            den = math.sqrt(abs(a[ci["TT"]]*a[ci["EE"]]))
            row["delta_TE_norm"] = (b[si["TE"]]-a[ci["TE"]])/den if den > 0 else ""
        outrows.append(row)

    print(f"\n{label}: tracker-smooth minus coordinate")
    for lo, hi in [(2,29),(30,499),(500,1499),(1500,3000)]:
        sub = [r for r in outrows if lo <= r["ell"] <= hi]
        parts = []
        for x in ("TT","EE","phiphi"):
            key = f"frac_{x}"
            vals = [abs(float(r[key])) for r in sub if key in r and r[key] != ""]
            if vals:
                rms = math.sqrt(sum(v*v for v in vals)/len(vals))
                parts.append(f"{x}: max={100*max(vals):.3f}%, rms={100*rms:.3f}%")
        vals = [abs(float(r["delta_TE_norm"])) for r in sub if r.get("delta_TE_norm","") != ""]
        if vals:
            rms = math.sqrt(sum(v*v for v in vals)/len(vals))
            parts.append(f"TE/sqrt(TT EE): max={100*max(vals):.3f}%, rms={100*rms:.3f}%")
        print(f"  ell {lo:4d}-{hi:4d}: " + "; ".join(parts))

    # TT peak locations in broad, non-overlapping acoustic bands.
    if "TT" in common:
        bands = [(100,350),(350,650),(650,950),(950,1250),(1250,1600),(1600,2000),(2000,2500)]
        print("  TT broad-band peak locations (coordinate -> tracker-smooth):")
        for j,(lo,hi) in enumerate(bands,1):
            ls = [l for l in ells if lo <= l <= hi]
            if not ls:
                continue
            lc = max(ls, key=lambda l: cby[l][ci["TT"]])
            lsmo = max(ls, key=lambda l: sby[l][si["TT"]])
            print(f"    peak {j}: ell {lc} -> {lsmo}  (Delta ell={lsmo-lc:+d})")
    return outrows


def main():
    allrows = []
    allrows += compare("Kp", "output/sdmc_kp_coordinate_", "output/sdmc_kp_tracker_smooth_")
    allrows += compare("NKp-v2", "output/sdmc_nkp_v2_coordinate_", "output/sdmc_nkp_v2_tracker_smooth_")
    out = "output/sdmc_tracker_smooth_cl_comparison.csv"
    keys = []
    for r in allrows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(out,"w",newline="",encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(allrows)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
