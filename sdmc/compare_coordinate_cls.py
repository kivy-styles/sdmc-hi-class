#!/usr/bin/env python3
"""Compare two CLASS CMB spectra without external Python dependencies.

This is a coordinate-control diagnostic only.  It does not claim to be the
full SDMC result; it isolates the response to the different frozen physical
matter densities before the tracker/No-Slip/terminal sectors are switched on.
"""
from __future__ import annotations

import csv
import math
import os
import sys


def read_class_table(path: str):
    columns = None
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                # CLASS header usually contains tokens such as 1:l 2:TT 3:EE ...
                toks = s.lstrip("#").split()
                named = []
                for tok in toks:
                    if ":" in tok and tok.split(":", 1)[0].isdigit():
                        named.append(tok.split(":", 1)[1])
                if len(named) >= 2:
                    columns = named
                continue
            vals = [float(x) for x in s.split()]
            rows.append(vals)

    if not rows:
        raise RuntimeError(f"No numerical rows in {path}")
    if columns is None or len(columns) != len(rows[0]):
        # Conservative fallback for standard scalar CMB CLASS output.
        default = ["l", "TT", "EE", "TE", "BB", "phiphi", "Tphi", "Ephi"]
        columns = default[: len(rows[0])]
    return columns, rows


def index_map(columns):
    return {name: i for i, name in enumerate(columns)}


def find_spectrum(prefix: str):
    candidates = [
        prefix + "cl_lensed.dat",
        prefix + "cl.dat",
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    raise FileNotFoundError("No CLASS C_l file found among: " + ", ".join(candidates))


def main():
    kp_prefix = sys.argv[1] if len(sys.argv) > 1 else "output/sdmc_kp_coordinate_"
    v2_prefix = sys.argv[2] if len(sys.argv) > 2 else "output/sdmc_nkp_v2_coordinate_"
    out = sys.argv[3] if len(sys.argv) > 3 else "output/sdmc_coordinate_cl_comparison.csv"

    kp_path = find_spectrum(kp_prefix)
    v2_path = find_spectrum(v2_prefix)
    ck, rk = read_class_table(kp_path)
    cv, rv = read_class_table(v2_path)
    ik, iv = index_map(ck), index_map(cv)

    common = [x for x in ("TT", "EE", "TE", "phiphi") if x in ik and x in iv]
    by_l_k = {int(round(r[ik["l"]])): r for r in rk}
    by_l_v = {int(round(r[iv["l"]])): r for r in rv}
    ells = sorted(set(by_l_k) & set(by_l_v))

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        fields = ["ell"]
        for x in common:
            fields += [f"Kp_{x}", f"NKp_v2_{x}", f"delta_{x}", f"frac_{x}"]
        if all(x in common for x in ("TT", "EE", "TE")):
            fields += ["delta_TE_over_sqrt_Kp_TT_EE"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for ell in ells:
            a, b = by_l_k[ell], by_l_v[ell]
            row = {"ell": ell}
            for x in common:
                va, vb = a[ik[x]], b[iv[x]]
                d = vb - va
                row[f"Kp_{x}"] = va
                row[f"NKp_v2_{x}"] = vb
                row[f"delta_{x}"] = d
                row[f"frac_{x}"] = d / va if va != 0.0 else ""
            if all(x in common for x in ("TT", "EE", "TE")):
                denom = math.sqrt(abs(a[ik["TT"]] * a[ik["EE"]]))
                row["delta_TE_over_sqrt_Kp_TT_EE"] = (
                    (b[iv["TE"]] - a[ik["TE"]]) / denom if denom > 0 else ""
                )
            w.writerow(row)

    bands = [(2, 29), (30, 499), (500, 1499), (1500, 3000)]
    print(f"Kp spectrum:     {kp_path}")
    print(f"NKp-v2 spectrum: {v2_path}")
    print("Coordinate-control differences (NKp-v2 minus Kp):")
    for lo, hi in bands:
        subset = [l for l in ells if lo <= l <= hi]
        if not subset:
            continue
        parts = []
        for x in ("TT", "EE"):
            if x not in common:
                continue
            vals = []
            for ell in subset:
                a, b = by_l_k[ell], by_l_v[ell]
                va = a[ik[x]]
                if va != 0.0:
                    vals.append(abs((b[iv[x]] - va) / va))
            if vals:
                parts.append(f"max |d{x}/{x}|={100*max(vals):.3f}%")
        print(f"  ell {lo:4d}-{hi:4d}: " + ", ".join(parts))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
