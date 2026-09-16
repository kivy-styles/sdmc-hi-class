#!/usr/bin/env python3
"""Audit the designer clustered SDMC tracker background tables.

Checks that the hi_class scalar sound speed follows the requested 0.003 target
and reports the kinetic coefficient D across the acoustic era and the z~20
handoff.  This is a numerical closure/stability diagnostic, not a likelihood.
"""
from __future__ import annotations

import glob
import math
import os
import re


def find_background(prefix):
    cands = sorted(glob.glob(prefix + "*_background.dat")) + sorted(glob.glob(prefix + "background.dat"))
    if not cands:
        raise FileNotFoundError(prefix)
    return cands[0]


def parse_header(line):
    # CLASS normally separates numbered titles with tabs.  Keep a regex
    # fallback for titles containing spaces such as "kin (D)".
    s = line.lstrip("#").strip()
    parts = [p.strip() for p in s.split("\t") if p.strip()]
    out = {}
    for p in parts:
        m = re.match(r"^(\d+)\s*:\s*(.*)$", p)
        if m:
            out[m.group(2).strip()] = int(m.group(1)) - 1
    return out


def read_background(path):
    cols = {}
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                h = parse_header(s)
                if h:
                    cols.update(h)
                continue
            rows.append([float(x) for x in s.split()])
    if not rows:
        raise RuntimeError(f"No numerical rows in {path}")
    return cols, rows


def find_col(cols, exact=None, contains=None):
    if exact:
        for name, idx in cols.items():
            if name == exact:
                return idx, name
    if contains:
        for name, idx in cols.items():
            if contains.lower() in name.lower():
                return idx, name
    raise KeyError(f"Could not find column exact={exact!r} contains={contains!r}; available={list(cols)}")


def audit(label, prefix, targets):
    path = find_background(prefix)
    cols, rows = read_background(path)
    iz, zn = find_col(cols, exact="z", contains="z")
    ics, csn = find_col(cols, exact="c_s^2", contains="c_s^2")
    iD, Dn = find_col(cols, exact="kin (D)", contains="kin (D)")
    ik, kn = find_col(cols, exact="kineticity_smg", contains="kineticity")

    finite = [r for r in rows if math.isfinite(r[ics]) and math.isfinite(r[iD])]
    mincs = min(finite, key=lambda r: r[ics])
    minD = min(finite, key=lambda r: r[iD])
    maxD = max(finite, key=lambda r: r[iD])

    print(f"\n{label}: {path}")
    print(f"  columns: z={zn!r}, cs2={csn!r}, D={Dn!r}, alphaK={kn!r}")
    print(f"  global min c_s^2 = {mincs[ics]:.8g} at z={mincs[iz]:.6g}")
    print(f"  global min D     = {minD[iD]:.8g} at z={minD[iz]:.6g}")
    print(f"  global max D     = {maxD[iD]:.8g} at z={maxD[iz]:.6g}")
    print("  selected epochs:")
    for zt in targets:
        r = min(finite, key=lambda q: abs(math.log1p(max(q[iz],0.0))-math.log1p(zt)))
        print(
            f"    z~{zt:10.3g}: actual z={r[iz]:11.5g}, "
            f"c_s^2={r[ics]:.8g}, D={r[iD]:.8g}, alpha_K={r[ik]:.8g}"
        )


def main():
    targets = [1.0e6, 3400.0, 1089.0, 100.0, 30.0, 21.0, 10.0, 1.0, 0.0]
    audit("Kp", "output/sdmc_kp_tracker_clustered_cs0003_", targets)
    audit("NKp-v2", "output/sdmc_nkp_v2_tracker_clustered_cs0003_", targets)


if __name__ == "__main__":
    main()
