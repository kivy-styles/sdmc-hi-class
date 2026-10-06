#!/usr/bin/env python3
"""Scan the late NKp-v2 No-Slip scalar sound speed with the full Boltzmann solver.

This is a perturbation-closure diagnostic, not a likelihood.  The frozen
NKp-v2 background, smooth endpoint-clean F/alpha_M profile, alpha_B=-2alpha_M,
and primordial parameters are held fixed.  Only the late target c_s^2 is
varied.  This tests whether the large late scalar response of the benchmark
c_s^2=0.1 run is tied to the kinetic handoff rather than the background.

For each successful run the script records native stability minima, the
endpoint Planck mass, alpha_M extrema, sigma8 at z=0 and z=10 computed directly
from P(k), and a few CMB/lensing ratios relative to the already-computed
NKp-v2 clustered-tracker control.
"""
from __future__ import annotations

import csv
import glob
import math
import os
from pathlib import Path
import subprocess

TEMPLATE = Path("sdmc/config/nkp_v2_full_clustered_noslip.ini")
REFERENCE_CL = Path("output/sdmc_nkp_v2_tracker_clustered_cs0003_00_cl_lensed.dat")
OUTCSV = Path("output/sdmc_nkp_v2_noslip_cs2_scan.csv")
VALUES = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]


def replace_line(text: str, key: str, value: str) -> str:
    lines = text.splitlines()
    hit = 0
    out = []
    for line in lines:
        if line.strip().startswith(key + " ="):
            out.append(f"{key} = {value}")
            hit += 1
        else:
            out.append(line)
    if hit != 1:
        raise RuntimeError(f"expected exactly one {key!r} line, found {hit}")
    return "\n".join(out) + "\n"


def read_table(path: Path):
    rows = []
    header = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if ":" in s and "1:" in s:
                    header = s
                continue
            rows.append([float(x) for x in s.split()])
    if not rows:
        raise RuntimeError(f"no numerical rows in {path}")
    return header, rows


def parse_numbered_header(line: str):
    import re
    s = line.lstrip("#").strip()
    marks = list(re.finditer(r"(\d+)\s*:\s*", s))
    cols = {}
    for j, m in enumerate(marks):
        end = marks[j + 1].start() if j + 1 < len(marks) else len(s)
        cols[s[m.end():end].strip()] = int(m.group(1)) - 1
    return cols


def sigma8(pk_rows):
    vals = []
    for k, p, *_ in pk_rows:
        x = 8.0 * k
        if abs(x) < 1.0e-8:
            w = 1.0
        else:
            w = 3.0 * (math.sin(x) - x * math.cos(x)) / (x ** 3)
        vals.append((math.log(k), k ** 3 * p * w * w / (2.0 * math.pi ** 2)))
    integ = 0.0
    for (x0, y0), (x1, y1) in zip(vals[:-1], vals[1:]):
        integ += 0.5 * (y0 + y1) * (x1 - x0)
    return math.sqrt(max(integ, 0.0))


def nearest(rows, ell, col):
    r = min(rows, key=lambda q: abs(q[0] - ell))
    return r[col]


def find_one(pattern: str) -> Path:
    c = sorted(glob.glob(pattern))
    if len(c) != 1:
        raise RuntimeError(f"expected one file for {pattern}, found {len(c)}")
    return Path(c[0])


def cleanup(prefix: str):
    for p in glob.glob(prefix + "*"):
        try:
            os.remove(p)
        except OSError:
            pass


def main():
    template = TEMPLATE.read_text(encoding="utf-8")
    _, refcl = read_table(REFERENCE_CL)

    fields = [
        "cs2_late", "status", "min_cs2", "z_min_cs2", "min_D", "z_min_D",
        "F0", "Fmax", "z_Fmax", "alphaM0", "alphaM_max", "z_alphaM_max",
        "alphaM_min", "z_alphaM_min", "sigma8_z0", "sigma8_z10",
        "growth_sigma8_0_over_10", "phiphi_ratio_l100", "phiphi_ratio_l500",
        "phiphi_ratio_l1000", "phiphi_ratio_l1500", "TT_ratio_l200",
        "TT_ratio_l1000", "TT_ratio_l2000", "error"
    ]
    results = []

    for cs2 in VALUES:
        tag = (f"{cs2:.3g}").replace(".", "p")
        root = f"output/sdmc_nkp_v2_noslip_cs2_{tag}_"
        cleanup(root)
        ini = Path(f"output/sdmc_nkp_v2_noslip_cs2_{tag}.ini")

        text = template
        text = replace_line(
            text,
            "parameters_smg",
            f"0.001, 0.003, {cs2:.12g}, 21.45883112, 21.10, 0.5",
        )
        text = replace_line(text, "z_pk", "0,10")
        text = replace_line(text, "root", root)
        text = replace_line(text, "background_verbose", "0")
        text = replace_line(text, "thermodynamics_verbose", "0")
        text = replace_line(text, "perturbations_verbose", "0")
        text = replace_line(text, "input_verbose", "0")
        text = replace_line(text, "output_verbose", "0")
        ini.write_text(text, encoding="utf-8")

        row = {k: "" for k in fields}
        row["cs2_late"] = f"{cs2:.12g}"
        try:
            cp = subprocess.run(
                ["./class", str(ini)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=180,
            )
            if cp.returncode != 0:
                row["status"] = "FAIL"
                tail = " | ".join(cp.stdout.strip().splitlines()[-5:])
                row["error"] = tail[:1000]
                results.append(row)
                print(f"c_s2_late={cs2:g}: FAIL {tail}")
                continue

            bg = find_one(root + "*_background.dat")
            cl = find_one(root + "*_cl_lensed.dat")
            pk0 = find_one(root + "*_z1_pk.dat")
            pk10 = find_one(root + "*_z2_pk.dat")

            hbg, bgrows = read_table(bg)
            cols = parse_numbered_header(hbg)
            iz = cols["z"]
            iF = cols["M*^2_smg"]
            iam = cols["M2_running_smg"]
            ics = cols["c_s^2"]
            iD = cols["kin (D)"]

            z0row = min(bgrows, key=lambda q: abs(q[iz]))
            mincs = min(bgrows, key=lambda q: q[ics])
            minD = min(bgrows, key=lambda q: q[iD])
            maxF = max(bgrows, key=lambda q: q[iF])
            maxam = max(bgrows, key=lambda q: q[iam])
            minam = min(bgrows, key=lambda q: q[iam])

            _, clrows = read_table(cl)
            _, p0 = read_table(pk0)
            _, p10 = read_table(pk10)
            s0 = sigma8(p0)
            s10 = sigma8(p10)

            row.update({
                "status": "OK",
                "min_cs2": f"{mincs[ics]:.12g}",
                "z_min_cs2": f"{mincs[iz]:.12g}",
                "min_D": f"{minD[iD]:.12g}",
                "z_min_D": f"{minD[iz]:.12g}",
                "F0": f"{z0row[iF]:.12g}",
                "Fmax": f"{maxF[iF]:.12g}",
                "z_Fmax": f"{maxF[iz]:.12g}",
                "alphaM0": f"{z0row[iam]:.12g}",
                "alphaM_max": f"{maxam[iam]:.12g}",
                "z_alphaM_max": f"{maxam[iz]:.12g}",
                "alphaM_min": f"{minam[iam]:.12g}",
                "z_alphaM_min": f"{minam[iz]:.12g}",
                "sigma8_z0": f"{s0:.12g}",
                "sigma8_z10": f"{s10:.12g}",
                "growth_sigma8_0_over_10": f"{s0/s10:.12g}",
            })
            for e in (100, 500, 1000, 1500):
                row[f"phiphi_ratio_l{e}"] = f"{nearest(clrows,e,5)/nearest(refcl,e,5):.12g}"
            for e in (200, 1000, 2000):
                row[f"TT_ratio_l{e}"] = f"{nearest(clrows,e,1)/nearest(refcl,e,1):.12g}"

            results.append(row)
            print(
                f"c_s2_late={cs2:g}: OK min_cs2={mincs[ics]:.6g} "
                f"minD={minD[iD]:.6g} sigma8={s0:.6g} "
                f"growth10to0={s0/s10:.6g} phi1000/ref={float(row['phiphi_ratio_l1000']):.6g}"
            )
        finally:
            # Keep only the compact scan table.  The production benchmark run
            # already preserves its complete spectra/background separately.
            cleanup(root)
            try:
                ini.unlink()
            except OSError:
                pass

    OUTCSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTCSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)
    print(f"Wrote {OUTCSV}")


if __name__ == "__main__":
    main()
