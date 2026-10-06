#!/usr/bin/env python3
"""Scan scalar perturbation initial conditions for the NKp-v2 No-Slip branch.

All rows use the same frozen NKp-v2 background and EFT functions. Only
pert_initial_conditions_smg is changed. This is a solver/closure diagnostic,
not a likelihood fit. Its purpose is to determine whether the large late-time
response of the fully-dynamical No-Slip run is an initial-condition artifact
or a robust property of the implemented scalar mode.
"""
from __future__ import annotations

import csv
import glob
import math
import os
import subprocess
from pathlib import Path

TEMPLATE = Path("sdmc/config/nkp_v2_full_clustered_noslip.ini")
REF = Path("output/sdmc_nkp_v2_tracker_clustered_cs0003_00_cl_lensed.dat")
OUT = Path("output/sdmc_nkp_v2_noslip_ic_scan.csv")
ICS = ["ext_field_attr", "zero", "single_clock", "gravitating_attr", "kin_only"]


def replace_line(text: str, key: str, value: str) -> str:
    lines = text.splitlines()
    out = []
    n = 0
    for line in lines:
        if line.strip().startswith(key + " ="):
            out.append(f"{key} = {value}")
            n += 1
        else:
            out.append(line)
    if n != 1:
        raise RuntimeError(f"expected one {key}, found {n}")
    return "\n".join(out) + "\n"


def set_or_insert(text: str, key: str, value: str, before_key: str = "method_qs_smg") -> str:
    lines = text.splitlines()
    matches = [i for i, line in enumerate(lines) if line.strip().startswith(key + " =")]
    if len(matches) == 1:
        lines[matches[0]] = f"{key} = {value}"
    elif len(matches) == 0:
        pos = next((i for i, line in enumerate(lines)
                    if line.strip().startswith(before_key + " =")), None)
        if pos is None:
            raise RuntimeError(f"could not find insertion point {before_key}")
        lines.insert(pos, f"{key} = {value}")
    else:
        raise RuntimeError(f"expected at most one {key}, found {len(matches)}")
    return "\n".join(lines) + "\n"


def table(path: str | Path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                rows.append([float(x) for x in s.split()])
    return rows


def find1(pattern: str) -> str:
    x = sorted(glob.glob(pattern))
    if len(x) != 1:
        raise RuntimeError(f"{pattern}: found {len(x)}")
    return x[0]


def sigma8(rows) -> float:
    vals = []
    for r in rows:
        k, p = r[0], r[1]
        x = 8.0 * k
        w = 1.0 if abs(x) < 1e-8 else 3.0 * (math.sin(x) - x * math.cos(x)) / x**3
        vals.append((math.log(k), k**3 * p * w*w / (2.0 * math.pi**2)))
    q = sum(0.5 * (b[1] + a[1]) * (b[0] - a[0]) for a, b in zip(vals[:-1], vals[1:]))
    return math.sqrt(max(q, 0.0))


def nearest(rows, ell: int, col: int) -> float:
    return min(rows, key=lambda r: abs(r[0] - ell))[col]


def clean(prefix: str) -> None:
    for p in glob.glob(prefix + "*"):
        try:
            os.remove(p)
        except OSError:
            pass


def main() -> None:
    base = TEMPLATE.read_text(encoding="utf-8")
    ref = table(REF)
    fields = [
        "initial_condition", "status", "sigma8_z0", "sigma8_z10", "growth_0_over_10",
        "phiphi_ratio_l100", "phiphi_ratio_l500", "phiphi_ratio_l1000", "phiphi_ratio_l1500",
        "TT_ratio_l200", "TT_ratio_l1000", "TT_ratio_l2000", "error"
    ]
    rows = []

    for ic in ICS:
        prefix = f"output/sdmc_nkp_v2_noslip_ic_{ic}_"
        clean(prefix)
        ini = Path(f"output/sdmc_nkp_v2_noslip_ic_{ic}.ini")
        log = Path(f"output/sdmc_nkp_v2_noslip_ic_{ic}.log")

        text = set_or_insert(base, "pert_initial_conditions_smg", ic)
        text = replace_line(text, "method_qs_smg", "fully_dynamic")
        text = replace_line(text, "z_pk", "0,10")
        text = replace_line(text, "root", prefix)
        for k in ("input_verbose", "background_verbose", "thermodynamics_verbose", "output_verbose"):
            text = replace_line(text, k, "0")
        text = replace_line(text, "perturbations_verbose", "1")
        ini.write_text(text, encoding="utf-8")

        row = {k: "" for k in fields}
        row["initial_condition"] = ic
        try:
            cp = subprocess.run(
                ["./class", str(ini)], text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=240
            )
            log.write_text(cp.stdout, encoding="utf-8")
            if cp.returncode:
                row["status"] = "FAIL"
                tail = [x for x in cp.stdout.strip().splitlines() if x.strip()][-14:]
                row["error"] = " | ".join(tail)[:2400]
                print(ic, "FAIL", row["error"])
                rows.append(row)
                continue

            cl = table(find1(prefix + "*_cl_lensed.dat"))
            p0 = table(find1(prefix + "*_z1_pk.dat"))
            p10 = table(find1(prefix + "*_z2_pk.dat"))
            s0, s10 = sigma8(p0), sigma8(p10)
            row.update(
                status="OK",
                sigma8_z0=f"{s0:.12g}",
                sigma8_z10=f"{s10:.12g}",
                growth_0_over_10=f"{s0/s10:.12g}",
            )
            for ell in (100, 500, 1000, 1500):
                row[f"phiphi_ratio_l{ell}"] = f"{nearest(cl, ell, 5) / nearest(ref, ell, 5):.12g}"
            for ell in (200, 1000, 2000):
                row[f"TT_ratio_l{ell}"] = f"{nearest(cl, ell, 1) / nearest(ref, ell, 1):.12g}"
            print(ic, "OK", "sigma8", s0, "growth", s0/s10,
                  "phi1000/ref", row["phiphi_ratio_l1000"])
            rows.append(row)
        except subprocess.TimeoutExpired as exc:
            row["status"] = "TIMEOUT"
            row["error"] = f"timeout after {exc.timeout} s"
            print(ic, "TIMEOUT")
            rows.append(row)
        finally:
            clean(prefix)
            try:
                ini.unlink()
            except OSError:
                pass

    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
