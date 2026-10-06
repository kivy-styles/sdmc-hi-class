#!/usr/bin/env python3
"""Scan a general-braiding completion of the frozen NKp-v2 background.

NOTE: this file keeps its historical filename because it already existed in the
fork.  It is now a *runtime diagnostic*, not a source patch applied during the
initial build.

The endpoint-clean B0 No-Slip completion (alpha_B=-2 alpha_M) is stable at the
background level but produces excessive fully-dynamical growth and CMB lensing.
This scan asks a narrower question before changing the NKp-v2 background:

    can the same F(a), alpha_M(a), expansion history, tracker and kinetic
    normalization admit a healthier Horndeski braiding relation

        alpha_B = beta_B alpha_M ?

For every beta_B the script changes only the single braiding coefficient in the
already-patched gravity model, recompiles incrementally, retains native hi_class
stability rejection, and runs full-dynamic TT/TE/EE+lensing+P(k).  The kinetic
D(a) seed remains the B0 value, so this is a *completion search*, not a refit.
Any viable beta found here must subsequently be re-closed/re-shot with its own
native c_s^2 numerator before it can become a candidate model.

The source file and executable are restored to exact No-Slip beta_B=-2 at exit.
"""
from __future__ import annotations

import csv
import glob
import math
import os
import subprocess
from pathlib import Path

SOURCE = Path("gravity_smg/gravity_models_smg.c")
TEMPLATE = Path("sdmc/config/nkp_v2_full_clustered_noslip.ini")
REF = Path("output/sdmc_nkp_v2_tracker_clustered_cs0003_00_cl_lensed.dat")
OUT = Path("output/sdmc_nkp_v2_braid_scan.csv")

BETAS = [-4.0, -3.0, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25,
         -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
ANCHOR = "double ab = -2.*am;"

FIELDS = [
    "beta_B", "status", "returncode", "min_cs2", "min_D", "F0",
    "alphaM0", "alphaM_max", "z_alphaM_max", "beta1_peak",
    "sigma8_z0", "sigma8_z10", "growth_0_over_10",
    "phiphi_ratio_l100", "phiphi_ratio_l500", "phiphi_ratio_l1000",
    "phiphi_ratio_l1500", "TT_ratio_l200", "TT_ratio_l1000",
    "TT_ratio_l2000", "error"
]


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


def clean(prefix: str) -> None:
    for p in glob.glob(prefix + "*"):
        try:
            os.remove(p)
        except OSError:
            pass


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


def background_metrics(path: str | Path):
    b = table(path)
    # background.dat zero-based columns in this fork:
    # z=0, M*^2=22, alpha_M=27, c_s^2=29, D=30.
    z0 = min(b, key=lambda r: abs(r[0]))
    ammax = max(b, key=lambda r: r[27])
    return {
        "min_cs2": min(r[29] for r in b),
        "min_D": min(r[30] for r in b),
        "F0": z0[22],
        "alphaM0": z0[27],
        "alphaM_max": ammax[27],
        "z_alphaM_max": ammax[0],
    }


def rebuild() -> tuple[int, str]:
    cp = subprocess.run(
        ["make", "-j2", "class"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=240
    )
    return cp.returncode, cp.stdout


def main() -> None:
    original = SOURCE.read_text(encoding="utf-8")
    if original.count(ANCHOR) != 1:
        raise RuntimeError(f"expected exactly one beta=-2 braiding anchor, found {original.count(ANCHOR)}")

    base = TEMPLATE.read_text(encoding="utf-8")
    ref = table(REF)
    rows = []

    try:
        for beta in BETAS:
            row = {k: "" for k in FIELDS}
            row["beta_B"] = f"{beta:.12g}"
            tag = (f"{beta:+.2f}").replace("+", "p").replace("-", "m").replace(".", "p")
            prefix = f"output/sdmc_nkp_v2_braid_{tag}_"
            ini = Path(f"output/sdmc_nkp_v2_braid_{tag}.ini")
            log = Path(f"output/sdmc_nkp_v2_braid_{tag}.log")
            clean(prefix)

            replacement = f"double ab = ({beta:.17g})*am; /* general-braiding diagnostic */"
            SOURCE.write_text(original.replace(ANCHOR, replacement), encoding="utf-8")

            try:
                rc, buildlog = rebuild()
                if rc:
                    row["status"] = "BUILD_FAIL"
                    row["returncode"] = str(rc)
                    row["error"] = " | ".join(buildlog.strip().splitlines()[-20:])[:4000]
                    rows.append(row)
                    print(beta, "BUILD_FAIL")
                    continue

                text = replace_line(base, "method_qs_smg", "fully_dynamic")
                text = replace_line(text, "z_pk", "0,10")
                text = replace_line(text, "root", prefix)
                for key in ("input_verbose", "background_verbose", "thermodynamics_verbose", "perturbations_verbose", "output_verbose"):
                    text = replace_line(text, key, "0")
                ini.write_text(text, encoding="utf-8")

                cp = subprocess.run(
                    ["./class", str(ini)], text=True, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, timeout=300
                )
                log.write_text(cp.stdout, encoding="utf-8")
                row["returncode"] = str(cp.returncode)
                if cp.returncode:
                    row["status"] = "FAIL"
                    row["error"] = " | ".join(cp.stdout.strip().splitlines()[-24:])[:4000]
                    rows.append(row)
                    print(beta, "FAIL", row["error"][:300])
                    continue

                cl = table(find1(prefix + "*_cl_lensed.dat"))
                p0 = table(find1(prefix + "*_z1_pk.dat"))
                p10 = table(find1(prefix + "*_z2_pk.dat"))
                s0, s10 = sigma8(p0), sigma8(p10)
                m = background_metrics(find1(prefix + "*_background.dat"))

                row.update(
                    status="OK",
                    sigma8_z0=f"{s0:.12g}",
                    sigma8_z10=f"{s10:.12g}",
                    growth_0_over_10=f"{s0/s10:.12g}",
                    beta1_peak=f"{(2.0+beta)*m['alphaM_max']:.12g}",
                    **{k: f"{v:.12g}" for k, v in m.items()},
                )
                for ell in (100, 500, 1000, 1500):
                    row[f"phiphi_ratio_l{ell}"] = f"{nearest(cl, ell, 5)/nearest(ref, ell, 5):.12g}"
                for ell in (200, 1000, 2000):
                    row[f"TT_ratio_l{ell}"] = f"{nearest(cl, ell, 1)/nearest(ref, ell, 1):.12g}"

                print(beta, "OK", "min_cs2", row["min_cs2"], "sigma8", row["sigma8_z0"],
                      "phi1000/ref", row["phiphi_ratio_l1000"])
                rows.append(row)

            except subprocess.TimeoutExpired as exc:
                row["status"] = "TIMEOUT"
                row["error"] = f"timeout after {exc.timeout} s"
                rows.append(row)
                print(beta, "TIMEOUT")
            finally:
                clean(prefix)
                try:
                    ini.unlink()
                except OSError:
                    pass

    finally:
        # Leave both source and executable in the exact beta=-2 state expected
        # by all subsequent/future frozen-B0 diagnostics.
        SOURCE.write_text(original, encoding="utf-8")
        rc, buildlog = rebuild()
        if rc:
            raise RuntimeError("failed to restore beta=-2 executable: " + " | ".join(buildlog.strip().splitlines()[-20:]))

    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
