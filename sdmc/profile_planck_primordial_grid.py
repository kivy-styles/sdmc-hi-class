#!/usr/bin/env python3
"""Profile A_s and n_s against Planck 2018 high-l TTTEEE Plik-lite.

This is a deterministic diagnostic grid, not a cosmological MCMC. It answers
one specific question: how much of the current SDMC high-l mismatch is
removable by the standard primordial amplitude and tilt while holding each
branch's background and early-sector parameters fixed.

The first pass used a broad common 5x5 grid.  This version keeps that broad
mode available, but by default uses branch-specific refined 5x5 grids centred
on the minima found in the first pass.  A_planck is profiled at every point
with the standard 0.0025 Gaussian calibration prior.  tau_reio remains fixed.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native

TCMB = 2.7255
CAL_SIGMA = 0.0025


def replace_scalar(text: str, key: str, value: str) -> str:
    pat = re.compile(rf"(?m)^\s*{re.escape(key)}\s*=.*$")
    if not pat.search(text):
        raise ValueError(f"missing {key} in base ini")
    return pat.sub(f"{key} = {value}", text)


def load_class_lensed(path: Path) -> np.ndarray:
    a = np.loadtxt(path)
    ell, tt, ee, te = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    conv = (TCMB * 1.0e6) ** 2
    return np.column_stack([ell, tt * conv, te * conv, ee * conv])


def profiled_chi2(lik, cl: np.ndarray):
    def like(A):
        return float(lik.chi_squared(cl, A_planck=float(A)))
    def total(A):
        return like(A) + ((float(A)-1.0)/CAL_SIGMA)**2
    opt = minimize_scalar(total, bounds=(0.96, 1.04), method="bounded",
                          options={"xatol": 1e-10})
    A = float(opt.x)
    c_like = like(A)
    c_prior = ((A-1.0)/CAL_SIGMA)**2
    return A, c_like, c_prior, c_like+c_prior


def run_case(lik, label: str, ini_path: Path, lnAs_grid, ns_grid, work: Path):
    base = ini_path.read_text()
    rows = []
    best = None
    counter = 0
    for lnAs in lnAs_grid:
        As = math.exp(float(lnAs)) / 1.0e10
        for ns in ns_grid:
            counter += 1
            stem = f"grid_{label}_{counter:03d}"
            root = work / (stem + "_")
            ini = base
            ini = replace_scalar(ini, "A_s", f"{As:.12e}")
            ini = replace_scalar(ini, "n_s", f"{float(ns):.9f}")
            ini = replace_scalar(ini, "root", str(root))
            if re.search(r"(?m)^\s*write parameters\s*=", ini):
                ini = replace_scalar(ini, "write parameters", "no")
            tmp_ini = work / f"{stem}.ini"
            tmp_ini.write_text(ini)
            proc = subprocess.run(["./class", str(tmp_ini)], stdout=subprocess.DEVNULL,
                                  stderr=subprocess.PIPE, text=True)
            if proc.returncode != 0:
                row = {"case": label, "ln10As": lnAs, "A_s": As, "n_s": ns,
                       "status": "FAIL", "A_planck": "", "chi2_like": "",
                       "chi2_cal_prior": "", "chi2_total": "",
                       "error": proc.stderr[-1000:]}
                rows.append(row)
                continue
            clfile = Path(str(root) + "00_cl_lensed.dat")
            cl = load_class_lensed(clfile)
            Acal, clike, cprior, ctot = profiled_chi2(lik, cl)
            row = {"case": label, "ln10As": float(lnAs), "A_s": As,
                   "n_s": float(ns), "status": "OK", "A_planck": Acal,
                   "chi2_like": clike, "chi2_cal_prior": cprior,
                   "chi2_total": ctot, "error": ""}
            rows.append(row)
            if best is None or ctot < best["chi2_total"]:
                best = dict(row)
                best["spectrum"] = str(clfile)
                best_copy = work / f"best_{label}_cl_lensed.dat"
                best_copy.write_bytes(clfile.read_bytes())
                best["best_spectrum"] = str(best_copy)
            for p in work.glob(stem + "_00_*"):
                try: p.unlink()
                except OSError: pass
            try: tmp_ini.unlink()
            except OSError: pass
            print(f"{label:24s} ln10As={lnAs:.4f} ns={ns:.4f} "
                  f"chi2={ctot:.3f} Acal={Acal:.6f}", flush=True)
    return rows, best


def grids_for(label: str, broad: bool):
    if broad:
        return (np.array([2.98, 3.02, 3.06, 3.10, 3.14]),
                np.array([0.94, 0.955, 0.970, 0.985, 1.000]))
    # Centres estimated from a local quadratic fit to the broad-grid results.
    if label == "LCDM_control":
        return (np.array([3.032, 3.038, 3.044, 3.050, 3.056]),
                np.array([0.955, 0.958, 0.961, 0.964, 0.967]))
    if label == "Kp_tracker_clustered":
        return (np.array([3.052, 3.058, 3.064, 3.070, 3.076]),
                np.array([0.970, 0.974, 0.978, 0.982, 0.986]))
    if label == "NKp_v2_tracker_clustered":
        return (np.array([3.066, 3.072, 3.078, 3.084, 3.090]),
                np.array([0.966, 0.970, 0.974, 0.978, 0.982]))
    raise ValueError(f"no refined grid defined for {label}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packages", required=True)
    ap.add_argument("--output", default="output/sdmc_planck_primordial_grid.csv")
    ap.add_argument("--best-output", default="output/sdmc_planck_primordial_best.csv")
    ap.add_argument("--broad", action="store_true",
                    help="use the original common broad 5x5 grid")
    ap.add_argument("cases", nargs="+", help="label=base.ini")
    args = ap.parse_args()

    lik = TTTEEE_lite_native(packages_path=args.packages)
    work = Path("output/planck_primordial_grid")
    work.mkdir(parents=True, exist_ok=True)
    allrows, bestrows = [], []
    for item in args.cases:
        label, raw = item.split("=", 1)
        lnAs_grid, ns_grid = grids_for(label, args.broad)
        rows, best = run_case(lik, label, Path(raw), lnAs_grid, ns_grid, work)
        allrows.extend(rows)
        if best is not None:
            bestrows.append(best)

    out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(allrows[0].keys()))
        w.writeheader(); w.writerows(allrows)
    bout = Path(args.best_output)
    with bout.open("w", newline="") as f:
        fields = list(bestrows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(bestrows)

    print("\nBEST GRID POINTS")
    for r in bestrows:
        print(f"{r['case']:24s} chi2={r['chi2_total']:.6f} "
              f"ln10As={r['ln10As']:.4f} n_s={r['n_s']:.4f} "
              f"A_planck={r['A_planck']:.7f}")
    print(out)
    print(bout)


if __name__ == "__main__":
    main()
