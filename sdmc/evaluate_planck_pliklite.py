#!/usr/bin/env python3
"""Evaluate Planck 2018 high-l TTTEEE Plik-lite on CLASS spectra.

This is a real nuisance-marginalized Planck *high-l* likelihood evaluation,
using Cobaya's native Python reimplementation of official plik_lite.  It is
not yet the complete Planck baseline because low-l TT/EE and lensing are not
included here.

CLASS files in this project use `format = class`, so the CMB columns are
DIMENSIONLESS D_l = l(l+1) C_l / (2 pi) ordered

    l, TT, EE, TE, BB, ...

Plik-lite's `chi_squared` helper expects

    l, TT, TE, EE

in D_l units of micro-K^2.  Hence temperature/polarization spectra are
multiplied by (T_CMB * 1e6)^2 before evaluation.

For every spectrum we report three values:
  * chi2_A1: likelihood chi2 at A_planck=1;
  * chi2_like_profile: likelihood-only profile over A_planck;
  * chi2_with_cal_prior: profile of likelihood + Planck calibration prior
        ((A_planck-1)/0.0025)^2.

The final quantity is useful for model comparison because it includes the
standard calibration prior explicitly, while keeping that prior separate in
the output.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import math
import numpy as np
from scipy.optimize import minimize_scalar

from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import (
    TTTEEE_lite_native,
)

TCMB = 2.7255
CAL_SIGMA = 0.0025


def load_class_lensed(path: Path) -> np.ndarray:
    a = np.loadtxt(path)
    if a.ndim != 2 or a.shape[1] < 4:
        raise ValueError(f"unexpected CLASS lensed Cl layout: {path} {a.shape}")
    ell = a[:, 0]
    tt = a[:, 1]
    ee = a[:, 2]
    te = a[:, 3]
    conv = (TCMB * 1.0e6) ** 2
    # Cobaya PlanckPlikLite.chi_squared expects L, D_TT, D_TE, D_EE.
    return np.column_stack([ell, tt * conv, te * conv, ee * conv])


def evaluate(lik, cl: np.ndarray):
    def like_chi2(A):
        return float(lik.chi_squared(cl, A_planck=float(A)))

    c1 = like_chi2(1.0)
    opt_like = minimize_scalar(like_chi2, bounds=(0.97, 1.03), method="bounded",
                               options={"xatol": 1e-10})

    def total(A):
        A = float(A)
        return like_chi2(A) + ((A - 1.0) / CAL_SIGMA) ** 2

    opt_prior = minimize_scalar(total, bounds=(0.97, 1.03), method="bounded",
                                options={"xatol": 1e-10})
    A = float(opt_prior.x)
    c_like = like_chi2(A)
    c_prior = ((A - 1.0) / CAL_SIGMA) ** 2
    return {
        "chi2_A1": c1,
        "A_like_profile": float(opt_like.x),
        "chi2_like_profile": float(opt_like.fun),
        "A_with_cal_prior": A,
        "chi2_like_at_cal_profile": c_like,
        "chi2_cal_prior": c_prior,
        "chi2_with_cal_prior": c_like + c_prior,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packages", required=True)
    ap.add_argument("--output", default="output/sdmc_planck_pliklite.csv")
    ap.add_argument("spectra", nargs="+", help="label=path to CLASS cl_lensed.dat")
    args = ap.parse_args()

    lik = TTTEEE_lite_native(packages_path=args.packages)
    rows = []
    for item in args.spectra:
        if "=" not in item:
            raise ValueError(f"expected label=path, got {item}")
        label, raw = item.split("=", 1)
        path = Path(raw)
        cl = load_class_lensed(path)
        r = {"case": label, "file": str(path)}
        r.update(evaluate(lik, cl))
        rows.append(r)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    print("Planck 2018 high-l TTTEEE Plik-lite native")
    print("NOTE: high-l only; low-l TT/EE and lensing not included.")
    for r in rows:
        print(
            f"{r['case']:24s}  chi2(A=1)={r['chi2_A1']:.6f}  "
            f"chi2(profile+calprior)={r['chi2_with_cal_prior']:.6f}  "
            f"A_planck={r['A_with_cal_prior']:.8f}"
        )
    print(out)


if __name__ == "__main__":
    main()
