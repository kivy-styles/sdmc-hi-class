#!/usr/bin/env python3
"""Evaluate a common Planck 2018 baseline-component screen on CLASS spectra.

Components:
  * high-l TTTEEE Plik-lite native;
  * low-l TT native;
  * low-l EE native;
  * Planck 2018 lensing native.

A single A_planck is shared by all components and is profiled jointly with the
standard Gaussian calibration prior. This is more complete than the project's
high-l-only Plik-lite statistic, but is still a baseline-component screen:
high-l remains nuisance-marginalized Plik-lite, not the full multifrequency
Plik nuisance likelihood.
"""
from __future__ import annotations

import argparse, csv, math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT as PlanckLowTT
from cobaya.likelihoods.planck_2018_lowl.EE import EE as PlanckLowEE
from cobaya.likelihoods.planck_2018_lensing import native as PlanckLensing

TCMB=2.7255
CAL_SIGMA=0.0025

def load_class(path: Path):
    a=np.loadtxt(path)
    ell=a[:,0].astype(int)
    lmax=int(ell.max())
    conv=(TCMB*1e6)**2
    tt=np.zeros(lmax+1); ee=np.zeros(lmax+1); te=np.zeros(lmax+1); pp=np.zeros(lmax+1)
    tt[ell]=a[:,1]*conv
    ee[ell]=a[:,2]*conv
    te[ell]=a[:,3]*conv
    # In the class-format cl_lensed.dat file used by this repository, the
    # written phiphi column is one factor L(L+1) below Cobaya's
    # provider.get_Cl(ell_factor=True)["pp"] convention. Multiplying by
    # L(L+1) reproduces the ~1e-7 Planck lensing bandpower scale.
    L=ell.astype(float)
    pp[ell]=a[:,5]*L*(L+1.0)
    high=np.column_stack([ell,tt[ell],te[ell],ee[ell]])
    return high, {"tt":tt,"te":te,"ee":ee,"pp":pp}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--packages",required=True)
    ap.add_argument("--output",required=True)
    ap.add_argument("spectra",nargs="+",help="label=CLASS_cl_lensed.dat")
    args=ap.parse_args()

    high=TTTEEE_lite_native(packages_path=args.packages)
    lowT=PlanckLowTT(packages_path=args.packages)
    lowE=PlanckLowEE(packages_path=args.packages)
    lens=PlanckLensing(packages_path=args.packages)

    rows=[]
    for item in args.spectra:
        label,raw=item.split("=",1)
        high_arr,dls=load_class(Path(raw))

        def components(A):
            A=float(A)
            h=float(high.chi_squared(high_arr,A_planck=A))
            lt=float(-2.0*lowT.log_likelihood(dls["tt"],calib=A))
            le=float(-2.0*lowE.log_likelihood(dls["ee"],calib=A))
            lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
            ll=float(-2.0*lens.log_likelihood(dls,**lp))
            pr=((A-1.0)/CAL_SIGMA)**2
            return h,lt,le,ll,pr,h+lt+le+ll+pr

        c1=components(1.0)
        opt=minimize_scalar(lambda A: components(A)[-1],bounds=(0.97,1.03),
                            method="bounded",options={"xatol":1e-10})
        A=float(opt.x)
        h,lt,le,ll,pr,total=components(A)
        rows.append({
            "case":label,"file":raw,
            "A_planck":A,
            "chi2_highl_pliklite":h,
            "chi2_lowl_TT":lt,
            "chi2_lowl_EE":le,
            "chi2_lensing":ll,
            "chi2_cal_prior":pr,
            "chi2_components_total":total,
            "chi2_components_A1":c1[-1],
            "lensing_calibration_param":getattr(lens,"calibration_param",None) or "",
        })

    out=Path(args.output); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    print("Planck 2018 baseline-component native screen")
    print("high-l Plik-lite + low-T + low-E + lensing; one shared A_planck")
    print("NOTE: this is not the full multifrequency high-l Plik nuisance likelihood.")
    for r in rows:
        print(
            f"{r['case']:24s} total={r['chi2_components_total']:.6f} "
            f"A={r['A_planck']:.8f} highl={r['chi2_highl_pliklite']:.6f} "
            f"lowT={r['chi2_lowl_TT']:.6f} lowE={r['chi2_lowl_EE']:.6f} "
            f"lens={r['chi2_lensing']:.6f} cal={r['chi2_cal_prior']:.6f}"
        )
    print(out)

if __name__=="__main__":
    main()
