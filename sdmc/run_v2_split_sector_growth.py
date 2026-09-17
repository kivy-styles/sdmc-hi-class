#!/usr/bin/env python3
"""Phenomenological split-sector growth diagnostic for NKp-v2.

This is deliberately NOT a covariant one-field completion.  It answers a
narrow question raised by the full hi_class run: what happens if the validated
early clustered tracker fixes the primordial/early normalization, while the
late structural release acts only through the exact No-Slip response

    mu(a) = Sigma(a) = 1/F(a),   eta = 1,

without inheriting the propagating early Horndeski scalar perturbation?

The script reads the branch-complete NKp-v2 background produced by hi_class,
uses its tabulated F=M_*^2/M_Pl^2 and H(a), and solves the standard
scale-independent linear growth equation from z=20 to z=0,

    D_NN + [2 + d ln H/dN] D_N - 3 Omega_m mu D / 2 = 0,

for both mu=1 and mu=1/F.  The two runs use identical initial conditions, so
all quoted ratios isolate the late No-Slip response on the SAME background.

For a useful absolute normalization screen only, sigma8(z=10) is read from the
full Boltzmann run.  Since the late response is already weak at z=10 this is a
controlled diagnostic, not a replacement for a true split-sector Boltzmann
implementation.
"""
from __future__ import annotations

from pathlib import Path
import csv
import math
import numpy as np

BG = Path("output/sdmc_nkp_v2_full_noslip_00_background.dat")
PK10 = Path("output/sdmc_nkp_v2_full_noslip_00_z7_pk.dat")
OUT = Path("output/sdmc_nkp_v2_split_sector_growth.csv")


def sigma8_from_pk(path: Path) -> float:
    a = np.loadtxt(path)
    k = a[:, 0]
    pk = a[:, 1]
    x = 8.0 * k
    W = np.ones_like(x)
    m = np.abs(x) > 1.0e-10
    W[m] = 3.0 * (np.sin(x[m]) - x[m] * np.cos(x[m])) / x[m] ** 3
    integrand = k**3 * pk * W**2 / (2.0 * math.pi**2)
    return float(np.sqrt(np.trapezoid(integrand, x=np.log(k))))


def prepare_background(path: Path):
    a = np.loadtxt(path)
    # CLASS columns in the long background table:
    # 1:z, 4:H, 18:Omega_m(z), 23:M*^2_smg.
    z = a[:, 0]
    H = a[:, 3]
    Om = a[:, 17]
    F = a[:, 22]
    N = -np.log1p(z)
    order = np.argsort(N)
    N, H, Om, F = N[order], H[order], Om[order], F[order]

    N_start = -math.log(21.0)  # z=20, after the early tracker handoff.
    keep = N >= N_start
    N, H, Om, F = N[keep], H[keep], Om[keep], F[keep]
    h = np.gradient(np.log(H), N)
    return N, h, Om, F


def integrate_growth(Ntab, htab, Omtab, Ftab, noslip: bool):
    N0 = float(Ntab[0])
    nstep = 8000
    grid = np.linspace(N0, 0.0, nstep + 1)
    dN = float(grid[1] - grid[0])
    D = np.empty_like(grid)
    V = np.empty_like(grid)
    D[0] = 1.0
    V[0] = 1.0

    def rhs(N, d, v):
        h = float(np.interp(N, Ntab, htab))
        om = float(np.interp(N, Ntab, Omtab))
        F = float(np.interp(N, Ntab, Ftab))
        mu = 1.0 / F if noslip else 1.0
        return v, -(2.0 + h) * v + 1.5 * om * mu * d

    for i in range(nstep):
        x = grid[i]
        d, v = D[i], V[i]
        k1d, k1v = rhs(x, d, v)
        k2d, k2v = rhs(x + 0.5*dN, d + 0.5*dN*k1d, v + 0.5*dN*k1v)
        k3d, k3v = rhs(x + 0.5*dN, d + 0.5*dN*k2d, v + 0.5*dN*k2v)
        k4d, k4v = rhs(x + dN, d + dN*k3d, v + dN*k3v)
        D[i+1] = d + dN*(k1d + 2*k2d + 2*k3d + k4d)/6.0
        V[i+1] = v + dN*(k1v + 2*k2v + 2*k3v + k4v)/6.0
    return grid, D, V


def sample(grid, arr, z):
    N = -math.log1p(z)
    return float(np.interp(N, grid, arr))


def main():
    if not BG.exists():
        raise SystemExit(f"missing required background: {BG}")
    Ntab, htab, Omtab, Ftab = prepare_background(BG)
    Ng, Dg, Vg = integrate_growth(Ntab, htab, Omtab, Ftab, False)
    Ns, Ds, Vs = integrate_growth(Ntab, htab, Omtab, Ftab, True)

    sig10 = sigma8_from_pk(PK10) if PK10.exists() else float("nan")
    Dg10, Ds10 = sample(Ng, Dg, 10.0), sample(Ns, Ds, 10.0)
    split_sig0 = sig10 * Ds[-1] / Ds10 if math.isfinite(sig10) else float("nan")
    gr_sig0 = sig10 * Dg[-1] / Dg10 if math.isfinite(sig10) else float("nan")

    rows = []
    for z in (20.0, 10.0, 5.0, 3.0, 2.0, 1.0, 0.5, 0.0):
        dg, vg = sample(Ng, Dg, z), sample(Ng, Vg, z)
        ds, vs = sample(Ns, Ds, z), sample(Ns, Vs, z)
        N = -math.log1p(z)
        F = float(np.interp(N, Ntab, Ftab))
        rows.append({
            "z": z,
            "F": F,
            "mu_split": 1.0/F,
            "delta_D_percent": 100.0*(ds/dg - 1.0),
            "delta_fD_percent": 100.0*(vs/vg - 1.0),
            "f_GR": vg/dg,
            "f_split": vs/ds,
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    print("NKp-v2 split-sector late-growth diagnostic")
    print(f"  sigma8(z=10) Boltzmann anchor = {sig10:.9f}")
    print(f"  sigma8(0) same-background GR  = {gr_sig0:.9f}")
    print(f"  sigma8(0) split No-Slip       = {split_sig0:.9f}")
    print(f"  split delta D(0)              = {rows[-1]['delta_D_percent']:.6f}%")
    print(f"  split delta fD(0)             = {rows[-1]['delta_fD_percent']:.6f}%")
    print(f"  f_split(0)                    = {rows[-1]['f_split']:.9f}")
    print(f"  wrote {OUT}")


if __name__ == "__main__":
    main()
