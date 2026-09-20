#!/usr/bin/env python3
"""
Matched linear / band-limited 3-D density-field comparison between flat LCDM
and the exact-leader SDMC r032 realization.

The two fields use:
  * identical primordial Gaussian random phases,
  * identical box/grid realization,
  * identical H0, omega_b, omega_cdm, A_s, n_s and tau_reio inputs,
  * model-specific linear matter power spectra P(k,z).

This is NOT a nonlinear N-body calculation.  The maps are deliberately
band-limited to k <= k_cut and should be interpreted as a controlled linear /
mildly-nonlinear structure comparison.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


ZPK = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]


def find_one(root: Path, basename: str) -> Path:
    matches = list(root.rglob(basename))
    if len(matches) != 1:
        raise FileNotFoundError(f"expected one {basename!r} below {root}, found {len(matches)}")
    return matches[0]


def load_pk(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    a = np.loadtxt(path)
    if a.ndim != 2 or a.shape[1] < 2:
        raise ValueError(f"bad P(k) table: {path}")
    k = np.asarray(a[:, 0], float)
    p = np.asarray(a[:, 1], float)
    q = np.isfinite(k) & np.isfinite(p) & (k > 0.0) & (p > 0.0)
    k, p = k[q], p[q]
    order = np.argsort(k)
    return k[order], p[order]


def log_interp(k: np.ndarray, kt: np.ndarray, pt: np.ndarray) -> np.ndarray:
    lk = np.log(np.clip(k, kt[0], kt[-1]))
    lp = np.interp(lk, np.log(kt), np.log(pt))
    out = np.exp(lp)
    out[(k < kt[0]) | (k > kt[-1])] = 0.0
    return out


def sigma_r_from_pk(k: np.ndarray, p: np.ndarray, R: float = 8.0) -> float:
    q = (k > 1e-4) & np.isfinite(p) & (p > 0)
    k = k[q]
    p = p[q]
    x = k * R
    W = np.ones_like(x)
    nz = x != 0
    W[nz] = 3.0 * (np.sin(x[nz]) - x[nz] * np.cos(x[nz])) / x[nz] ** 3
    y = k**3 * p * W**2 / (2.0 * np.pi**2)
    return float(np.sqrt(np.trapz(y, x=np.log(k))))


def make_kgrid(n: int, box: float):
    dx = box / n
    kx = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    kz = 2.0 * np.pi * np.fft.rfftfreq(n, d=dx)
    kk = np.sqrt(kx[:, None, None] ** 2 + ky[None, :, None] ** 2 + kz[None, None, :] ** 2)
    return kk, dx


def field_from_pk(
    white_k: np.ndarray,
    kk: np.ndarray,
    dx: float,
    ktab: np.ndarray,
    ptab: np.ndarray,
    kcut: float,
    n: int,
) -> np.ndarray:
    pk = log_interp(kk, ktab, ptab)
    # Smooth eighth-order low-pass keeps the visual field within the regime
    # where the linear Boltzmann spectra are the intended description.
    filt = np.exp(-((kk / kcut) ** 8))
    filt[kk == 0] = 0.0
    cell_volume = dx**3
    amp = np.sqrt(pk / cell_volume) * filt
    dk = white_k * amp
    return np.fft.irfftn(dk, s=(n, n, n)).real


def shell_power(delta: np.ndarray, box: float, nbins: int = 28):
    n = delta.shape[0]
    dx = box / n
    vol = box**3
    dk = np.fft.rfftn(delta)
    kx = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    kz = 2.0 * np.pi * np.fft.rfftfreq(n, d=dx)
    kk = np.sqrt(kx[:, None, None] ** 2 + ky[None, :, None] ** 2 + kz[None, None, :] ** 2)
    pmode = (dx**6) * np.abs(dk) ** 2 / vol

    q = (kk > 0) & np.isfinite(pmode)
    kvals = kk[q]
    pvals = pmode[q]
    edges = np.geomspace(kvals.min(), kvals.max(), nbins + 1)
    ib = np.digitize(kvals, edges) - 1
    kout, pout, nm = [], [], []
    for i in range(nbins):
        m = ib == i
        if np.count_nonzero(m) < 8:
            continue
        kout.append(float(np.exp(np.mean(np.log(kvals[m])))))
        pout.append(float(np.mean(pvals[m])))
        nm.append(int(np.count_nonzero(m)))
    return np.asarray(kout), np.asarray(pout), np.asarray(nm)


def slab(delta: np.ndarray, thickness: int = 8) -> np.ndarray:
    n = delta.shape[2]
    c = n // 2
    h = max(1, thickness // 2)
    return np.mean(delta[:, :, c - h:c + h], axis=2)


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    aa = a.ravel() - np.mean(a)
    bb = b.ravel() - np.mean(b)
    den = np.sqrt(np.dot(aa, aa) * np.dot(bb, bb))
    return float(np.dot(aa, bb) / den) if den > 0 else float("nan")


def make_frame(
    out: Path,
    z: float,
    lcdm: np.ndarray,
    sdmc: np.ndarray,
    thickness: int,
) -> Path:
    a = slab(lcdm, thickness)
    b = slab(sdmc, thickness)
    sig = float(np.std(a))
    d = (b - a) / max(sig, 1e-30)

    vmax = float(np.quantile(np.abs(np.concatenate([a.ravel(), b.ravel()])), 0.995))
    dmax = float(np.quantile(np.abs(d.ravel()), 0.995))
    vmax = max(vmax, 1e-6)
    dmax = max(dmax, 1e-6)

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2), constrained_layout=True)
    im0 = ax[0].imshow(a.T, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax[0].set_title(r"Matched $\Lambda$CDM")
    im1 = ax[1].imshow(b.T, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax[1].set_title("SDMC r032 exact leader")
    im2 = ax[2].imshow(d.T, origin="lower", cmap="PuOr_r", vmin=-dmax, vmax=dmax)
    ax[2].set_title(r"$(\delta_{\rm SDMC}-\delta_{\Lambda CDM})/\sigma_{\Lambda CDM}$")
    for x in ax:
        x.set_xticks([])
        x.set_yticks([])
    fig.colorbar(im1, ax=ax[:2], shrink=0.78, label=r"projected linear $\delta$")
    fig.colorbar(im2, ax=ax[2], shrink=0.78, label="normalized residual")
    fig.suptitle(f"Same phases, band-limited linear density field — z={z:g}")
    path = out / f"density_compare_z{str(z).replace('.', 'p')}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sdmc-root", required=True, type=Path)
    ap.add_argument("--lcdm-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--grid", type=int, default=128)
    ap.add_argument("--box", type=float, default=500.0, help="box side in Mpc/h")
    ap.add_argument("--kcut", type=float, default=0.25, help="linear low-pass scale in h/Mpc")
    ap.add_argument("--seed", type=int, default=20260920)
    ap.add_argument("--slab", type=int, default=8)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    n = int(args.grid)
    if n < 32 or n > 256:
        raise ValueError("grid must be between 32 and 256")
    if args.box <= 0 or args.kcut <= 0:
        raise ValueError("box and kcut must be positive")

    sdmc_files = [find_one(args.sdmc_root, f"covariant_observable_00_z{i}_pk.dat") for i in range(1, 8)]
    lcdm_files = [find_one(args.lcdm_root, f"matched_lcdm_00_z{i}_pk.dat") for i in range(1, 8)]

    rng = np.random.default_rng(args.seed)
    white = rng.normal(size=(n, n, n))
    white -= white.mean()
    white_k = np.fft.rfftn(white)
    kk, dx = make_kgrid(n, args.box)

    metrics = []
    target_ratios: Dict[float, Tuple[np.ndarray, np.ndarray]] = {}
    frame_paths = []

    for z, fs, fl in zip(ZPK, sdmc_files, lcdm_files):
        ks, ps = load_pk(fs)
        kl, pl = load_pk(fl)

        lcdm = field_from_pk(white_k, kk, dx, kl, pl, args.kcut, n)
        sdmc = field_from_pk(white_k, kk, dx, ks, ps, args.kcut, n)
        diff = sdmc - lcdm

        common_min = max(ks.min(), kl.min(), 0.01)
        common_max = min(ks.max(), kl.max(), args.kcut)
        kval = np.geomspace(common_min, common_max, 300)
        psl = log_interp(kval, ks, ps)
        pll = log_interp(kval, kl, pl)
        ratio = psl / np.maximum(pll, 1e-300)
        target_ratios[z] = (kval, ratio)

        row = {
            "z": z,
            "sigma8_target_lcdm": sigma_r_from_pk(kl, pl),
            "sigma8_target_sdmc": sigma_r_from_pk(ks, ps),
            "sigma8_ratio_sdmc_over_lcdm": sigma_r_from_pk(ks, ps) / sigma_r_from_pk(kl, pl),
            "field_rms_lcdm": float(np.std(lcdm)),
            "field_rms_sdmc": float(np.std(sdmc)),
            "field_rms_diff": float(np.std(diff)),
            "field_corrcoef": safe_corr(lcdm, sdmc),
            "mean_pk_ratio_0p02_0p20": float(np.mean(ratio[(kval >= 0.02) & (kval <= min(0.20, common_max))])),
        }
        for kp in (0.02, 0.05, 0.10, 0.20):
            if kp <= common_max:
                row[f"pk_ratio_k{str(kp).replace('.', 'p')}"] = float(np.interp(kp, kval, ratio))
        metrics.append(row)

        kgl, pgl, nml = shell_power(lcdm, args.box)
        kgs, pgs, nms = shell_power(sdmc, args.box)
        np.savetxt(
            args.out / f"realized_power_z{str(z).replace('.', 'p')}.csv",
            np.column_stack([kgl, pgl, np.interp(kgl, kgs, pgs)]),
            delimiter=",",
            header="k_hMpc,P_realized_LCDM,P_realized_SDMC",
            comments="",
        )

        frame_paths.append(make_frame(args.out, z, lcdm, sdmc, args.slab))

        # Free the two real-space cubes before the next redshift.
        del lcdm, sdmc, diff

    fields = sorted({k for row in metrics for k in row.keys()}, key=lambda x: (x != "z", x))
    with (args.out / "matched_simulation_metrics.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(metrics)

    # Power-spectrum ratio figure.
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    for z in ZPK:
        k, r = target_ratios[z]
        ax.semilogx(k, r, label=f"z={z:g}")
    ax.axhline(1.0, lw=1.0, color="black", alpha=0.6)
    ax.set_xlim(0.01, args.kcut)
    ax.set_xlabel(r"$k\;[h\,{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$P_{\rm SDMC}/P_{\Lambda CDM}$")
    ax.set_title("Exact-leader SDMC / matched ΛCDM linear power ratio")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out / "pk_ratio_evolution.png", dpi=180)
    plt.close(fig)

    # sigma8 history from the actual target spectra.
    zs = np.asarray([m["z"] for m in metrics])
    s_l = np.asarray([m["sigma8_target_lcdm"] for m in metrics])
    s_s = np.asarray([m["sigma8_target_sdmc"] for m in metrics])
    order = np.argsort(zs)
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.plot(zs[order], s_l[order], marker="o", label="matched ΛCDM")
    ax.plot(zs[order], s_s[order], marker="o", label="SDMC")
    ax.set_xlabel("redshift z")
    ax.set_ylabel(r"$\sigma_8(z)$ from linear $P(k,z)$")
    ax.set_title("Matched growth-amplitude comparison")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out / "sigma8_evolution.png", dpi=180)
    plt.close(fig)

    # Animated comparison from high z to today.
    ordered_frames = [frame_paths[ZPK.index(z)] for z in sorted(ZPK, reverse=True)]
    imgs = [Image.open(p).convert("P", palette=Image.Palette.ADAPTIVE) for p in ordered_frames]
    imgs[0].save(
        args.out / "matched_density_evolution.gif",
        save_all=True,
        append_images=imgs[1:],
        duration=850,
        loop=0,
        optimize=False,
    )
    for im in imgs:
        im.close()

    summary = [
        "Matched LCDM–SDMC linear 3-D simulation",
        "==========================================",
        f"seed = {args.seed}",
        f"grid = {n}^3",
        f"box = {args.box:g} Mpc/h",
        f"k_cut = {args.kcut:g} h/Mpc",
        "same Gaussian phases used for both cosmologies at every redshift",
        "interpretation = linear/band-limited; not a nonlinear N-body calculation",
        "",
        "z    sigma8_LCDM   sigma8_SDMC   ratio      corr(field)   rms(diff)",
    ]
    for m in sorted(metrics, key=lambda r: r["z"]):
        summary.append(
            f"{m['z']:4.1f}  {m['sigma8_target_lcdm']:.8f}  "
            f"{m['sigma8_target_sdmc']:.8f}  {m['sigma8_ratio_sdmc_over_lcdm']:.8f}  "
            f"{m['field_corrcoef']:.9f}  {m['field_rms_diff']:.8e}"
        )
    (args.out / "simulation_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary), flush=True)


if __name__ == "__main__":
    main()
