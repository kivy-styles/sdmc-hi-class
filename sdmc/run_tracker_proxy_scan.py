#!/usr/bin/env python3
"""Scan the HyRec-aware early-tracker background proxy.

This is a calibration diagnostic, not the final SDMC tracker implementation.
It varies only the Doran-Robbers early fraction Omega_e used as a proxy and
records the thermodynamic milestones returned by hi_class.  The purpose is to
ask whether one constant early fraction can simultaneously reproduce the
frozen SDMC recombination ruler r_s(z*) and the raw drag ruler r_d.
"""
from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
import tempfile

GRID = [0.0, 0.004, 0.006, 0.007, 0.0075, 0.008, 0.0085,
        0.009, 0.0095, 0.010, 0.0105, 0.011, 0.012]

BRANCHES = {
    "kp": {
        "template": "sdmc/config/kp_tracker_background_proxy.ini",
        "target_rs": 143.779000,
        "target_rd": 146.518000,
    },
    "nkp_v2": {
        "template": "sdmc/config/nkp_v2_tracker_background_proxy.ini",
        "target_rs": 142.688139,
        "target_rd": 145.395978,
    },
}

PATTERNS = {
    "age_Gyr": r"age = ([0-9.eE+-]+) Gyr",
    "z_eq": r"radiation/matter equality at z = ([0-9.eE+-]+)",
    "z_rec": r"recombination \(maximum of visibility function\) at z = ([0-9.eE+-]+)",
    "rs_rec_Mpc": r"with comoving sound horizon = ([0-9.eE+-]+) Mpc",
    "theta_s_100": r"sound horizon angle 100\*theta_s = ([0-9.eE+-]+)",
    "z_star": r"Thomson optical depth crosses one at z_\* = ([0-9.eE+-]+)",
    "theta_star_100": r"giving an angle 100\*theta_\* = ([0-9.eE+-]+)",
    "z_drag": r"baryon drag stops at z = ([0-9.eE+-]+)",
    "rd_Mpc": r"with comoving sound horizon rs = ([0-9.eE+-]+) Mpc",
}


def parse_value(text: str, pattern: str) -> float:
    m = re.search(pattern, text)
    if not m:
        raise RuntimeError(f"Could not parse pattern {pattern!r}")
    return float(m.group(1))


def run_one(class_exe: str, branch: str, omega_e: float) -> dict[str, float | str]:
    cfg = BRANCHES[branch]
    with open(cfg["template"], "r", encoding="utf-8") as f:
        text = f.read()

    text, n = re.subn(
        r"(?m)^expansion_smg\s*=\s*0\.70\s*,\s*-1\.0\s*,\s*[0-9.eE+-]+\s*$",
        f"expansion_smg = 0.70, -1.0, {omega_e:.8f}",
        text,
    )
    if n != 1:
        raise RuntimeError(f"Expected one expansion_smg line in {cfg['template']}, replaced {n}")

    # Keep scan products out of the permanent output namespace.  We only need stdout.
    text = re.sub(
        r"(?m)^root\s*=.*$",
        f"root = output/sdmc_tracker_scan_{branch}_{omega_e:.5f}_",
        text,
    )
    text = re.sub(r"(?m)^write background\s*=.*$", "write background = no", text)
    text = re.sub(r"(?m)^write thermodynamics\s*=.*$", "write thermodynamics = no", text)

    with tempfile.NamedTemporaryFile("w", suffix=".ini", delete=False, encoding="utf-8") as tf:
        tf.write(text)
        path = tf.name

    try:
        proc = subprocess.run([class_exe, path], text=True, capture_output=True)
    finally:
        os.unlink(path)

    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"hi_class failed for {branch} Omega_e={omega_e}")

    row: dict[str, float | str] = {"branch": branch, "Omega_e": omega_e}
    for key, pat in PATTERNS.items():
        row[key] = parse_value(proc.stdout, pat)
    row["target_rs_Mpc"] = cfg["target_rs"]
    row["target_rd_Mpc"] = cfg["target_rd"]
    row["delta_rs_Mpc"] = float(row["rs_rec_Mpc"]) - float(cfg["target_rs"])
    row["delta_rd_Mpc"] = float(row["rd_Mpc"]) - float(cfg["target_rd"])
    row["chi2_ruler_equal_weight"] = (
        float(row["delta_rs_Mpc"]) ** 2 + float(row["delta_rd_Mpc"]) ** 2
    )
    return row


def linear_root(rows: list[dict], key: str) -> float | None:
    ordered = sorted(rows, key=lambda r: float(r["Omega_e"]))
    for a, b in zip(ordered, ordered[1:]):
        ya, yb = float(a[key]), float(b[key])
        if ya == 0:
            return float(a["Omega_e"])
        if ya * yb <= 0 and yb != ya:
            xa, xb = float(a["Omega_e"]), float(b["Omega_e"])
            return xa + (0.0 - ya) * (xb - xa) / (yb - ya)
    return None


def main() -> None:
    class_exe = sys.argv[1] if len(sys.argv) > 1 else "./class"
    out = sys.argv[2] if len(sys.argv) > 2 else "output/sdmc_tracker_proxy_scan.csv"
    rows: list[dict] = []
    for branch in BRANCHES:
        for omega_e in GRID:
            row = run_one(class_exe, branch, omega_e)
            rows.append(row)
            print(
                f"{branch:7s} Omega_e={omega_e:0.5f} "
                f"rs={float(row['rs_rec_Mpc']):.6f} "
                f"rd={float(row['rd_Mpc']):.6f} "
                f"drs={float(row['delta_rs_Mpc']):+.6f} "
                f"drd={float(row['delta_rd_Mpc']):+.6f}"
            )

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fields = list(rows[0].keys())
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print("\nInterpolated constant-fraction proxy matches:")
    for branch in BRANCHES:
        rr = [r for r in rows if r["branch"] == branch]
        ers = linear_root(rr, "delta_rs_Mpc")
        erd = linear_root(rr, "delta_rd_Mpc")
        best = min(rr, key=lambda r: float(r["chi2_ruler_equal_weight"]))
        print(
            f"  {branch}: Omega_e(rs target)={ers if ers is not None else float('nan'):.6f}; "
            f"Omega_e(rd target)={erd if erd is not None else float('nan'):.6f}; "
            f"grid best={float(best['Omega_e']):.5f} "
            f"(drs={float(best['delta_rs_Mpc']):+.4f} Mpc, "
            f"drd={float(best['delta_rd_Mpc']):+.4f} Mpc)"
        )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
