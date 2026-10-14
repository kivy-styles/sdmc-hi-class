#!/usr/bin/env python3
"""Instrument hi_class and decompose the native scalar cs2 numerator.

This is a diagnostic only.  It exposes the exact internal non-SMG density and
pressure used by gravity_functions_smg.c, then reconstructs every term in the
native cs2 numerator for selected SDMC transition-bridge runs.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path

OUT = Path("output")
OMEGA_CDM = 0.12000047526014168
OMEGA_FLD = 3.0549311e-8
EXPANSION = "0.716246932, 20.0, 21.26, 0.5, 0.021223342, 0.25, 0.0116038358, 1.5"
CASES = [
    ("baseline_A0", "0.1, 1.0e-4, 0.0, 15.0, -1.0, 1.0, 1.0"),
    ("F_bump_bestwide", "0.1, 1.0e-4, 0.01, 15.0, -1.0, 1.4, 1.0"),
]


def insert_after_line(text: str, needle: str, additions: list[str]) -> str:
    pos = text.find(needle)
    if pos < 0:
        raise RuntimeError(f"anchor not found: {needle}")
    line_start = text.rfind("\n", 0, pos) + 1
    indent = text[line_start:pos]
    line_end = text.find("\n", pos)
    if line_end < 0:
        line_end = len(text)
        tail = "\n"
    else:
        tail = ""
    block = "".join(indent + line + "\n" for line in additions)
    return text[: line_end + 1] + block + text[line_end + 1 :] + tail


def instrument() -> None:
    p = Path("gravity_smg/background_smg.c")
    s = p.read_text(encoding="utf-8")
    title = 'class_store_columntitle(titles,"rho_tot_wo_smg_dbg",_TRUE_);'
    data = 'class_store_double(dataptr,pvecback[pba->index_bg_rho_tot_wo_smg],_TRUE_,storeidx);'
    if title not in s:
        s = insert_after_line(
            s,
            'class_store_columntitle(titles,"cs2num_p",_TRUE_);',
            [
                'class_store_columntitle(titles,"rho_tot_wo_smg_dbg",_TRUE_);',
                'class_store_columntitle(titles,"p_tot_wo_smg_dbg",_TRUE_);',
            ],
        )
    if data not in s:
        s = insert_after_line(
            s,
            'class_store_double(dataptr,pvecback[pba->index_bg_cs2num_prime_smg],_TRUE_,storeidx);',
            [
                'class_store_double(dataptr,pvecback[pba->index_bg_rho_tot_wo_smg],_TRUE_,storeidx);',
                'class_store_double(dataptr,pvecback[pba->index_bg_p_tot_wo_smg],_TRUE_,storeidx);',
            ],
        )
    p.write_text(s, encoding="utf-8")
    print("Exact rho_tot_wo_smg/p_tot_wo_smg diagnostics exposed.")


def ini_text(pars: str, root: str) -> str:
    lines = [
        "H0 = 70.8514",
        "omega_b = 0.02239952",
        f"omega_cdm = {OMEGA_CDM:.17g}",
        "N_ncdm = 0",
        "N_ur = 3.046",
        "T_cmb = 2.7255",
        "YHe = 0.2453",
        "Omega_Lambda = 0",
        f"Omega_fld = {OMEGA_FLD:.12g}",
        "fluid_equation_of_state = SDMC_TRACKER",
        "cs2_fld = 0.003",
        "use_ppf = no",
        "Omega_smg = -1",
        "gravity_model = sdmc_kp_transition_bridge",
        f"parameters_smg = {pars}",
        "expansion_model = sdmc_full",
        f"expansion_smg = {EXPANSION}",
        "pert_initial_conditions_smg = zero",
        "method_qs_smg = fully_dynamic",
        "skip_stability_tests_smg = yes",
        "output_background_smg = 3",
        "write background = yes",
        f"root = {root}",
        "format = class",
        "write parameters = yes",
        "input_verbose = 0",
        "background_verbose = 0",
        "thermodynamics_verbose = 0",
        "perturbations_verbose = 0",
        "output_verbose = 0",
    ]
    return "\n".join(lines) + "\n"


def read_table(path: Path) -> list[dict[str, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = [x for x in lines if x.startswith("#") and "1:z" in x]
    if not headers:
        raise RuntimeError(f"could not find numbered header in {path}")
    header = headers[-1].lstrip("#").strip()
    marks = list(re.finditer(r"(\d+)\s*:\s*", header))
    names: list[str] = []
    for i, mark in enumerate(marks):
        end = marks[i + 1].start() if i + 1 < len(marks) else len(header)
        names.append(header[mark.end() : end].strip())
    rows: list[dict[str, float]] = []
    for line in lines:
        if line.strip() and not line.startswith("#"):
            vals = [float(v) for v in line.split()]
            rows.append(dict(zip(names, vals)))
    return rows


def negative_intervals(rows: list[dict[str, float]]) -> str:
    rows = sorted(rows, key=lambda x: x["z"])
    intervals: list[str] = []
    start = None
    for i, row in enumerate(rows):
        isneg = row["cs2num"] < 0.0
        if isneg and start is None:
            start = row["z"]
        if start is not None and ((not isneg) or i == len(rows) - 1):
            end = rows[i - 1]["z"] if not isneg else row["z"]
            intervals.append(f"{start:.9g}:{end:.9g}")
            start = None
    return "|".join(intervals)


def run() -> None:
    OUT.mkdir(exist_ok=True)
    all_rows: list[dict[str, float | str]] = []
    summary: list[dict[str, float | str | int]] = []

    for name, pars in CASES:
        root = f"output/{name}_"
        ini = OUT / f"{name}.ini"
        ini.write_text(ini_text(pars, root), encoding="utf-8")
        cp = subprocess.run(
            ["./class", str(ini)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=90,
        )
        (OUT / f"{name}.log").write_text(cp.stdout, encoding="utf-8")
        bg = OUT / f"{name}_00_background.dat"
        if cp.returncode != 0 or not bg.exists():
            summary.append({"case": name, "status": "RUN_FAIL", "returncode": cp.returncode})
            print(name, "failed", cp.returncode, "background", bg.exists())
            continue

        case_rows: list[dict[str, float | str]] = []
        for r in read_table(bg):
            z = r["z"]
            a = 1.0 / (1.0 + z)
            H = r["H [1/Mpc]"]
            M2 = r["M*^2_smg"]
            bra = r["braiding_smg"]
            run_m = r["M2_running_smg"]
            ten = r["tensor_excess_smg"]
            beh = r["beyond_horndeski_smg"]
            rho_smg = r["(.)rho_smg"]
            p_smg = r["(.)p_smg"]
            rho_wo = r["rho_tot_wo_smg_dbg"]
            p_wo = r["p_tot_wo_smg_dbg"]
            bra_p = r["braiding_prime_smg"]

            # Native hi_class cs2 numerator.  The SDMC bridge is Horndeski,
            # so beyond_horndeski_smg == 0 and its derivative is identically 0.
            t_mix = 0.5 * (2.0 - bra) * (
                bra + 2.0 * beh + 2.0 * run_m + 2.0 * beh * run_m
                - 2.0 * ten + bra * ten
            )
            t_smg = 1.5 * (2.0 - bra) * (1.0 + beh) * (rho_smg + p_smg) / (H * H)
            t_non = -1.5 * (1.0 + beh) * (
                2.0 + 2.0 * beh - 2.0 * M2 + bra * M2
            ) / M2 * (rho_wo + p_wo) / (H * H)
            t_brap = (1.0 + beh) * bra_p / (a * H)
            t_behp = 0.0
            total = t_mix + t_smg + t_non + t_brap + t_behp

            rec: dict[str, float | str] = {
                "case": name,
                "z": z,
                "a": a,
                "cs2num": r["cs2num"],
                "D": r["kin (D)"],
                "cs2": r["c_s^2"],
                "term_mix": t_mix,
                "term_smg_rhop": t_smg,
                "term_non_smg": t_non,
                "term_bra_prime": t_brap,
                "term_beh_prime": t_behp,
                "sum_terms": total,
                "closure_error": total - r["cs2num"],
                "M2": M2,
                "alphaM": run_m,
                "alphaB": bra,
                "alphaB_prime_tau": bra_p,
                "rho_smg": rho_smg,
                "p_smg": p_smg,
                "rho_wo": rho_wo,
                "p_wo": p_wo,
            }
            case_rows.append(rec)
            all_rows.append(rec)

        for label, zmax in (("zlt100", 100.0), ("all", float("inf"))):
            sub = [r for r in case_rows if float(r["z"]) <= zmax]
            if not sub:
                continue
            min_num = min(sub, key=lambda r: float(r["cs2num"]))
            min_cs = min(sub, key=lambda r: float(r["cs2"]))
            summary.append(
                {
                    "case": name,
                    "status": "OK_" + label,
                    "returncode": 0,
                    "min_cs2num": min_num["cs2num"],
                    "z_at_min_cs2num": min_num["z"],
                    "min_cs2": min_cs["cs2"],
                    "z_at_min_cs2": min_cs["z"],
                    "negative_intervals": negative_intervals(sub),
                    "max_closure_error": max(abs(float(r["closure_error"])) for r in sub),
                }
            )

    fields = [
        "case", "z", "a", "cs2num", "D", "cs2", "term_mix", "term_smg_rhop",
        "term_non_smg", "term_bra_prime", "term_beh_prime", "sum_terms", "closure_error",
        "M2", "alphaM", "alphaB", "alphaB_prime_tau", "rho_smg", "p_smg", "rho_wo", "p_wo",
    ]
    with (OUT / "sdmc_cs2_term_audit_v2.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)

    summary_fields = [
        "case", "status", "returncode", "min_cs2num", "z_at_min_cs2num",
        "min_cs2", "z_at_min_cs2", "negative_intervals", "max_closure_error",
    ]
    with (OUT / "sdmc_cs2_term_summary_v2.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows([{k: row.get(k, "") for k in summary_fields} for row in summary])

    print("SUMMARY")
    for row in summary:
        print(row)
    print(
        "MAX_CLOSURE_ERROR_ALL",
        max(abs(float(r["closure_error"])) for r in all_rows) if all_rows else None,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("instrument", "run"))
    args = ap.parse_args()
    if args.mode == "instrument":
        instrument()
    else:
        run()


if __name__ == "__main__":
    main()
