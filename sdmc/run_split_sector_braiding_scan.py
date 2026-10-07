#!/usr/bin/env python3
"""Scan the split-sector late Horndeski braiding relation.

This is a diagnostic around the exact-No-Slip split-sector completion.  The
background and the early tracker fluid are held fixed.  Only the late scalar
relation

    alpha_B = beta_B alpha_M

is varied.  For each beta_B we re-close alpha_K from the native hi_class sound-
speed numerator already used by the split-sector patch, retain the native
stability checks, and run the full-dynamic TT/TE/EE+lensing calculation when
stable.

The scan is not allowed to redefine either frozen branch.  beta_B=-2 is the
exact-No-Slip point and is kept explicitly as the reference test.  Any other
stable beta_B is only evidence about completion freedom.
"""
from __future__ import annotations

import csv
import glob
import math
import os
import re
import subprocess
from pathlib import Path

SOURCE = Path("gravity_smg/gravity_models_smg.c")
TEMPLATES = {
    "Kp_core": Path("sdmc/config/kp_core_split_sector.ini"),
    "NKp_v2": Path("sdmc/config/nkp_v2_split_sector.ini"),
}
REFS = {
    "Kp_core": Path("output/sdmc_kp_tracker_clustered_cs0003_00_cl_lensed.dat"),
    "NKp_v2": Path("output/sdmc_nkp_v2_tracker_clustered_cs0003_00_cl_lensed.dat"),
}

BETAS = [-6.0,-5.0,-4.0,-3.5,-3.0,-2.75,-2.5,-2.25,-2.10,-2.05,-2.00,
         -1.95,-1.90,-1.75,-1.50,-1.00,-0.50,0.0,0.50,1.0,2.0,3.0,4.0,5.0,6.0]

FIELDS = [
    "case","beta_B","status","returncode","min_cs2","min_D","F0",
    "alphaM0","alphaM_min","alphaM_max","z_alphaM_max",
    "TT_ratio_l200","TT_ratio_l1000","TT_ratio_l2000",
    "EE_ratio_l200","EE_ratio_l1000","EE_ratio_l2000",
    "phiphi_ratio_l100","phiphi_ratio_l500","phiphi_ratio_l1000","phiphi_ratio_l1500",
    "error"
]


def set_line(text: str, key: str, value: str) -> str:
    pat = re.compile(rf"(?m)^\s*{re.escape(key)}\s*=.*$")
    if not pat.search(text):
        raise RuntimeError(f"missing config key {key}")
    return pat.sub(f"{key} = {value}", text)


def clean(prefix: str) -> None:
    for p in glob.glob(prefix + "*"):
        try:
            os.remove(p)
        except OSError:
            pass


def table(path: str | Path):
    rows=[]
    with open(path,encoding="utf-8") as f:
        for line in f:
            s=line.strip()
            if s and not s.startswith("#"):
                rows.append([float(x) for x in s.split()])
    return rows


def find1(pattern: str) -> str:
    xs=sorted(glob.glob(pattern))
    if len(xs)!=1:
        raise RuntimeError(f"{pattern}: found {len(xs)}")
    return xs[0]


def nearest(rows, ell: int, col: int) -> float:
    return min(rows,key=lambda r:abs(r[0]-ell))[col]


def background_metrics(path: str | Path):
    b=table(path)
    z0=min(b,key=lambda r:abs(r[0]))
    ammax=max(b,key=lambda r:r[27])
    return {
        "min_cs2":min(r[29] for r in b),
        "min_D":min(r[30] for r in b),
        "F0":z0[22],
        "alphaM0":z0[27],
        "alphaM_min":min(r[27] for r in b),
        "alphaM_max":ammax[27],
        "z_alphaM_max":ammax[0],
    }


def rebuild() -> tuple[int,str]:
    cp=subprocess.run(["make","-j2","class"],text=True,stdout=subprocess.PIPE,
                      stderr=subprocess.STDOUT,timeout=240)
    return cp.returncode,cp.stdout


def main():
    original=SOURCE.read_text(encoding="utf-8")
    # These two anchors occur once per late split-sector closure (v2 and Kp).
    if original.count("double bra = -2.*am;") != 2:
        raise RuntimeError("expected two split-sector alpha_B=-2 alpha_M anchors")
    if original.count("double bra_N = -2.*amp;") != 2:
        raise RuntimeError("expected two split-sector alpha_B'=-2 alpha_M' anchors")

    refs={k:table(v) for k,v in REFS.items()}
    rows=[]
    try:
        for beta in BETAS:
            mod=original.replace("double bra = -2.*am;",
                                 f"double bra = ({beta:.17g})*am; /* split-sector beta scan */")
            mod=mod.replace("double bra_N = -2.*amp;",
                            f"double bra_N = ({beta:.17g})*amp; /* split-sector beta scan */")
            SOURCE.write_text(mod,encoding="utf-8")
            rc,buildlog=rebuild()
            if rc:
                for case in TEMPLATES:
                    row={k:"" for k in FIELDS}; row.update(case=case,beta_B=beta,status="BUILD_FAIL",returncode=rc,
                        error=" | ".join(buildlog.strip().splitlines()[-20:])[:4000])
                    rows.append(row)
                print(beta,"BUILD_FAIL",flush=True)
                continue

            for case,tpl in TEMPLATES.items():
                row={k:"" for k in FIELDS}; row.update(case=case,beta_B=f"{beta:.12g}")
                tag=(f"{beta:+.2f}").replace("+","p").replace("-","m").replace(".","p")
                prefix=f"output/sdmc_splitbraid_{case}_{tag}_"
                ini=Path(f"output/sdmc_splitbraid_{case}_{tag}.ini")
                log=Path(f"output/sdmc_splitbraid_{case}_{tag}.log")
                clean(prefix)
                try:
                    text=tpl.read_text(encoding="utf-8")
                    text=set_line(text,"root",prefix)
                    text=set_line(text,"P_k_max_h/Mpc","0.5")
                    text=set_line(text,"z_pk","0")
                    for key in ("input_verbose","background_verbose","thermodynamics_verbose","perturbations_verbose","spectra_verbose","lensing_verbose","output_verbose"):
                        text=set_line(text,key,"0")
                    ini.write_text(text,encoding="utf-8")
                    cp=subprocess.run(["./class",str(ini)],text=True,stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT,timeout=300)
                    log.write_text(cp.stdout,encoding="utf-8")
                    row["returncode"]=cp.returncode
                    if cp.returncode:
                        row["status"]="FAIL"
                        row["error"]=" | ".join(cp.stdout.strip().splitlines()[-24:])[:4000]
                        rows.append(row)
                        print(case,beta,"FAIL",flush=True)
                        continue

                    bmet=background_metrics(find1(prefix+"*_background.dat"))
                    cl=table(find1(prefix+"*_cl_lensed.dat"))
                    ref=refs[case]
                    row.update(status="OK",**{k:f"{v:.12g}" for k,v in bmet.items()})
                    for ell in (200,1000,2000):
                        row[f"TT_ratio_l{ell}"]=f"{nearest(cl,ell,1)/nearest(ref,ell,1):.12g}"
                        row[f"EE_ratio_l{ell}"]=f"{nearest(cl,ell,3)/nearest(ref,ell,3):.12g}"
                    for ell in (100,500,1000,1500):
                        row[f"phiphi_ratio_l{ell}"]=f"{nearest(cl,ell,5)/nearest(ref,ell,5):.12g}"
                    rows.append(row)
                    print(case,beta,"OK","min_cs2",row["min_cs2"],"phi1000/ref",row["phiphi_ratio_l1000"],flush=True)
                except subprocess.TimeoutExpired as exc:
                    row["status"]="TIMEOUT"; row["error"]=f"timeout after {exc.timeout}s"; rows.append(row)
                finally:
                    clean(prefix)
                    try: ini.unlink()
                    except OSError: pass
    finally:
        SOURCE.write_text(original,encoding="utf-8")
        rc,buildlog=rebuild()
        if rc:
            raise RuntimeError("failed to restore exact No-Slip source: "+" | ".join(buildlog.strip().splitlines()[-20:]))

    out=Path("output/sdmc_split_sector_braiding_scan.csv")
    with out.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    print("Wrote",out)


if __name__=="__main__":
    main()
