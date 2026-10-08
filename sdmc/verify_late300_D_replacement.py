#!/usr/bin/env python3
"""Strict provenance and free-evolution checks for *real* late300 D replacements.

Unlike old replay work, inspect source target INI -> actual target history
-> reconstructed G2 coefficients -> free Horndeski evolution. This verifies
data flow and consistency only, not independent microscopic origins.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import numpy as np

def ini_values(path,key):
    hits=[l.split("=",1)[1].strip() for l in path.read_text().splitlines()
          if l.split("#",1)[0].strip().startswith(key+" ") and "=" in l]
    if len(hits)!=1: raise ValueError(f"Expected exactly one {key} in {path}: {len(hits)}")
    return [float(x) for x in hits[0].split(",")]

def bgread(path):
    lines=path.read_text().splitlines()
    h=[s for s in lines if s.startswith("#") and re.search(r"1\s*:",s)]
    if not h:raise RuntimeError(f"No numbered background header: {path}")
    h=h[-1].lstrip("#").strip()
    marks=list(re.finditer(r"(\d+)\s*:\s*",h))
    cols=[h[m.end():marks[i+1].start() if i+1<len(marks) else len(h)].strip()
          for i,m in enumerate(marks)]
    a=np.atleast_2d(np.loadtxt(path))
    return {key:a[:,i] for i,key in enumerate(cols)}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def close(x,y,atol=3e-7,rtol=2e-6):
    return abs(x-y) <= atol+rtol*abs(y)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--case",required=True)
    ap.add_argument("--output",default="output/D_replacement_verdict.json")
    ap.add_argument("--validate-free",action="store_true")
    args=ap.parse_args()
    target=Path("output/linear_cov_target.ini")
    AF, zcF, wF, D0, pD, DF=ini_values(target,"parameters_smg")
    zcD=3.927876388467848
    wD=0.33114133956842123
    y=bgread(Path("output/linear_cov_target_00_background.dat"))
    z=y["z"]; N=-np.log1p(z)
    Sd=.5*(1+np.tanh((N+np.log1p(zcD))/(2*wD)))
    needD=DF + D0*Sd**pD
    measuredD=y["kin (D)"]
    maxerr=float(np.max(np.abs(measuredD-needD)))
    if maxerr>3e-6:raise AssertionError(f"{args.case}: target parameter did not propagate: max |D_native-D_configured|={maxerr}")
    z0=int(np.argmin(np.abs(z)))
    if abs(z[z0])>1e-7:raise RuntimeError("Missing z=0 node")
    report={
      "case":args.case, "D0_correction":D0, "D_floor_baseline":DF,
      "pD":pD,"F_zc":zcF,"F_width":wF,"D_zc_historic":zcD,"D_width_historic":wD,
      "target_D0_today":float(measuredD[z0]),
      "analytic_D0_today":float(needD[z0]),
      "target_D_early_plateau":float(measuredD[np.argmax(z)]),
      "analytic_D_early_plateau":float(needD[np.argmax(z)]),
      "Dtarget_closure_max_abs":maxerr,
      "target_hubble_today":float(y["H [1/Mpc]"][z0]),
      "target_F_today":float(y["M*^2_smg"][z0]),
      "target_min_D":float(np.min(measuredD)),
      "target_min_cs2":float(np.min(y["c_s^2"])),
      "target_max_cs2":float(np.max(y["c_s^2"])),
      "target_background_sha256":sha(Path("output/linear_cov_target_00_background.dat")),
      "target_spectrum_sha256":sha(Path("output/linear_cov_target_00_cl_lensed.dat")),
      "interpretation":"D0 is a correction to Dfloor; replacing parameters in reconstruction inputs is a consistency test, not a first-principles prediction."
    }
    if args.validate_free:
      x=bgread(Path("output/linear_cov_free_00_background.dat"))
      zv=x["z"];Nv=-np.log1p(zv)
      o=np.argsort(N);Nsort=N[o]
      Href=np.interp(Nv,Nsort,y["H [1/Mpc]"][o])
      Fref=np.interp(Nv,Nsort,y["M*^2_smg"][o])
      Dref=np.interp(Nv,Nsort,y["kin (D)"][o])
      ii=zv<=100
      report.update({
        "free_max_abs_dH_H_z100":float(np.max(np.abs(x["H [1/Mpc]"][ii]/Href[ii]-1))),
        "free_max_abs_dF_F_z100":float(np.max(np.abs(x["M*^2_smg"][ii]/Fref[ii]-1))),
        "free_max_abs_D_err_z100":float(np.max(np.abs(x["kin (D)"][ii]-Dref[ii]))),
        "free_min_D_z100":float(np.min(x["kin (D)"][ii])),
        "free_min_cs2_z100":float(np.min(x["c_s^2"][ii])),
        "free_max_cs2_z100":float(np.max(x["c_s^2"][ii])),
        "free_spectrum_sha256":sha(Path("output/linear_cov_free_00_cl_lensed.dat")),
        "free_background_sha256":sha(Path("output/linear_cov_free_00_background.dat")),
      })
      if not 0. < report["free_min_cs2_z100"] <= 1.:
        raise AssertionError(f"{args.case}: reconstructed scalar gradient instability")
      if report["free_max_cs2_z100"]>1.0:
        raise AssertionError(f"{args.case}: reconstructed scalar superluminality")
      if report["free_min_D_z100"]<=0:raise AssertionError(f"{args.case}: scalar kinetic ghost")
      if report["free_max_abs_dH_H_z100"]>1e-3 or report["free_max_abs_dF_F_z100"]>1e-3:
        raise AssertionError(f"{args.case}: free action does not reproduce target")
      action_summary=Path("output/linear_cov_reconstruction_summary.txt").read_text()
      found=dict(re.findall(r"^(D0|Dfloor|Dclosure_max_rel)=(\S+)",action_summary,re.M))
      if not {"D0","Dfloor","Dclosure_max_rel"}<=set(found):
        raise RuntimeError("Action summary has no D0/Dfloor provenance")
      if not close(float(found["D0"]),D0,atol=1e-12) or not close(float(found["Dfloor"]),DF,atol=1e-12):
        raise AssertionError(f"{args.case}: reconstructed action used different D0/Dfloor: {found}")
      report["reconstructed_action_D0_from_summary"]=float(found["D0"])
      report["reconstructed_action_Dfloor_from_summary"]=float(found["Dfloor"])
      report["action_relative_D_closure"]=float(found["Dclosure_max_rel"])
      if report["action_relative_D_closure"]>1e-6:raise AssertionError("Action fails kinetic closure")
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("TRUE_D_REPLACEMENT_VERIFIED",json.dumps(report,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
