#!/usr/bin/env python3
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument("--loga-stop",type=float,required=True)
args=ap.parse_args()
if not (0.0 < args.loga_stop <= 10.0):
    raise SystemExit("require 0 < loga-stop <= 10")

p=Path("source/perturbations.c")
s=p.read_text()
old='''  /* Future SDMC audit: retain the observational conformal_age for CMB
     calculations, but allow non-CMB perturbation diagnostics to follow the
     already-computed background table beyond a=1. */
  double tau_end_future_audit =
    (ppt->has_cmb == _TRUE_) ? pba->conformal_age
                             : pba->tau_table[pba->bt_size-1];
'''
new=f'''  /* Future SDMC audit: retain the observational conformal_age for CMB
     calculations.  For non-CMB diagnostics, stop at an interior requested
     log(a) rather than forcing the perturbation solver onto the final
     background-table boundary. */
  double tau_end_future_audit = pba->conformal_age;
  if (ppt->has_cmb == _FALSE_) {{
    class_call(background_tau_of_z(pba,
                                   exp(-({args.loga_stop:.17e}))-1.,
                                   &tau_end_future_audit),
               pba->error_message,
               ppt->error_message);
  }}
'''
if old not in s:
    if "stop at an interior requested" in s:
        print("FUTURE_PERTURBATION_LOGA_STOP already installed")
    else:
        raise RuntimeError("future perturbation endpoint anchor not found; apply patch_future_perturbation_endpoint.py first")
else:
    s=s.replace(old,new,1)
    p.write_text(s)
print("FUTURE_PERTURBATION_LOGA_STOP",args.loga_stop)
