#!/usr/bin/env python3
"""
Patch the CLASS background integration endpoint into the future for a dedicated
SDMC background-only audit.

Upstream CLASS/hi_class fixes log(a/a0)_final=0, i.e. today.  This experiment
changes only that endpoint.  It is intentionally applied in the workflow after
all accepted z>=0 products have been generated and is reverted before the
ordinary structural-clock basin test.
"""
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument("--loga-final",type=float,default=5.0)
args=ap.parse_args()
if not (0. < args.loga_final <= 10.):
    raise SystemExit("require 0 < loga-final <= 10 for this audit")

p=Path("source/background.c")
s=p.read_text()
old='  loga_final = 0.; // with our conventions, loga is in fact log(a/a_0); we integrate until today, when log(a/a_0) = 0'
new=(f'  loga_final = {args.loga_final:.17e}; '
     '// SDMC FUTURE AUDIT: extend background table beyond a/a0=1')
if old not in s:
    if new in s:
        print("FUTURE_BACKGROUND_ENDPOINT already installed")
        raise SystemExit(0)
    raise RuntimeError("unique CLASS loga_final anchor not found")
if s.count(old)!=1:
    raise RuntimeError(f"expected one loga_final anchor, found {s.count(old)}")
p.write_text(s.replace(old,new,1))
print("FUTURE_BACKGROUND_ENDPOINT",args.loga_final)
