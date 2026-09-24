#!/usr/bin/env python3
"""
Patch the CLASS background integration endpoint into the future for a dedicated
SDMC background-only audit.

Upstream CLASS/hi_class fixes log(a/a0)_final=0, i.e. today.  This experiment
changes that endpoint and releases the background redshift table from the
standard z>=0 clamp so future rows carry their physical -1<z<0 values.  It is
intentionally applied in the workflow only after all accepted z>=0 products
have been generated.
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
s=s.replace(old,new,1)

# Standard CLASS clamps all post-present background rows to z=0 because the
# production table ends at today.  A genuine future table must instead remain
# one-to-one in log(a), z and tau; otherwise the z/tau spline receives a long
# repeated z=0 plateau and the SMG post-processing interpolation becomes
# ill-defined.
zold='  pba->z_table[index_loga] = MAX(0.,1./a-1.);'
znew='  pba->z_table[index_loga] = 1./a-1.; /* SDMC FUTURE AUDIT: allow -1<z<0 */'
if zold in s:
    if s.count(zold)!=1:
        raise RuntimeError(f"expected one redshift clamp anchor, found {s.count(zold)}")
    s=s.replace(zold,znew,1)
elif znew not in s:
    raise RuntimeError("future redshift-table anchor not found")

p.write_text(s)
print("FUTURE_BACKGROUND_ENDPOINT",args.loga_final)
print("FUTURE_BACKGROUND_REDSHIFT_UNCLAMP installed")
