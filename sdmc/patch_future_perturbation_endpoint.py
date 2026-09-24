#!/usr/bin/env python3
"""
Patch the isolated future-audit build so the perturbation source sampling and
mode integration continue to the final future background table entry instead
of stopping at the preserved observational 'today'.

This patch is intended only after patch_future_background_endpoint.py and
patch_future_today_bookkeeping.py.  It does not alter the background table or
the accepted z>=0 cosmology.

CMB runs are deliberately left anchored to the observational conformal age.
For non-CMB future diagnostics, the endpoint becomes tau_table[bt_size-1].
"""
from pathlib import Path

P=Path("source/perturbations.c")
s=P.read_text()
start=s.index("int perturbations_timesampling_for_sources(")
next_start=s.index("\nint perturbations_get_k_list(",start)
block=s[start:next_start]

decl="""  double tau_mid;

  double timescale_source;"""
newdecl="""  double tau_mid;
  /* Future SDMC audit: retain the observational conformal_age for CMB
     calculations, but allow non-CMB perturbation diagnostics to follow the
     already-computed background table beyond a=1. */
  double tau_end_future_audit =
    (ppt->has_cmb == _TRUE_) ? pba->conformal_age
                             : pba->tau_table[pba->bt_size-1];

  double timescale_source;"""
if newdecl not in block:
    if block.count(decl)!=1:
        raise RuntimeError("could not locate tau declaration anchor")
    block=block.replace(decl,newdecl,1)

n=block.count("pba->conformal_age")
# One occurrence remains intentionally inside the definition above.
block=block.replace("pba->conformal_age","tau_end_future_audit")
# Undo the self-reference introduced inside the definition.
block=block.replace(
  "(ppt->has_cmb == _TRUE_) ? tau_end_future_audit",
  "(ppt->has_cmb == _TRUE_) ? pba->conformal_age",
  1
)

s=s[:start]+block+s[next_start:]
P.write_text(s)
print("FUTURE_PERTURBATION_ENDPOINT_PATCH")
print("replaced_conformal_age_references",n-1)
