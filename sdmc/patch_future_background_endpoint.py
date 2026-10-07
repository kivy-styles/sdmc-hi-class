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
import argparse,re

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
if old in s:
    if s.count(old)!=1:
        raise RuntimeError(f"expected one loga_final anchor, found {s.count(old)}")
    s=s.replace(old,new,1)
elif new in s:
    print("FUTURE_BACKGROUND_ENDPOINT already installed")
else:
    # Allow a later audit step to lengthen/shorten an already isolated future
    # endpoint without requiring a fresh checkout.
    pat=(r'  loga_final = [-+0-9.eE]+; '
         r'// SDMC FUTURE AUDIT: extend background table beyond a/a0=1')
    ms=list(re.finditer(pat,s))
    if len(ms)!=1:
        raise RuntimeError(
          f"unique CLASS future loga_final anchor not found (found {len(ms)})")
    s=s[:ms[0].start()]+new+s[ms[0].end():]

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


# ndf15 can finish at loga_final with the final requested output point missed
# by a floating-point comparison in its output loop.  That is harmless when
# loga_final=0 in the stock code because the endpoint is usually hit exactly,
# but with an arbitrary future endpoint it can leave the final tau/z/background
# row uninitialised.  The SMG post-processing then sees tau_max=0 and fails.
#
# Refresh the final row explicitly from the final integration state.  This does
# not change the trajectory; it only stores the already-computed endpoint.
post_old='''  evolver_ndf15_abstol = 1e-15;

  /** - recover some quantities today */
'''
post_new='''  evolver_ndf15_abstol = 1e-15;

  /* SDMC FUTURE AUDIT: explicitly store the already integrated final row.
     The ndf15 output loop can miss an arbitrary positive loga_final by a
     floating-point comparison, leaving tau_table[bt_size-1] uninitialised. */
  pba->loga_table[pba->bt_size-1] = loga_final;
  pba->z_table[pba->bt_size-1] = 1./exp(loga_final)-1.;
  pba->tau_table[pba->bt_size-1] =
    pvecback_integration[pba->index_bi_tau];
  class_call_except(
    background_functions(
      pba,
      exp(loga_final),
      pvecback_integration,
      long_info,
      pba->background_table+(pba->bt_size-1)*pba->bg_size),
    pba->error_message,
    pba->error_message,
    background_free_noinput(pba);
    free(pvecback);
    free(pvecback_integration);
    free(used_in_output);
  );

  printf("SDMC_FUTURE_ENDPOINT_ROW loga=%e z=%e tau=%e\\n",
         pba->loga_table[pba->bt_size-1],
         pba->z_table[pba->bt_size-1],
         pba->tau_table[pba->bt_size-1]);

  /** - recover some quantities today */
'''
if post_old not in s:
    if "SDMC_FUTURE_ENDPOINT_ROW" not in s:
        raise RuntimeError("could not locate post-evolver endpoint anchor")
else:
    s=s.replace(post_old,post_new,1)

p.write_text(s)
print("FUTURE_BACKGROUND_ENDPOINT",args.loga_final)
print("FUTURE_BACKGROUND_REDSHIFT_UNCLAMP installed")
print("FUTURE_BACKGROUND_FINAL_ROW_REFRESH installed")
