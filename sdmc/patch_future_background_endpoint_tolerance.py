#!/usr/bin/env python3
from pathlib import Path

p=Path("source/background.c")
s=p.read_text()
old='''  class_test(loga > pba->loga_table[pba->bt_size-1],
             pba->error_message,
             "out of range: a/a_0 = %e > a_max/a_0 = %e\\n",1./(1.+z),exp(pba->loga_table[pba->bt_size-1]));
'''
new='''  /* SDMC FUTURE AUDIT: floating-point endpoint guard.  Arbitrary positive
     future endpoints can round a reconstructed log(a) a few ulps above the
     final tabulated value.  Clamp only this numerically indistinguishable
     overshoot; genuine out-of-range requests remain errors. */
  if ((loga > pba->loga_table[pba->bt_size-1]) &&
      (loga-pba->loga_table[pba->bt_size-1] < 1.e-12)) {
    loga = pba->loga_table[pba->bt_size-1];
  }

  class_test(loga > pba->loga_table[pba->bt_size-1],
             pba->error_message,
             "out of range: a/a_0 = %e > a_max/a_0 = %e\\n",1./(1.+z),exp(pba->loga_table[pba->bt_size-1]));
'''
if new in s:
    print("FUTURE_BACKGROUND_ENDPOINT_TOLERANCE already installed")
elif old in s:
    s=s.replace(old,new,1)
else:
    raise RuntimeError("background_at_z endpoint anchor not found")
p.write_text(s)
print("FUTURE_BACKGROUND_ENDPOINT_TOLERANCE installed")
