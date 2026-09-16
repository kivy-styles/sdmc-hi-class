#!/usr/bin/env python3
"""Align the NKp-v2 B0 shooting closure with hi_class's native Horndeski c_s^2 numerator.

For alpha_H=alpha_T=0 and exact No-Slip alpha_B=-2 alpha_M, the native
hi_class expression reduces algebraically to

  N_s^hi = -2(1+alpha_M) h - 2 alpha_M'
           - (3 Omega_m + 4 Omega_r)/F,

where prime denotes d/d ln a, h=d ln H/d ln a and F=M_*^2/M_Pl^2.

The first implementation used an extra '+alpha_M' inside the first bracket,
which belongs to a different convention/rearrangement and generated a
spurious native gradient-instability immediately after the B0 activation.
This patch removes that extra term both from the shooting ODE and from the
runtime N_s audit.  It does not alter the frozen background.
"""
from pathlib import Path

p = Path("gravity_smg/gravity_models_smg.c")
text = p.read_text(encoding="utf-8")
old_rhs = "*amp = -.5*(Ns + 2.*(1.+am)*(h+am) + (3.*Om+4.*Or)/F);"
new_rhs = "*amp = -.5*(Ns + 2.*(1.+am)*h + (3.*Om+4.*Or)/F);"
old_ns = "double Ns_actual = -2.*(1.+am)*(h+am)-2.*amp-(3.*Om+4.*Or)/F;"
new_ns = "double Ns_actual = -2.*(1.+am)*h-2.*amp-(3.*Om+4.*Or)/F;"

for old, new in ((old_rhs, new_rhs), (old_ns, new_ns)):
    if new in text:
        continue
    if text.count(old) != 1:
        raise RuntimeError(f"expected one native No-Slip patch anchor: {old}")
    text = text.replace(old, new)

p.write_text(text, encoding="utf-8")
print("SDMC NKp-v2 native hi_class No-Slip convention patch complete.")
