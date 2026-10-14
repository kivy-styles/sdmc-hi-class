#!/usr/bin/env python3
"""Numerical regularization for the decoupled tail of sdmc_early_cs.

The physical tracker kineticity is O(1-10) at early times and tends rapidly to
zero after the frozen handoff.  Cubic-spline interpolation of an exactly
vanishing asymptote can undershoot by ~1e-15 and trigger hi_class's ghost test.
This patch imposes alpha_K >= 1e-8 only in that asymptotic decoupled tail.
The floor is many orders of magnitude below the physical early-time tracker
kineticity and is therefore a numerical regularizer, not a fitted parameter.
"""
from pathlib import Path

p = Path("gravity_smg/gravity_models_smg.c")
text = p.read_text(encoding="utf-8")
old = """    pvecback[pba->index_bg_kineticity_smg] = 3.*enthalpy_smg/(H2*cs2_target);\n"""
new = """    double alphaK_sdmc = 3.*enthalpy_smg/(H2*cs2_target);\n    if (alphaK_sdmc < 1.e-8) alphaK_sdmc = 1.e-8;\n    pvecback[pba->index_bg_kineticity_smg] = alphaK_sdmc;\n"""
if new in text:
    print("sdmc_early_cs alpha_K floor already present")
elif old in text:
    p.write_text(text.replace(old, new, 1), encoding="utf-8")
    print("Applied sdmc_early_cs alpha_K numerical floor: 1e-8")
else:
    raise RuntimeError("sdmc_early_cs kineticity assignment not found")
