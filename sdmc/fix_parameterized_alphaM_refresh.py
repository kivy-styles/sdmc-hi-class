#!/usr/bin/env python3
"""Refresh numerically reconstructed alpha_M in pvecback before As/cs2 are built.

In parameterized hi_class runs with field_evolution_smg = false and
M2_evolution_smg = false, the first-pass background explicitly resets
M2_running_smg to zero.  During background_solve_smg the numerical value
alpha_M = d(delta M2)/d ln a / M2 is reconstructed, but upstream code only
copies that value into the background table.  The live pvecback used
immediately afterwards by gravity_functions_As_from_alphas_smg therefore
still contains the stale zero value.  This makes cs2num, the A/lambda
combinations, Geff and slip internally inconsistent with the final stored
M2_running_smg column.

This patch writes the reconstructed alpha_M into pvecback first, then copies
the same value to the background table.  No frozen SDMC background function
or parameter is changed.
"""
from pathlib import Path

p = Path("gravity_smg/background_smg.c")
s = p.read_text(encoding="utf-8")

old = """\tif (pba->field_evolution_smg == _FALSE_ && pba->M2_evolution_smg == _FALSE_){
\t\tdouble alpha_M = pvecback_derivs[pba->index_bg_delta_M2_smg]/pvecback[pba->index_bg_M2_smg];
\t\tcopy_to_background_table_smg(pba, i, pba->index_bg_M2_running_smg, alpha_M);
\t}
"""

new = """\tif (pba->field_evolution_smg == _FALSE_ && pba->M2_evolution_smg == _FALSE_){
\t\tdouble alpha_M = pvecback_derivs[pba->index_bg_delta_M2_smg]/pvecback[pba->index_bg_M2_smg];
\t\t/* Keep the live background vector synchronized with the table before
\t\t * gravity_functions_As_from_alphas_smg() is called in this same loop.
\t\t * Otherwise derived functions are evaluated with the stale first-pass
\t\t * M2_running_smg (zero in parameterized M2 runs). */
\t\tpvecback[pba->index_bg_M2_running_smg] = alpha_M;
\t\tcopy_to_background_table_smg(pba, i, pba->index_bg_M2_running_smg, alpha_M);
\t}
"""

if new in s:
    print("parameterized alpha_M refresh fix already present")
    raise SystemExit(0)

count = s.count(old)
if count != 1:
    raise RuntimeError(f"expected one alpha_M refresh anchor, found {count}")

p.write_text(s.replace(old, new, 1), encoding="utf-8")
print("parameterized alpha_M refresh fix applied")
