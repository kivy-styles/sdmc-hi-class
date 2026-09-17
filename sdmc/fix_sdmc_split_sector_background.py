#!/usr/bin/env python3
"""Fix background bookkeeping for the SDMC split-sector Boltzmann test.

Applied *after* apply_sdmc_split_sector_patch.py.

The first split-sector draft moved the early tracker perturbations into a
standard CLASS fluid, but sdmc_full still built its target H(a) from rho_tot
that already contained that fluid.  This counted the tracker carrier twice in
the target background.  In addition, CLASS's early-radiation consistency test
only counts components accumulated into rho_de; the custom tracker fluid was
not included there even though it is part of the SDMC structural sector.

This patch does two bookkeeping corrections without changing the frozen target
H(a):

1. For fluid_equation_of_state=SDMC_TRACKER, include rho_fld in rho_de for
   CLASS's background fractions/initial radiation check.
2. When sdmc_full constructs rho_late, remove the explicit tracker carrier from
   the base density/pressure first.  The desired tracker factor 1/(1-f) is then
   applied once, and rho_smg is the residual needed after the separately
   conserved tracker fluid has been added.

For all non-split configurations the original sdmc_full equations are
unchanged.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: split bookkeeping fix already present")
        return
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{path}: expected one fix anchor, found {n}\n{old}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: split bookkeeping fixed")


# Treat the explicit tracker carrier as structural/DE for the purpose of the
# CLASS Omega_de bookkeeping and the early radiation-dominance check.
replace_once(
    "source/background.c",
    '''    rho_tot += pvecback[pba->index_bg_rho_fld];
    p_tot += w_fld * pvecback[pba->index_bg_rho_fld];
    dp_dloga += (a*dw_over_da-3*(1+w_fld)*w_fld)*pvecback[pba->index_bg_rho_fld];''',
    '''    rho_tot += pvecback[pba->index_bg_rho_fld];
    p_tot += w_fld * pvecback[pba->index_bg_rho_fld];
    dp_dloga += (a*dw_over_da-3*(1+w_fld)*w_fld)*pvecback[pba->index_bg_rho_fld];
    if (pba->fluid_equation_of_state == SDMC_TRACKER)
      rho_de += pvecback[pba->index_bg_rho_fld];'''
)


# In sdmc_full, construct the frozen late+tracker target from the ordinary
# matter/radiation base, not from a rho_tot that already contains the separate
# tracker carrier.
replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    double rho_X = Omega_X0*pow(pba->H0,2);
    double rho_flat = rho_tot + rho_X;
    double rho_non_N = -3.*(rho_tot + p_tot);
    double rho_flat_N = rho_non_N;''',
    '''    double rho_X = Omega_X0*pow(pba->H0,2);

    /* Split-sector runs carry the early tracker explicitly as rho_fld.
     * Remove that carrier before constructing the frozen sdmc_full target,
     * otherwise the tracker is counted both here and through 1/(1-f). */
    double rho_base = rho_tot;
    double p_base = p_tot;
    if (pba->has_fld == _TRUE_ &&
        pba->fluid_equation_of_state == SDMC_TRACKER) {
      rho_base -= pvecback[pba->index_bg_rho_fld];
      p_base -= pvecback[pba->index_bg_w_fld]
                *pvecback[pba->index_bg_rho_fld];
    }

    double rho_flat = rho_base + rho_X;
    double rho_base_N = -3.*(rho_base + p_base);
    double rho_all_non_smg_N = -3.*(rho_tot + p_tot);
    double rho_flat_N = rho_base_N;'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    double rho_smg = rho_target-rho_tot;
    double rho_smg_N = rho_target_N-rho_non_N;
    double p_smg = -rho_smg-rho_smg_N/3.;''',
    '''    double rho_smg = rho_target-rho_tot;
    double rho_smg_N = rho_target_N-rho_all_non_smg_N;
    double p_smg = -rho_smg-rho_smg_N/3.;'''
)

print("SDMC split-sector background bookkeeping fix complete.")
