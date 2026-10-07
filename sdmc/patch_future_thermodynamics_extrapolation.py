#!/usr/bin/env python3
"""
Add a minimal post-today thermodynamics continuation for the isolated future
perturbation audit.

CLASS thermodynamics tables are observational tables with z>=0.  The future
background audit reaches -1<z<0, so a direct perturbation evolution would
otherwise query outside the thermodynamics interpolation range.

For z<0 this patch:
  * evaluates the existing thermodynamics table at z=0,
  * keeps the present ionization fraction,
  * dilutes Thomson opacity as a^-2,
  * cools baryons adiabatically as T_b ~ a^-2,
  * scales w_b and c_b^2 as a^-2,
  * sets visibility-source functions to zero in the future.

This is a controlled late-time diagnostic continuation, not a prediction of
future atomic/astrophysical thermodynamics.
"""
from pathlib import Path

P=Path("source/thermodynamics.c")
s=P.read_text()
anchor="""  /* The fact that z is in the pre-computed range 0 <= z <= z_initial will be checked in the interpolation routines below. Before
     trying to interpolate, allow the routine to deal with the case z > z_initial: then, all relevant quantities can be extrapolated
     using simple analytic approximations */

  if (z >= pth->z_table[pth->tt_size-1]) {
"""
insert="""  /* Future-only SDMC perturbation audit.  The standard thermodynamics table
     ends at z=0.  For -1<z<0, first load the z=0 state and then apply a
     transparent late-time dilution/cooling continuation. */
  if (z < 0.) {
    double z_future = z;
    double afac = 1. + z_future; /* = 1/a */
    int index_today = 0;
    double Tb0, wb0, cb20, dk0;
    double aH, aH_prime;

    class_test(afac <= 0.,
               pth->error_message,
               "future thermodynamics requested at z=%e <= -1",z_future);

    class_call(thermodynamics_at_z(pba,
                                   pth,
                                   0.,
                                   inter_normal,
                                   &index_today,
                                   pvecback,
                                   pvecthermo),
               pth->error_message,
               pth->error_message);

    Tb0 = pvecthermo[pth->index_th_Tb];
    wb0 = pvecthermo[pth->index_th_wb];
    cb20 = pvecthermo[pth->index_th_cb2];
    dk0 = pvecthermo[pth->index_th_dkappa];

    pvecthermo[pth->index_th_dkappa] = dk0*afac*afac;
    pvecthermo[pth->index_th_ddkappa] =
      -2.*pvecback[pba->index_bg_H]/afac*
      pvecthermo[pth->index_th_dkappa];
    pvecthermo[pth->index_th_dddkappa] =
      (pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]/afac
       -pvecback[pba->index_bg_H_prime])*2./afac*
      pvecthermo[pth->index_th_dkappa];

    /* No future line-of-sight visibility source is required for this
       non-CMB diagnostic. */
    pvecthermo[pth->index_th_exp_m_kappa] = 1.;
    pvecthermo[pth->index_th_g] = 0.;
    pvecthermo[pth->index_th_dg] = 0.;
    pvecthermo[pth->index_th_ddg] = 0.;

    /* Adiabatically decoupled non-relativistic baryon gas: T_b ~ a^-2. */
    pvecthermo[pth->index_th_Tb] = Tb0*afac*afac;
    pvecthermo[pth->index_th_dTb] = 2.*Tb0*afac;
    pvecthermo[pth->index_th_wb] = wb0*afac*afac;
    pvecthermo[pth->index_th_cb2] = cb20*afac*afac;

    if (pth->compute_cb2_derivatives == _TRUE_) {
      aH = pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H];
      aH_prime =
        pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H_prime]
        + aH*aH;
      pvecthermo[pth->index_th_dcb2] =
        -2.*aH*pvecthermo[pth->index_th_cb2];
      pvecthermo[pth->index_th_ddcb2] =
        (-2.*aH_prime+4.*aH*aH)*pvecthermo[pth->index_th_cb2];
    }

    pvecthermo[pth->index_th_rate] =
      fabs(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);
    *last_index = 0;
    return _SUCCESS_;
  }

  if (z >= pth->z_table[pth->tt_size-1]) {
"""
if "Future-only SDMC perturbation audit" not in s:
    if s.count(anchor)!=1:
        raise RuntimeError("thermodynamics insertion anchor not unique")
    s=s.replace(anchor,insert,1)
P.write_text(s)
print("FUTURE_THERMODYNAMICS_EXTRAPOLATION_PATCH")
