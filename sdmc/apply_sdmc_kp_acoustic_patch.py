#!/usr/bin/env python3
"""Patch CLASS with the frozen SDMC Kp terminal acoustic operator.

This is a controlled Boltzmann diagnostic, not a post-hoc C_ell warp.  The
operator acts only on the photon-baryon acoustic monopole/dipole sector during
the frozen drag-tail window.  It implements the manuscript oscillator form

  Theta0'' + [... - C'/C] Theta0' + C^2 k^2/[3(1+R_Hu)] Theta0 = S,

by (i) multiplying the photon pressure-gradient term by C^2 and (ii) adding
(C'/C) theta to the photon and baryon Euler equations.  Higher photon
multipoles keep the physical free-streaming speed; their polarization response
therefore follows dynamically from the modified dipole/quadrupole source.

The window is a smooth two-edge top hat

  W(z)=1/4 [1+tanh((z_on-z)/dz)] [1+tanh((z-z_off)/dz)],
  C(z)=1-(1-C_min) W(z).

Frozen values are C_min=0.539812 and z_off=1059.0.  z_on and dz are supplied
by the dedicated workflow.  The three z_on values used there are fixed by the
1.208 Mpc residual-tail integral for dz=0.5,1,2; they are not fitted to Planck.

A guard aborts if the material part of the window overlaps CLASS tight
coupling.  That makes this first implementation auditable: either the whole
operator is evolved by the full photon hierarchy, or the run explicitly fails
instead of silently using an inconsistent TCA formula.
"""
from pathlib import Path
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--dz", type=float, required=True)
ap.add_argument("--zon", type=float, required=True)
ap.add_argument("--zoff", type=float, default=1059.0)
ap.add_argument("--cmin", type=float, default=0.539812)
a = ap.parse_args()
if not (0.0 < a.cmin <= 1.0):
    raise SystemExit("cmin must be in (0,1]")
if a.dz <= 0 or a.zon <= a.zoff:
    raise SystemExit("require dz>0 and zon>zoff")

p = Path("source/perturbations.c")
text = p.read_text(encoding="utf-8")
marker = "SDMC_KP_TERMINAL_ACOUSTIC_OPERATOR"
if marker in text:
    raise SystemExit("Kp acoustic operator already patched")


def replace_once(old: str, new: str) -> None:
    global text
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"expected one perturbations.c anchor, found {n}:\n{old[:300]}")
    text = text.replace(old, new, 1)

# Local variables in perturbations_derivs only (the delta_idm continuation makes
# this anchor unique relative to the TCA helper).
old = """  double delta_g=0.,theta_g=0.,shear_g=0.;
  double delta_b,theta_b;
  double delta_idm = 0., theta_idm = 0.;
"""
new = """  double delta_g=0.,theta_g=0.,shear_g=0.;
  double delta_b,theta_b;
  /* SDMC_KP_TERMINAL_ACOUSTIC_OPERATOR */
  double sdmc_kp_C=1., sdmc_kp_C2=1., sdmc_kp_dlnC_dtau=0., sdmc_kp_W=0.;
  double delta_idm = 0., theta_idm = 0.;
"""
replace_once(old, new)

# Evaluate the frozen terminal window after H(a) is available.  Since
# dz/dtau=-H in CLASS units, d ln C/dtau = -H d ln C/dz.
old = """  a = pvecback[pba->index_bg_a];
  a2 = a*a;
  a_prime_over_a = pvecback[pba->index_bg_H] * a;

  /* There are two definitions of R:
"""
new = f"""  a = pvecback[pba->index_bg_a];
  a2 = a*a;
  a_prime_over_a = pvecback[pba->index_bg_H] * a;

  /* SDMC frozen Kp terminal acoustic window. */
  {{
    const double kp_Cmin = {a.cmin:.15g};
    const double kp_zoff = {a.zoff:.15g};
    const double kp_zon  = {a.zon:.15g};
    const double kp_dz   = {a.dz:.15g};
    double kp_z = 1./a - 1.;
    double u_on = (kp_zon-kp_z)/kp_dz;
    double u_off = (kp_z-kp_zoff)/kp_dz;
    double t_on = tanh(u_on), t_off = tanh(u_off);
    double Awin = .5*(1.+t_on);
    double Bwin = .5*(1.+t_off);
    double dA_dz = -.5*(1.-t_on*t_on)/kp_dz;
    double dB_dz =  .5*(1.-t_off*t_off)/kp_dz;
    double dW_dz;
    sdmc_kp_W = Awin*Bwin;
    /* Compact the numerically irrelevant tanh tails. */
    if (sdmc_kp_W < 1.e-12) sdmc_kp_W = 0.;
    dW_dz = dA_dz*Bwin + Awin*dB_dz;
    if (sdmc_kp_W == 0.) dW_dz = 0.;
    sdmc_kp_C = 1.-(1.-kp_Cmin)*sdmc_kp_W;
    sdmc_kp_C2 = sdmc_kp_C*sdmc_kp_C;
    sdmc_kp_dlnC_dtau = pvecback[pba->index_bg_H]
                         *(1.-kp_Cmin)*dW_dz/sdmc_kp_C;

    /* The retained Kp window should lie after CLASS has left TCA.  Abort
       instead of silently applying an inconsistent tight-coupling formula. */
    if ((sdmc_kp_W > 1.e-3) &&
        (ppw->approx[ppw->index_ap_tca] == (int)tca_on)) {{
      class_stop(error_message,
                 \"SDMC Kp acoustic window overlaps tight coupling at z=%g (W=%g); extend the operator consistently into TCA before using this width\",
                 kp_z,sdmc_kp_W);
    }}
  }}

  /* There are two definitions of R:
"""
replace_once(old, new)

# Baryon Euler equation after TCA: the common acoustic clock receives C'/C.
old = """      dy[pv->index_pt_theta_b] =
        - a_prime_over_a*theta_b
        + metric_euler
        + k2*delta_p_b_over_rho_b
        + R*pvecthermo[pth->index_th_dkappa]*(theta_g-theta_b);
"""
new = """      dy[pv->index_pt_theta_b] =
        - a_prime_over_a*theta_b
        + metric_euler
        + k2*delta_p_b_over_rho_b
        + R*pvecthermo[pth->index_th_dkappa]*(theta_g-theta_b)
        + sdmc_kp_dlnC_dtau*theta_b;
"""
replace_once(old, new)

# The TCA branch is guarded above.  Keep its standard equations untouched.

# Photon dipole after TCA: C^2 pressure gradient and the -C'/C oscillator
# friction (implemented as +C'/C theta in the first-order Euler equation).
old = """        dy[pv->index_pt_theta_g] =
          k2*(delta_g/4.-s2_squared*y[pv->index_pt_shear_g])
          + metric_euler
          + photon_scattering_rate*(theta_b-theta_g);
"""
new = """        dy[pv->index_pt_theta_g] =
          k2*(sdmc_kp_C2*delta_g/4.-s2_squared*y[pv->index_pt_shear_g])
          + metric_euler
          + photon_scattering_rate*(theta_b-theta_g)
          + sdmc_kp_dlnC_dtau*theta_g;
"""
replace_once(old, new)

p.write_text(text, encoding="utf-8")
print(f"Applied Kp acoustic operator: Cmin={{a.cmin:.9g}}, zoff={{a.zoff:.9g}}, zon={{a.zon:.9g}}, dz={{a.dz:.9g}}")
