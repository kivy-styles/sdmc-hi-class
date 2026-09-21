#!/usr/bin/env python3
"""Install the Part-VI reciprocal terminal acoustic operator.

Implements the manuscript's phenomenological Euler completion:
  F_K = (C'/C) theta_common + (C^2-1) k^2 Theta0/(1+R_Hu)
with the SAME acceleration added to photon and baryon Euler equations.

CLASS uses R_CLASS=(4/3)rho_gamma/rho_b = 1/R_Hu, hence
  Theta0/(1+R_Hu) = [R_CLASS/(1+R_CLASS)] delta_gamma/4.

The activation coordinate is the manuscript's normalized residual-acoustic
coordinate, evaluated directly from CLASS's background sound horizon:
  x = (r_s(tau)-r_s_star)/(r_s_d-r_s_star).
"""
from pathlib import Path

def replace_once(path,old,new):
    p=Path(path); s=p.read_text()
    if new in s:
        return
    if s.count(old)!=1:
        raise RuntimeError(f"{path}: anchor count {s.count(old)} for {old[:100]!r}")
    p.write_text(s.replace(old,new,1))

# 1) Runtime parameters in perturbations structure.
replace_once(
    "include/perturbations.h",
    """  short has_perturbed_recombination;
  /** Neutrino contribution to tensors */""",
    """  short has_perturbed_recombination;

  /* SDMC Part-VI reciprocal terminal acoustic operator */
  int kp_terminal_enabled;
  double kp_terminal_Cmin;
  double kp_terminal_x1;
  double kp_terminal_x2;
  double kp_terminal_width;

  /** Neutrino contribution to tensors */"""
)

# 2) Defaults.
replace_once(
    "source/input.c",
    """  /** 2) Perturbed recombination */
  ppt->has_perturbed_recombination=_FALSE_;
  /** 3) Modes */""",
    """  /** 2) Perturbed recombination */
  ppt->has_perturbed_recombination=_FALSE_;

  /* SDMC reciprocal terminal operator: disabled is the exact CLASS null. */
  ppt->kp_terminal_enabled = 0;
  ppt->kp_terminal_Cmin = 1.0;
  ppt->kp_terminal_x1 = 0.0207;
  ppt->kp_terminal_x2 = 0.9793;
  ppt->kp_terminal_width = 0.012;

  /** 3) Modes */"""
)

# 3) Input parser.
replace_once(
    "source/input.c",
    """    /** 2) Perturbed recombination */
    /* Read */
    class_read_flag_or_deprecated("perturbed_recombination","perturbed recombination",ppt->has_perturbed_recombination);

    /** 3) Modes */""",
    """    /** 2) Perturbed recombination */
    /* Read */
    class_read_flag_or_deprecated("perturbed_recombination","perturbed recombination",ppt->has_perturbed_recombination);

    /* SDMC Part-VI reciprocal terminal acoustic operator. */
    class_read_int("kp_terminal_enabled",ppt->kp_terminal_enabled);
    class_read_double("kp_terminal_Cmin",ppt->kp_terminal_Cmin);
    class_read_double("kp_terminal_x1",ppt->kp_terminal_x1);
    class_read_double("kp_terminal_x2",ppt->kp_terminal_x2);
    class_read_double("kp_terminal_width",ppt->kp_terminal_width);
    class_test(ppt->kp_terminal_Cmin <= 0. || ppt->kp_terminal_Cmin > 1.,
               errmsg,
               "kp_terminal_Cmin must obey 0 < Cmin <= 1");
    class_test(ppt->kp_terminal_width <= 0.,
               errmsg,
               "kp_terminal_width must be positive");
    class_test(ppt->kp_terminal_x1 < 0. || ppt->kp_terminal_x2 > 1. ||
               ppt->kp_terminal_x1 >= ppt->kp_terminal_x2,
               errmsg,
               "kp terminal edges must obey 0 <= x1 < x2 <= 1");

    /** 3) Modes */"""
)

# 4) Helper at top of perturbations module.
replace_once(
    "source/perturbations.c",
    """#include "parallel.h"


/**""",
    r'''#include "parallel.h"

/*
 * SDMC Part-VI reciprocal terminal profile.
 *
 * The manuscript coordinate x is acoustic distance across the residual
 * star-to-drag interval.  CLASS already evolves r_s(tau), so no redshift or
 * conformal-time surrogate is needed.
 */
static void sdmc_kp_terminal_terms(
    struct background * pba,
    struct thermodynamics * pth,
    struct perturbations * ppt,
    double * pvecback,
    double R_class,
    double k2,
    double delta_g,
    double theta_g,
    double theta_b,
    double * C,
    double * dlnC,
    double * Fk) {

  *C = 1.;
  *dlnC = 0.;
  *Fk = 0.;

  if (ppt->kp_terminal_enabled == 0) return;

  double drs = pth->rs_d - pth->rs_star;
  if (!(drs > 0.)) return;

  double x = (pvecback[pba->index_bg_rs] - pth->rs_star)/drs;
  double w = ppt->kp_terminal_width;
  double u1 = (x-ppt->kp_terminal_x1)/w;
  double u2 = (x-ppt->kp_terminal_x2)/w;
  double t1 = tanh(u1);
  double t2 = tanh(u2);
  double F = 0.5*(t1-t2);
  double dFdx = 0.5*((1.-t1*t1)-(1.-t2*t2))/w;

  *C = 1. - (1.-ppt->kp_terminal_Cmin)*F;
  if (*C < 1.e-8) *C = 1.e-8;

  /*
   * dr_s/deta = c_s with the photon-baryon sound speed.
   * CLASS convention R_class=(4/3)rho_g/rho_b=1/R_Hu.
   */
  double cs2 = R_class/(3.*(1.+R_class));
  if (cs2 < 0.) cs2 = 0.;
  double dx_deta = sqrt(cs2)/drs;
  double Cp = -(1.-ppt->kp_terminal_Cmin)*dFdx*dx_deta;
  *dlnC = Cp/(*C);

  double theta_common =
    (R_class*theta_g + theta_b)/(1.+R_class);

  /*
   * Eq. (136), translated from manuscript R_Hu to CLASS R_class.
   * Theta_0 = delta_gamma/4.
   */
  *Fk = (*dlnC)*theta_common
    + ((*C)*(*C)-1.)*k2*R_class/(1.+R_class)*(delta_g/4.);
}


/**'''
)

# 5) Add the SAME acceleration to photon and baryon velocities only after the
# standard TCA/non-TCA equations have both been assembled. This leaves the
# standard leading slip closure intact by construction.
anchor=r'''    /** - ---> cdm */

    if (pba->has_cdm == _TRUE_) {'''
insert=r'''    /*
     * SDMC reciprocal terminal acoustic operator (Part VI Eq. 136).
     * Add the same acceleration to photon and baryon Euler equations.
     * During RSA the terminal profile is already far outside its support;
     * photons are not integrated there, so no terminal force is applied.
     */
    if ((ppt->kp_terminal_enabled != 0) &&
        (ppw->approx[ppw->index_ap_rsa] == (int)rsa_off)) {
      double kp_C, kp_dlnC, kp_Fk;
      sdmc_kp_terminal_terms(pba,pth,ppt,pvecback,R,k2,
                             delta_g,theta_g,theta_b,
                             &kp_C,&kp_dlnC,&kp_Fk);
      dy[pv->index_pt_theta_b] += kp_Fk;
      dy[pv->index_pt_theta_g] += kp_Fk;
    }

    /** - ---> cdm */

    if (pba->has_cdm == _TRUE_) {'''
replace_once("source/perturbations.c",anchor,insert)

print("installed SDMC reciprocal terminal acoustic operator")
