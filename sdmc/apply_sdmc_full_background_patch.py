#!/usr/bin/env python3
"""Add the branch-complete SDMC background on top of the tracker patch.

This incremental runtime patch is applied *after* apply_sdmc_tracker_patch.py.
It preserves the early-only `sdmc_tracker` control and adds a second expansion
model, `sdmc_full`, for the actual branch background

  H^2 = H_flat^2 (1 + delta_H)^2 / (1 - f_tr),

with

  delta_H(z) = A z exp(-z/tau_A) - B z^2 exp(-z/tau_B),
  f_tr(N)    = W(N) 3(1+w_b)/lambda_e^2.

The effective structural density is defined from the target H^2 itself and its
pressure is obtained from covariant background conservation.  Thus H and H'
are generated self-consistently rather than by multiplying H after the fact.

Parameter order for expansion_smg is

  Omega_X0, lambda_e, z_t, DeltaN, A, tau_A, B, tau_B.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: full-background patch already present")
        return
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one patch anchor, found {count}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: patched for sdmc_full")


# Register the full-background expansion model without altering the early-only
# sdmc_tracker control.
replace_once(
    "include/background.h",
    "enum expansion_model {lcdm, wowa, wowa_w, wede, sdmc_tracker};",
    "enum expansion_model {lcdm, wowa, wowa_w, wede, sdmc_tracker, sdmc_full};",
)

# Parse the eight branch-background parameters.
old = '''  if (strcmp(string1,"sdmc_tracker") == 0) {
    pba->expansion_model_smg = sdmc_tracker;
    flag2=_TRUE_;
    pba->parameters_size_smg = 4;
    pba->rho_evolution_smg=_FALSE_;
    class_read_list_of_doubles("expansion_smg",pba->parameters_smg,pba->parameters_size_smg);
  }

  class_test(flag2==_FALSE_,
             errmsg,
             "could not identify expansion_model value, check that it is either lcdm, wowa, wowa_w, wede, sdmc_tracker ...");'''
new = '''  if (strcmp(string1,"sdmc_tracker") == 0) {
    pba->expansion_model_smg = sdmc_tracker;
    flag2=_TRUE_;
    pba->parameters_size_smg = 4;
    pba->rho_evolution_smg=_FALSE_;
    class_read_list_of_doubles("expansion_smg",pba->parameters_smg,pba->parameters_size_smg);
  }

  if (strcmp(string1,"sdmc_full") == 0) {
    pba->expansion_model_smg = sdmc_full;
    flag2=_TRUE_;
    pba->parameters_size_smg = 8;
    pba->rho_evolution_smg=_FALSE_;
    class_read_list_of_doubles("expansion_smg",pba->parameters_smg,pba->parameters_size_smg);
  }

  class_test(flag2==_FALSE_,
             errmsg,
             "could not identify expansion_model value, check that it is either lcdm, wowa, wowa_w, wede, sdmc_tracker, sdmc_full ...");'''
replace_once("gravity_smg/gravity_models_smg.c", old, new)

# Insert the full branch background immediately after the early-only tracker.
old = '''    pvecback[pba->index_bg_rho_smg] = rho_X + rho_tr;
    pvecback[pba->index_bg_p_smg] = -rho_X + p_tr;
  }

  return _SUCCESS_;'''
new = '''    pvecback[pba->index_bg_rho_smg] = rho_X + rho_tr;
    pvecback[pba->index_bg_p_smg] = -rho_X + p_tr;
  }

  else if (pba->expansion_model_smg == sdmc_full){

    /* Branch-complete SDMC background.
     * expansion_smg = Omega_X0, lambda_e, z_t, DeltaN,
     *                 A, tau_A, B, tau_B
     *
     * First construct the late/intermediate branch history
     *
     *   H_late^2 = (rho_non + rho_X) [1+delta_H(z)]^2,
     *
     * then place the frozen early tracker on top,
     *
     *   H_target^2 = H_late^2/(1-f_tr).
     *
     * rho_smg is whatever density is required by this target H^2.  Its
     * pressure is then fixed from d rho_smg/dN + 3(rho_smg+p_smg)=0,
     * guaranteeing a background-consistent H'.
     */
    double Omega_X0 = pba->parameters_smg[0];
    double lambda_e = pba->parameters_smg[1];
    double z_t = pba->parameters_smg[2];
    double dN = pba->parameters_smg[3];
    double A_late = pba->parameters_smg[4];
    double tau_A = pba->parameters_smg[5];
    double B_late = pba->parameters_smg[6];
    double tau_B = pba->parameters_smg[7];

    class_test(lambda_e <= 0., pba->error_message,
               "sdmc_full requires lambda_e > 0");
    class_test(dN <= 0., pba->error_message,
               "sdmc_full requires DeltaN > 0");
    class_test(tau_A <= 0. || tau_B <= 0., pba->error_message,
               "sdmc_full requires positive late-profile decay scales");

    double z = 1./a - 1.;
    double N = log(a);
    double rho_X = Omega_X0*pow(pba->H0,2);
    double rho_flat = rho_tot + rho_X;
    double rho_non_N = -3.*(rho_tot + p_tot);
    double rho_flat_N = rho_non_N;

    double eA = exp(-z/tau_A);
    double eB = exp(-z/tau_B);
    double delta = A_late*z*eA - B_late*z*z*eB;
    double ddelta_dz = A_late*eA*(1.-z/tau_A)
                       - B_late*eB*(2.*z-z*z/tau_B);
    double delta_N = -(1.+z)*ddelta_dz;
    double g = 1.+delta;

    class_test(g <= 0., pba->error_message,
               "sdmc_full reached 1+delta_H <= 0 at z=%g",z);

    double rho_late = rho_flat*g*g;
    double rho_late_N = rho_flat_N*g*g + 2.*rho_flat*g*delta_N;

    /* Frozen lambda_e tracker and smooth handoff. */
    double Nt = -log(1.+z_t);
    double u = (N-Nt)/dN;
    double th = tanh(u);
    double W = 0.5*(1.-th);
    double W_N = -0.5*(1.-th*th)/dN;

    double wb = p_tot/rho_tot;
    double x_r = 3.*wb;
    if (x_r < 0.) x_r = 0.;
    if (x_r > 1.) x_r = 1.;
    double f0 = (3.+x_r)/(lambda_e*lambda_e);
    double f0_N = -x_r*(1.-x_r)/(lambda_e*lambda_e);
    double f = W*f0;
    double f_N = W_N*f0 + W*f0_N;

    class_test(f >= 1., pba->error_message,
               "sdmc_full tracker fraction reached f >= 1");

    double one_minus_f = 1.-f;
    double rho_target = rho_late/one_minus_f;
    double rho_target_N = rho_late_N/one_minus_f
                        + rho_late*f_N/(one_minus_f*one_minus_f);

    double rho_smg = rho_target-rho_tot;
    double rho_smg_N = rho_target_N-rho_non_N;
    double p_smg = -rho_smg-rho_smg_N/3.;

    pvecback[pba->index_bg_rho_smg] = rho_smg;
    pvecback[pba->index_bg_p_smg] = p_smg;
  }

  return _SUCCESS_;'''
replace_once("gravity_smg/gravity_models_smg.c", old, new)

print("SDMC full late+tracker background patch complete.")
