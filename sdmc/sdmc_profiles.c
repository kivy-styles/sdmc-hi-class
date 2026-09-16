#include "sdmc_profiles.h"
#include <math.h>

static const double SDMC_H0 = 70.8514;
static const double SDMC_OMEGA_R0 = 8.3268e-5;
static const double SDMC_OMEGA_B_H2 = 0.02239952;

/* Frozen Kp optimized exact-No-Slip release. */
static const double SDMC_KP_NOSLIP_AM = 0.02586;
static const double SDMC_KP_NOSLIP_ZC = 4.38;
static const double SDMC_KP_NOSLIP_W = 0.4467;

sdmc_background_params sdmc_params(sdmc_branch branch) {
  sdmc_background_params p;
  p.H0_km_s_Mpc = SDMC_H0;
  p.Omega_r0 = SDMC_OMEGA_R0;
  p.omega_b = SDMC_OMEGA_B_H2;
  p.lambda_e = 20.0;
  p.cphi2_e = 0.003;
  p.tracker_dN = 0.5;

  if (branch == SDMC_BRANCH_KP) {
    p.Omega_m0 = 0.2836698;
    p.r_s_Mpc = 143.779;
    p.r_d_Mpc = 145.2577;
    p.tracker_zt = 21.26;
  }
  else {
    p.Omega_m0 = 0.2925181427;
    p.r_s_Mpc = 142.688139;
    p.r_d_Mpc = 145.395978;
    p.tracker_zt = 21.10;
  }

  return p;
}

double sdmc_deltaH(double z, sdmc_branch branch) {
  if (branch == SDMC_BRANCH_KP) {
    return 0.874 * (
      0.024283 * z * exp(-z/0.25)
      - 0.0132767 * z * z * exp(-z/1.5)
    );
  }

  return 0.01105624999 * z * exp(-z/0.25)
       - 0.01951933685 * z * z * exp(-z/1.5);
}

double sdmc_E2_control(double z, sdmc_branch branch) {
  sdmc_background_params p = sdmc_params(branch);
  const double zp1 = 1.0 + z;
  const double Omega_x0 = 1.0 - p.Omega_m0 - p.Omega_r0;
  return p.Omega_m0 * zp1*zp1*zp1
       + p.Omega_r0 * zp1*zp1*zp1*zp1
       + Omega_x0;
}

double sdmc_E_late(double z, sdmc_branch branch) {
  return sqrt(sdmc_E2_control(z, branch)) * (1.0 + sdmc_deltaH(z, branch));
}

double sdmc_tracker_window_N(double N, sdmc_branch branch) {
  const sdmc_background_params p = sdmc_params(branch);
  const double Nt = -log(1.0 + p.tracker_zt);
  return 0.5 * (1.0 - tanh((N - Nt)/p.tracker_dN));
}

double sdmc_tracker_fraction(double z, sdmc_branch branch) {
  const sdmc_background_params p = sdmc_params(branch);
  const double zp1 = 1.0 + z;
  const double rho_m = p.Omega_m0 * zp1*zp1*zp1;
  const double rho_r = p.Omega_r0 * zp1*zp1*zp1*zp1;
  const double wb = (rho_r/3.0)/(rho_m + rho_r);
  return 3.0 * (1.0 + wb)/(p.lambda_e*p.lambda_e);
}

double sdmc_E_with_tracker(double z, sdmc_branch branch) {
  const double N = -log(1.0 + z);
  const double W = sdmc_tracker_window_N(N, branch);
  const double ftrk = sdmc_tracker_fraction(z, branch);
  const double f = W * ftrk;
  return sdmc_E_late(z, branch)/sqrt(1.0 - f);
}

static double sdmc_kp_release_S(double N) {
  const double Nc = -log(1.0 + SDMC_KP_NOSLIP_ZC);
  const double x = (N - Nc)/(2.0*SDMC_KP_NOSLIP_W);
  return 0.5*(1.0 + tanh(x));
}

double sdmc_kp_noslip_M2(double N, sdmc_branch branch) {
  if (branch != SDMC_BRANCH_KP) {
    return 1.0;
  }
  return exp(SDMC_KP_NOSLIP_AM * sdmc_kp_release_S(N));
}

double sdmc_kp_noslip_alphaM(double N, sdmc_branch branch) {
  double Nc, x, ch;
  if (branch != SDMC_BRANCH_KP) {
    return 0.0;
  }
  Nc = -log(1.0 + SDMC_KP_NOSLIP_ZC);
  x = (N - Nc)/(2.0*SDMC_KP_NOSLIP_W);
  ch = cosh(x);
  return SDMC_KP_NOSLIP_AM/(4.0*SDMC_KP_NOSLIP_W*ch*ch);
}

double sdmc_kp_noslip_alphaB(double N, sdmc_branch branch) {
  return -2.0 * sdmc_kp_noslip_alphaM(N, branch);
}

double sdmc_kp_terminal_C(double z, sdmc_branch branch,
                          double C_b, double z_on, double dz) {
  if (branch != SDMC_BRANCH_KP) {
    return 1.0;
  }

  /* W_K -> 1 below the upper edge and -> 0 above it. */
  const double Wk = 0.5 * (1.0 - tanh((z - z_on)/dz));
  return 1.0 - (1.0 - C_b) * Wk;
}
