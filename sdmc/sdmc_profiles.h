#ifndef SDMC_PROFILES_H
#define SDMC_PROFILES_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  SDMC_BRANCH_KP = 0,
  SDMC_BRANCH_NKP_V2 = 1
} sdmc_branch;

typedef struct {
  double H0_km_s_Mpc;
  double Omega_m0;
  double Omega_r0;
  double omega_b;
  double r_s_Mpc;
  double r_d_Mpc;
  double lambda_e;
  double cphi2_e;
  double tracker_zt;
  double tracker_dN;
} sdmc_background_params;

/* Frozen parameter registry. */
sdmc_background_params sdmc_params(sdmc_branch branch);

/* Late phenomenological correction defined by H/H_control = 1 + deltaH. */
double sdmc_deltaH(double z, sdmc_branch branch);

/* Flat matter+radiation+constant-X control E(z)^2. */
double sdmc_E2_control(double z, sdmc_branch branch);

/* Current background diagnostic before a fully dynamical scalar rewrite. */
double sdmc_E_late(double z, sdmc_branch branch);

/* Smooth early-tracker handoff: 1 at early times, 0 at late times. */
double sdmc_tracker_window_N(double N, sdmc_branch branch);

/* Exponential-tracker fractional density diagnostic 3(1+w_b)/lambda_e^2. */
double sdmc_tracker_fraction(double z, sdmc_branch branch);

/* Tracker-corrected background diagnostic used by the closure tests. */
double sdmc_E_with_tracker(double z, sdmc_branch branch);

/*
 * Frozen Kp exact-No-Slip late profile from Part VI.
 * N = ln(a), with a0=1.  The monotonic saturated release is
 *   ln M_*^2 = A_M S(N),
 *   S = [1+tanh((N-Nc)/(2w))]/2,
 * where A_M=0.02586, zc=4.38, w=0.4467.
 * This reproduces M_*0^2/M_Pl^2 ~=1.02559 and alpha_M,max~=0.01447.
 * NKp-v2 is deliberately not assigned this Kp profile.
 */
double sdmc_kp_noslip_M2(double N, sdmc_branch branch);
double sdmc_kp_noslip_alphaM(double N, sdmc_branch branch);
double sdmc_kp_noslip_alphaB(double N, sdmc_branch branch);

/* Kp terminal reciprocal-acoustic coefficient.  For NKp-v2 this is identically 1. */
double sdmc_kp_terminal_C(double z, sdmc_branch branch,
                          double C_b, double z_on, double dz);

#ifdef __cplusplus
}
#endif

#endif
