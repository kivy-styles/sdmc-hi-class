#include "sdmc_profiles.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

static int approx(double x, double y, double tol) {
  return fabs(x-y) <= tol;
}

int main(void) {
  sdmc_background_params kp = sdmc_params(SDMC_BRANCH_KP);
  sdmc_background_params v2 = sdmc_params(SDMC_BRANCH_NKP_V2);

  assert(approx(kp.H0_km_s_Mpc, 70.8514, 1e-12));
  assert(approx(kp.Omega_m0, 0.2836698, 1e-12));
  assert(approx(v2.Omega_m0, 0.2925181427, 1e-12));
  assert(approx(kp.r_d_Mpc, 145.2577, 1e-10));
  assert(approx(v2.r_d_Mpc, 145.395978, 1e-10));

  /* By construction both late corrections vanish at z=0. */
  assert(fabs(sdmc_deltaH(0.0, SDMC_BRANCH_KP)) < 1e-15);
  assert(fabs(sdmc_deltaH(0.0, SDMC_BRANCH_NKP_V2)) < 1e-15);
  assert(approx(sdmc_E_late(0.0, SDMC_BRANCH_KP), 1.0, 1e-12));
  assert(approx(sdmc_E_late(0.0, SDMC_BRANCH_NKP_V2), 1.0, 1e-12));

  /* Tracker window should be early-on and late-off. */
  assert(sdmc_tracker_window_N(-8.0, SDMC_BRANCH_KP) > 0.999);
  assert(sdmc_tracker_window_N(0.0, SDMC_BRANCH_KP) < 1e-4);
  assert(sdmc_tracker_window_N(-8.0, SDMC_BRANCH_NKP_V2) > 0.999);
  assert(sdmc_tracker_window_N(0.0, SDMC_BRANCH_NKP_V2) < 1e-4);

  /* The tracker fraction should approach the radiation and matter attractor values. */
  {
    double frad = sdmc_tracker_fraction(1e8, SDMC_BRANCH_KP);
    double fmat = sdmc_tracker_fraction(100.0, SDMC_BRANCH_KP);
    assert(frad > 0.0095 && frad < 0.0101); /* 4/lambda^2 = 0.01 */
    assert(fmat > 0.0074 && fmat < 0.0082); /* near 3/lambda^2 = 0.0075 */
  }

  /* Frozen Kp exact-No-Slip profile reproduces the Part VI benchmark. */
  {
    const double Nc = -log(1.0 + 4.38);
    const double M20 = sdmc_kp_noslip_M2(0.0, SDMC_BRANCH_KP);
    const double aMpeak = sdmc_kp_noslip_alphaM(Nc, SDMC_BRANCH_KP);
    const double aBpeak = sdmc_kp_noslip_alphaB(Nc, SDMC_BRANCH_KP);
    assert(approx(M20, 1.02559769, 2e-8));
    assert(approx(aMpeak, 0.01447280, 2e-8));
    assert(approx(aBpeak, -0.02894560, 4e-8));
    assert(sdmc_kp_noslip_M2(-10.0, SDMC_BRANCH_KP) < 1.000001);
  }

  /* Do not silently assign the Kp No-Slip profile to NKp-v2. */
  assert(approx(sdmc_kp_noslip_M2(0.0, SDMC_BRANCH_NKP_V2), 1.0, 1e-15));
  assert(approx(sdmc_kp_noslip_alphaM(0.0, SDMC_BRANCH_NKP_V2), 0.0, 1e-15));
  assert(approx(sdmc_kp_noslip_alphaB(0.0, SDMC_BRANCH_NKP_V2), 0.0, 1e-15));

  /* NKp-v2 has no terminal reciprocal-acoustic operator. */
  assert(approx(sdmc_kp_terminal_C(1090.0, SDMC_BRANCH_NKP_V2,
                                   0.539812, 1088.68, 0.5), 1.0, 1e-15));

  /* Kp operator remains bounded between C_b and unity. */
  {
    const double Cb = 0.539812;
    double C1 = sdmc_kp_terminal_C(1000.0, SDMC_BRANCH_KP, Cb, 1088.68, 0.5);
    double C2 = sdmc_kp_terminal_C(1200.0, SDMC_BRANCH_KP, Cb, 1088.68, 0.5);
    assert(C1 >= Cb && C1 <= 1.0);
    assert(C2 >= Cb && C2 <= 1.0);
    assert(C1 < C2);
  }

  printf("SDMC profile smoke test passed.\n");
  printf("Kp: Omega_m0=%.10f rs=%.6f rd=%.6f zt=%.3f M2_0=%.8f alphaM_peak=%.8f\n",
         kp.Omega_m0, kp.r_s_Mpc, kp.r_d_Mpc, kp.tracker_zt,
         sdmc_kp_noslip_M2(0.0, SDMC_BRANCH_KP),
         sdmc_kp_noslip_alphaM(-log(1.0+4.38), SDMC_BRANCH_KP));
  printf("v2: Omega_m0=%.10f rs=%.6f rd=%.6f zt=%.3f\n",
         v2.Omega_m0, v2.r_s_Mpc, v2.r_d_Mpc, v2.tracker_zt);
  return 0;
}
