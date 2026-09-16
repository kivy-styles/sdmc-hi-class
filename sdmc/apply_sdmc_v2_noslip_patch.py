#!/usr/bin/env python3
"""Add the NKp-v2 B0 No-Slip perturbation closure to hi_class.

This incremental runtime patch is applied after the tracker/clustering and
branch-complete background patches.  It adds a frozen-branch gravity model,
`sdmc_v2_noslip`, that reconstructs the endpoint-clean B0 profile by solving
its shooting ODE directly in N=ln(a), rather than fitting an arbitrary pulse.

The profile obeys

  F' = alpha_M F,
  N_s = -2(1+alpha_M)(h+alpha_M) - 2 alpha_M'
        - (3 Omega_m + 4 Omega_r)/F,

with constant positive N_s in the late B0 sector and exact No-Slip

  alpha_B = -2 alpha_M,  alpha_T = alpha_H = 0.

The start redshift is supplied as a parameter.  For the current exact hi_class
background, z_start=26.52404317 is the endpoint-clean re-shoot of the recovered
N_s=1e-3 B0 family; it gives F0 ~= 1.027893 and alpha_M(0) ~= 0.  The recovered
original table used z_start=26.499317 and F0=1.0278926, so the re-shoot is only
a numerical background-coordinate reconciliation, not a new physical fit.

To connect the early clustered scalar to the late No-Slip scalar without a
new timing parameter, the scalar sound speed is blended with the already
frozen tracker handoff W(N):

  c_s^2(N) = W c_e^2 + (1-W) c_l^2,
  D = N_s(actual)/c_s^2,
  alpha_K = D - 3 alpha_B^2/2.

Thus the acoustic era retains c_e^2=0.003, while the late B0 closure tends to
c_l^2=0.1 and D=0.01.  Stability is checked by native hi_class diagnostics.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: v2 No-Slip patch already present")
        return
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{path}: expected one patch anchor, found {n}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: patched for sdmc_v2_noslip")


# Register the gravity model.
replace_once(
    "include/background.h",
    "    constant_alphas, sdmc_early_cs,\n",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip,\n",
)

# Parse: N_s target, c_e^2, c_l^2, z_start, z_handoff, DeltaN_handoff.
old = '''  if (strcmp(string1,"sdmc_early_cs") == 0) {
     pba->gravity_model_smg = sdmc_early_cs;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 1;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0., errmsg,
                "sdmc_early_cs requires c_phi^2 > 0");
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
new = '''  if (strcmp(string1,"sdmc_early_cs") == 0) {
     pba->gravity_model_smg = sdmc_early_cs;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 1;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0., errmsg,
                "sdmc_early_cs requires c_phi^2 > 0");
   }

  if (strcmp(string1,"sdmc_v2_noslip") == 0) {
     pba->gravity_model_smg = sdmc_v2_noslip;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 6;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0., errmsg,
                "sdmc_v2_noslip requires N_s > 0");
     class_test(pba->parameters_2_smg[1] <= 0. || pba->parameters_2_smg[2] <= 0., errmsg,
                "sdmc_v2_noslip requires positive early and late c_s^2");
     class_test(pba->parameters_2_smg[3] <= 0. || pba->parameters_2_smg[4] <= 0. || pba->parameters_2_smg[5] <= 0., errmsg,
                "sdmc_v2_noslip requires positive start/handoff redshifts and handoff width");
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
replace_once("gravity_smg/gravity_models_smg.c", old, new)

# Insert compact ODE helpers before gravity_models_get_alphas_par_smg.
anchor = '''int gravity_models_get_alphas_par_smg(
                                      struct background *pba,
                                      double a,
                                      double * pvecback,
                                      double * pvecback_B
                                      ) {'''
helper = r'''/* -------------------------------------------------------------------------
 * SDMC NKp-v2 endpoint-clean B0 profile helpers.
 * ------------------------------------------------------------------------- */
static void sdmc_v2_background_for_b0(struct background *pba,
                                      double N,
                                      double *h,
                                      double *Om,
                                      double *Or) {
  double a = exp(N);
  double z = 1./a-1.;

  /* Matter/radiation coordinates used by the frozen NKp-v2 branch. */
  double Om0 = pba->Omega0_b + pba->Omega0_cdm;
  double Or0 = pba->Omega0_g + pba->Omega0_ur;
  double rm = Om0*pow(a,-3.);
  double rr = Or0*pow(a,-4.);

  /* The B0 closure is defined on sdmc_full. */
  double Ox = pba->parameters_smg[0];
  double lambda_e = pba->parameters_smg[1];
  double zt = pba->parameters_smg[2];
  double dNt = pba->parameters_smg[3];
  double A = pba->parameters_smg[4];
  double tauA = pba->parameters_smg[5];
  double B = pba->parameters_smg[6];
  double tauB = pba->parameters_smg[7];

  double flat = rm+rr+Ox;
  double flat_N = -3.*rm-4.*rr;

  double eA = exp(-z/tauA);
  double eB = exp(-z/tauB);
  double delta = A*z*eA-B*z*z*eB;
  double ddelta_dz = A*eA*(1.-z/tauA)-B*eB*(2.*z-z*z/tauB);
  double delta_N = -(1.+z)*ddelta_dz;

  double xr = rr/(rm+rr);
  double f0 = (3.+xr)/(lambda_e*lambda_e);
  double f0_N = -xr*(1.-xr)/(lambda_e*lambda_e);
  double Nt = -log(1.+zt);
  double u = (N-Nt)/dNt;
  double th = tanh(u);
  double W = .5*(1.-th);
  double W_N = -.5*(1.-th*th)/dNt;
  double f = W*f0;
  double f_N = W_N*f0+W*f0_N;

  double g = 1.+delta;
  double E2 = flat*g*g/(1.-f);
  double lnE2_N = flat_N/flat + 2.*delta_N/g + f_N/(1.-f);

  *h = .5*lnE2_N;
  *Om = rm/E2;
  *Or = rr/E2;
}

static void sdmc_v2_b0_rhs(struct background *pba,
                           double N,
                           double Ns,
                           double F,
                           double am,
                           double *Fp,
                           double *amp) {
  double h, Om, Or;
  sdmc_v2_background_for_b0(pba,N,&h,&Om,&Or);
  *Fp = am*F;
  *amp = -.5*(Ns + 2.*(1.+am)*(h+am) + (3.*Om+4.*Or)/F);
}

static void sdmc_v2_b0_profile(struct background *pba,
                               double N,
                               double Ns,
                               double zstart,
                               double *F,
                               double *am,
                               double *amp) {
  double N0 = -log(1.+zstart);
  if (N <= N0) {
    *F = 1.;
    *am = 0.;
    *amp = 0.;
    return;
  }
  if (N > 0.) N = 0.;

  /* Fixed-step RK4 is deterministic and sufficiently fine for the smooth B0
   * profile.  dN<=0.002 gives sub-1e-6 agreement in alpha_M with the recovered
   * 14,000-point reference table over its active interval. */
  int nstep = (int)ceil((N-N0)/0.002);
  if (nstep < 1) nstep = 1;
  double dN = (N-N0)/nstep;
  double x=N0, f=1., m=0.;
  int i;
  for (i=0;i<nstep;i++) {
    double k1f,k1m,k2f,k2m,k3f,k3m,k4f,k4m;
    sdmc_v2_b0_rhs(pba,x,Ns,f,m,&k1f,&k1m);
    sdmc_v2_b0_rhs(pba,x+.5*dN,Ns,f+.5*dN*k1f,m+.5*dN*k1m,&k2f,&k2m);
    sdmc_v2_b0_rhs(pba,x+.5*dN,Ns,f+.5*dN*k2f,m+.5*dN*k2m,&k3f,&k3m);
    sdmc_v2_b0_rhs(pba,x+dN,Ns,f+dN*k3f,m+dN*k3m,&k4f,&k4m);
    f += dN*(k1f+2.*k2f+2.*k3f+k4f)/6.;
    m += dN*(k1m+2.*k2m+2.*k3m+k4m)/6.;
    x += dN;
  }
  *F=f;
  *am=m;
  {
    double dummy;
    sdmc_v2_b0_rhs(pba,N,Ns,f,m,&dummy,amp);
  }
}

int gravity_models_get_alphas_par_smg(
                                      struct background *pba,
                                      double a,
                                      double * pvecback,
                                      double * pvecback_B
                                      ) {'''
replace_once("gravity_smg/gravity_models_smg.c", anchor, helper)

# Add the No-Slip alpha closure after the early-c_s designer model.
old = '''  else if (pba->gravity_model_smg == sdmc_early_cs) {

    double cs2_target = pba->parameters_2_smg[0];
    double enthalpy_smg = pvecback[pba->index_bg_rho_smg]
                        + pvecback[pba->index_bg_p_smg];
    double H2 = rho_tot;

    /* The frozen tracker has positive enthalpy throughout its finite handoff.
     * Keep only a numerical floor far below the physical tracker signal.
     */
    if (enthalpy_smg < 1.e-30*H2)
      enthalpy_smg = 1.e-30*H2;

    pvecback[pba->index_bg_kineticity_smg] = 3.*enthalpy_smg/(H2*cs2_target);
    pvecback[pba->index_bg_braiding_smg] = 0.;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = 0.;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = 0.;
    pvecback[pba->index_bg_M2_smg] = 1.;
  }

  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
new = '''  else if (pba->gravity_model_smg == sdmc_early_cs) {

    double cs2_target = pba->parameters_2_smg[0];
    double enthalpy_smg = pvecback[pba->index_bg_rho_smg]
                        + pvecback[pba->index_bg_p_smg];
    double H2 = rho_tot;

    /* The frozen tracker has positive enthalpy throughout its finite handoff.
     * Keep only a numerical floor far below the physical tracker signal.
     */
    if (enthalpy_smg < 1.e-30*H2)
      enthalpy_smg = 1.e-30*H2;

    pvecback[pba->index_bg_kineticity_smg] = 3.*enthalpy_smg/(H2*cs2_target);
    pvecback[pba->index_bg_braiding_smg] = 0.;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = 0.;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = 0.;
    pvecback[pba->index_bg_M2_smg] = 1.;
  }

  else if (pba->gravity_model_smg == sdmc_v2_noslip) {

    double Ns_target = pba->parameters_2_smg[0];
    double cs2_early = pba->parameters_2_smg[1];
    double cs2_late = pba->parameters_2_smg[2];
    double zstart = pba->parameters_2_smg[3];
    double zhand = pba->parameters_2_smg[4];
    double dNhand = pba->parameters_2_smg[5];
    double N = log(a);
    double F, am, amp;

    sdmc_v2_b0_profile(pba,N,Ns_target,zstart,&F,&am,&amp);

    double ab = -2.*am;
    double H = pvecback[pba->index_bg_H];
    double h = pvecback[pba->index_bg_H_prime]/(a*H*H);
    double rho_m = pvecback[pba->index_bg_rho_b]
                 + pvecback[pba->index_bg_rho_cdm];
    double rho_r = pvecback[pba->index_bg_rho_g]
                 + pvecback[pba->index_bg_rho_ur];
    double Om = rho_m/(H*H);
    double Or = rho_r/(H*H);
    double Ns_actual = -2.*(1.+am)*(h+am)-2.*amp-(3.*Om+4.*Or)/F;

    double Nh = -log(1.+zhand);
    double Wh = .5*(1.-tanh((N-Nh)/dNhand));
    double cs2_target = Wh*cs2_early+(1.-Wh)*cs2_late;
    double D = Ns_actual/cs2_target;
    double ak = D-1.5*ab*ab;

    class_test(Ns_actual <= 0., pba->error_message,
               "sdmc_v2_noslip reached N_s=%g <= 0 at z=%g",Ns_actual,1./a-1.);
    class_test(D <= 0., pba->error_message,
               "sdmc_v2_noslip reached D=%g <= 0 at z=%g",D,1./a-1.);

    pvecback[pba->index_bg_kineticity_smg] = ak;
    pvecback[pba->index_bg_braiding_smg] = ab;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = am;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = F-1.;
    pvecback[pba->index_bg_M2_smg] = F;
  }

  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
replace_once("gravity_smg/gravity_models_smg.c", old, new)

print("SDMC NKp-v2 B0 No-Slip patch complete.")
