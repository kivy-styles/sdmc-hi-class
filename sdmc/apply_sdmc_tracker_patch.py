#!/usr/bin/env python3
"""Patch hi_class with the frozen SDMC early tracker and designer clustering.

The patch is deliberately kept as a small, auditable layer while the SDMC
Boltzmann implementation is being validated.  It adds:

1. `sdmc_tracker`, a non-dynamical expansion model implementing the frozen
   lambda_e tracker with the adopted late handoff; and
2. `sdmc_early_cs`, a designer EFT perturbation closure that keeps
   alpha_B=alpha_M=alpha_T=alpha_H=0 and chooses alpha_K so that hi_class's
   scalar propagation speed equals a requested constant c_phi^2.

For the tracker background

    f_tr(N) = W(N) * 3(1+w_b)/lambda_e^2,

with W the already adopted late handoff.  The Hubble correction is implemented
as the equivalent separately conserved tracker density on top of a constant
late structural density.  The tracker pressure is obtained from d rho_tr/dN,
so H and H' remain mutually consistent rather than inserting a bare H rescale.

For the designer early clustering closure, with M_*^2=1 and vanishing
braiding/running/tensor excess, hi_class reduces to

    c_s^2 = [3 (rho_smg+p_smg)/H^2] / alpha_K.

Hence choosing

    alpha_K = 3 (rho_smg+p_smg) / (H^2 c_phi^2)

realizes the requested sound speed directly.  The constant late structural
piece cancels from rho_smg+p_smg, so the kineticity follows the dynamical
tracker rather than becoming large at late times.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: patch already present")
        return
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one patch anchor, found {count}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: patched")


# 1) Register the expansion model.
replace_once(
    "include/background.h",
    "enum expansion_model {lcdm, wowa, wowa_w, wede};",
    "enum expansion_model {lcdm, wowa, wowa_w, wede, sdmc_tracker};",
)

# 2) Register the designer early-clustering gravity model.
replace_once(
    "include/background.h",
    "enum gravity_model {propto_omega, propto_omega_bh, propto_scale,\n    constant_alphas,\n",
    "enum gravity_model {propto_omega, propto_omega_bh, propto_scale,\n    constant_alphas, sdmc_early_cs,\n",
)

# 3) Parse the designer sound-speed target from parameters_smg.
anchor = '''  if (strcmp(string1,"constant_alphas") == 0) {
     pba->gravity_model_smg = constant_alphas;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 5;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
replacement = '''  if (strcmp(string1,"constant_alphas") == 0) {
     pba->gravity_model_smg = constant_alphas;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 5;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

  if (strcmp(string1,"sdmc_early_cs") == 0) {
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
replace_once("gravity_smg/gravity_models_smg.c", anchor, replacement)

# 4) Parse expansion_smg = Omega_X0, lambda_e, z_t, DeltaN.
anchor = '''  if (strcmp(string1,"wede") == 0) {
    //ILSWEDE
    pba->expansion_model_smg = wede;
    flag2=_TRUE_;
    pba->parameters_size_smg = 3;
    class_read_list_of_doubles("expansion_smg",pba->parameters_smg,pba->parameters_size_smg);
    class_read_double("wede_Omega_e_regularizer_smg",pba->wede_Omega_e_regularizer_smg);
    // \t//optimize the guessing BUG: eventually leads to problem in the MCMC, perhaps the guess is too good?
    // \tif(pba->tuning_index_smg == 0){
    // \t  pba->parameters_smg[0] = pba->Omega0_smg;
    // \t}
  }

  class_test(flag2==_FALSE_,
             errmsg,
             "could not identify expansion_model value, check that it is either lcdm, wowa, wowa_w, wede ...");'''
replacement = '''  if (strcmp(string1,"wede") == 0) {
    //ILSWEDE
    pba->expansion_model_smg = wede;
    flag2=_TRUE_;
    pba->parameters_size_smg = 3;
    class_read_list_of_doubles("expansion_smg",pba->parameters_smg,pba->parameters_size_smg);
    class_read_double("wede_Omega_e_regularizer_smg",pba->wede_Omega_e_regularizer_smg);
    // \t//optimize the guessing BUG: eventually leads to problem in the MCMC, perhaps the guess is too good?
    // \tif(pba->tuning_index_smg == 0){
    // \t  pba->parameters_smg[0] = pba->Omega0_smg;
    // \t}
  }

  if (strcmp(string1,"sdmc_tracker") == 0) {
    pba->expansion_model_smg = sdmc_tracker;
    flag2=_TRUE_;
    pba->parameters_size_smg = 4;
    pba->rho_evolution_smg=_FALSE_;
    class_read_list_of_doubles("expansion_smg",pba->parameters_smg,pba->parameters_size_smg);
  }

  class_test(flag2==_FALSE_,
             errmsg,
             "could not identify expansion_model value, check that it is either lcdm, wowa, wowa_w, wede, sdmc_tracker ...");'''
replace_once("gravity_smg/gravity_models_smg.c", anchor, replacement)

# 5) Add the homogeneous tracker density and the pressure required by
#    background conservation.  rho_tot/p_tot here exclude the smg sector.
anchor = '''  else if (pba->expansion_model_smg == wede){

    //Doran-Robbers model astro-ph/0601544
    //as implemented in Pettorino et al. 1301.5279
    //NOTE: rewrite the expressions integrating the equation of state

    double Om0 = pba->parameters_smg[0];
    double w0 = pba->parameters_smg[1];
    double Ome = pba->parameters_smg[2] + pba->wede_Omega_e_regularizer_smg;

    double Om = ((Om0 - Ome*(1.-pow(a,-3.*w0)))/(Om0 + (1.-Om0)*pow(a,3*w0)) + Ome*(1.-pow(a,-3*w0)));
    double dOm_da = (3*pow(a,-1 - 3*w0)*(-1 + Om0)*(-2*pow(a,3*w0)*(-1 + Om0)*Ome + Om0*Ome + pow(a,6*w0)*(Om0 - 2*Ome + Om0*Ome))*w0)/pow(-(pow(a,3*w0)*(-1 + Om0)) + Om0,2); //from Mathematica
    //I took a_eq = a*rho_r/rho_m, with rho_r = 3*p_tot_wo_smg
    double a_eq = 3.*a*p_tot/(pvecback[pba->index_bg_rho_b]+pvecback[pba->index_bg_rho_cdm]); //tested!
    double w = a_eq/(3.*(a+a_eq)) -a/(3.*(1-Om)*Om)*dOm_da;

    pvecback[pba->index_bg_rho_smg] = rho_tot*Om/(1.-Om);
    //pow(pba->H0,2)/pow(a,3)*Om*(Om-1.)/(Om0-1.)*(1.+a_eq/a)/(1.+a_eq); //this eq is from Pettorino et al, not working
    pvecback[pba->index_bg_p_smg] = w*pvecback[pba->index_bg_rho_smg];

//       if (a>0.9)
// \tprintf("a = %e, w = %f, Om_de = %e, rho_de/rho_t = %e \\n",a,w,Om,
// \t       pvecback[pba->index_bg_rho_smg]/(pvecback[pba->index_bg_rho_smg]+rho_tot));
  }

  return _SUCCESS_;'''
replacement = '''  else if (pba->expansion_model_smg == wede){

    //Doran-Robbers model astro-ph/0601544
    //as implemented in Pettorino et al. 1301.5279
    //NOTE: rewrite the expressions integrating the equation of state

    double Om0 = pba->parameters_smg[0];
    double w0 = pba->parameters_smg[1];
    double Ome = pba->parameters_smg[2] + pba->wede_Omega_e_regularizer_smg;

    double Om = ((Om0 - Ome*(1.-pow(a,-3.*w0)))/(Om0 + (1.-Om0)*pow(a,3*w0)) + Ome*(1.-pow(a,-3*w0)));
    double dOm_da = (3*pow(a,-1 - 3*w0)*(-1 + Om0)*(-2*pow(a,3*w0)*(-1 + Om0)*Ome + Om0*Ome + pow(a,6*w0)*(Om0 - 2*Ome + Om0*Ome))*w0)/pow(-(pow(a,3*w0)*(-1 + Om0)) + Om0,2); //from Mathematica
    //I took a_eq = a*rho_r/rho_m, with rho_r = 3*p_tot_wo_smg
    double a_eq = 3.*a*p_tot/(pvecback[pba->index_bg_rho_b]+pvecback[pba->index_bg_rho_cdm]); //tested!
    double w = a_eq/(3.*(a+a_eq)) -a/(3.*(1-Om)*Om)*dOm_da;

    pvecback[pba->index_bg_rho_smg] = rho_tot*Om/(1.-Om);
    //pow(pba->H0,2)/pow(a,3)*Om*(Om-1.)/(Om0-1.)*(1.+a_eq/a)/(1.+a_eq); //this eq is from Pettorino et al, not working
    pvecback[pba->index_bg_p_smg] = w*pvecback[pba->index_bg_rho_smg];

//       if (a>0.9)
// \tprintf("a = %e, w = %f, Om_de = %e, rho_de/rho_t = %e \\n",a,w,Om,
// \t       pvecback[pba->index_bg_rho_smg]/(pvecback[pba->index_bg_rho_smg]+rho_tot));
  }

  else if (pba->expansion_model_smg == sdmc_tracker){

    /* Frozen SDMC lambda_e tracker plus a constant late structural density.
     * expansion_smg = Omega_X0, lambda_e, z_t, DeltaN
     *
     * The tracker fraction is the manuscript attractor
     *   f0 = 3(1+w_b)/lambda_e^2,
     * multiplied by an early-on/late-off tanh handoff W(N).  Defining
     * rho_tr = f/(1-f) rho_base makes H/H_base = 1/sqrt(1-f) exactly.
     * The pressure below follows from d rho_tr/dN + 3(rho_tr+p_tr)=0.
     */
    double Omega_X0 = pba->parameters_smg[0];
    double lambda_e = pba->parameters_smg[1];
    double z_t = pba->parameters_smg[2];
    double dN = pba->parameters_smg[3];
    double rho_X = Omega_X0*pow(pba->H0,2);
    double rho_base = rho_tot + rho_X;
    double p_base = p_tot - rho_X;
    double N = log(a);
    double Nt = -log(1.+z_t);
    double u = (N-Nt)/dN;
    double th = tanh(u);
    double W = 0.5*(1.-th);
    double sech2 = 1.-th*th;
    double W_N = -0.5*sech2/dN;

    /* For the frozen runs used here the non-smg sector is matter+radiation.
     * x_r=3 w_b is the radiation share of that sector and x_r'=-x_r(1-x_r).
     */
    double wb = p_tot/rho_tot;
    double x_r = 3.*wb;
    if (x_r < 0.) x_r = 0.;
    if (x_r > 1.) x_r = 1.;
    double f0 = (3.+x_r)/(lambda_e*lambda_e);
    double f0_N = -x_r*(1.-x_r)/(lambda_e*lambda_e);
    double f = W*f0;
    double f_N = W_N*f0 + W*f0_N;

    class_test(lambda_e <= 0., pba->error_message,
               "sdmc_tracker requires lambda_e > 0");
    class_test(dN <= 0., pba->error_message,
               "sdmc_tracker requires DeltaN > 0");
    class_test(f >= 1., pba->error_message,
               "sdmc_tracker fraction reached f >= 1");

    double rho_tr = 0.;
    double p_tr = 0.;
    if (f > 1.e-14) {
      double q = f/(1.-f);
      double dlnq_dN = f_N/(f*(1.-f));
      double w_base = p_base/rho_base;
      double w_tr = w_base - dlnq_dN/3.;
      rho_tr = q*rho_base;
      p_tr = w_tr*rho_tr;
    }

    pvecback[pba->index_bg_rho_smg] = rho_X + rho_tr;
    pvecback[pba->index_bg_p_smg] = -rho_X + p_tr;
  }

  return _SUCCESS_;'''
replace_once("gravity_smg/gravity_models_smg.c", anchor, replacement)

# 6) Add the designer alpha_K closure.  In this no-braiding, constant-M2
#    limit hi_class has cs2num = 3(rho_smg+p_smg)/H^2 and D=alpha_K.
anchor = '''  else if (pba->gravity_model_smg == constant_alphas) {

    double c_k = pba->parameters_2_smg[0];
    double c_b = pba->parameters_2_smg[1];
    double c_m = pba->parameters_2_smg[2];
    double c_t = pba->parameters_2_smg[3];

    pvecback[pba->index_bg_kineticity_smg] = c_k;
    pvecback[pba->index_bg_braiding_smg] = c_b;
    pvecback[pba->index_bg_tensor_excess_smg] = c_t;
    pvecback[pba->index_bg_M2_running_smg] = c_m;
    pvecback[pba->index_bg_delta_M2_smg] = delta_M2; //M2-1
    pvecback[pba->index_bg_M2_smg] = 1.+delta_M2;
  }

  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
replacement = '''  else if (pba->gravity_model_smg == constant_alphas) {

    double c_k = pba->parameters_2_smg[0];
    double c_b = pba->parameters_2_smg[1];
    double c_m = pba->parameters_2_smg[2];
    double c_t = pba->parameters_2_smg[3];

    pvecback[pba->index_bg_kineticity_smg] = c_k;
    pvecback[pba->index_bg_braiding_smg] = c_b;
    pvecback[pba->index_bg_tensor_excess_smg] = c_t;
    pvecback[pba->index_bg_M2_running_smg] = c_m;
    pvecback[pba->index_bg_delta_M2_smg] = delta_M2; //M2-1
    pvecback[pba->index_bg_M2_smg] = 1.+delta_M2;
  }

  else if (pba->gravity_model_smg == sdmc_early_cs) {

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
replace_once("gravity_smg/gravity_models_smg.c", anchor, replacement)

print("SDMC tracker + designer early-clustering source patch complete.")
