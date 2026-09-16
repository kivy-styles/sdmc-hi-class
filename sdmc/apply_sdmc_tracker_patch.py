#!/usr/bin/env python3
"""Patch hi_class with the frozen SDMC lambda_e tracker expansion model.

The patch is deliberately kept as a small, auditable layer while the SDMC
Boltzmann implementation is being validated.  It adds one non-dynamical
expansion model, `sdmc_tracker`, leaving all existing hi_class models unchanged.

For a matter+radiation background the frozen tracker fraction is

    f_tr(N) = W(N) * 3(1+w_b)/lambda_e^2,

with W the already adopted late handoff.  The Hubble correction is implemented
as the equivalent separately conserved tracker density on top of a constant
late structural density.  The tracker pressure is obtained from d rho_tr/dN,
so H and H' remain mutually consistent rather than inserting a bare H rescale.
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

# 2) Parse expansion_smg = Omega_X0, lambda_e, z_t, DeltaN.
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

# 3) Add the homogeneous tracker density and the pressure required by
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

print("SDMC tracker source patch complete.")
