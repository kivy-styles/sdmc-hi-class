from pathlib import Path


def replace_once(path, old, new):
    p = Path(path)
    text = p.read_text()
    if new in text:
        return
    if text.count(old) != 1:
        raise RuntimeError(f"{path}: expected one anchor, got {text.count(old)}")
    p.write_text(text.replace(old, new))


replace_once(
    "include/background.h",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip,\n",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip,\n",
)

parser_anchor = '''  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
parser_new = '''  if (strcmp(string1,"sdmc_v3_native_noslip") == 0) {
     pba->gravity_model_smg = sdmc_v3_native_noslip;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 5;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||
                pba->parameters_2_smg[2] <= 0. || pba->parameters_2_smg[3] <= 0. ||
                pba->parameters_2_smg[4] <= 0., errmsg,
                "sdmc_v3_native_noslip requires positive A_F, zc, width, c_s^2 and D floor");
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
replace_once("gravity_smg/gravity_models_smg.c", parser_anchor, parser_new)

closure_anchor = '''  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
closure_new = r'''  else if (pba->gravity_model_smg == sdmc_v3_native_noslip) {

    /* NKp-v3 native-exact No-Slip family. The profile is judged directly
     * by hi_class's native Horndeski scalar-gradient numerator after the
     * split-sector bookkeeping, rather than by the older B0 designer ODE. */
    double AF = pba->parameters_2_smg[0];
    double zc = pba->parameters_2_smg[1];
    double width = pba->parameters_2_smg[2];
    double cs2_target = pba->parameters_2_smg[3];
    double Dfloor = pba->parameters_2_smg[4];
    double N = log(a);
    double Nc = -log(1.+zc);
    double x = (N-Nc)/(2.*width);
    double th = tanh(x);
    double ch = cosh(x);
    double sech2 = 1./(ch*ch);
    double S = .5*(1.+th);

    /* F=M_*^2, alpha_M=d ln F/dN, alpha_B=-2 alpha_M exactly. */
    double F = exp(AF*S);
    double am = AF/(4.*width)*sech2;
    double amp = -AF/(4.*width*width)*sech2*th;
    double bra = -2.*am;
    double bra_N = -2.*amp;

    double H = pvecback[pba->index_bg_H];
    double rho_o = pvecback[pba->index_bg_rho_tot_wo_smg];
    double p_o = pvecback[pba->index_bg_p_tot_wo_smg];
    double rho_s = pvecback[pba->index_bg_rho_smg];
    double p_s = pvecback[pba->index_bg_p_smg];

    /* Exact native hi_class Horndeski numerator for alpha_H=alpha_T=0. */
    double cs2num =
        .5*(2.-bra)*(bra+2.*am)
      + 1.5*(2.-bra)*(rho_s+p_s)/(H*H)
      - 1.5*(2.-2.*F+bra*F)/F*(rho_o+p_o)/(H*H)
      + bra_N;

    /* D changes only the denominator. Positive native numerator gives the
     * requested sound speed; the floor regularises the spectator limit but
     * is deliberately not allowed to hide a negative numerator. */
    double Dnative = cs2num/cs2_target;
    double D = Dnative;
    if (D < Dfloor) D = Dfloor;
    double ak = D-1.5*bra*bra;

    pvecback[pba->index_bg_kineticity_smg] = ak;
    pvecback[pba->index_bg_braiding_smg] = bra;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = am;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = F-1.;
    pvecback[pba->index_bg_M2_smg] = F;
  }

  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
replace_once("gravity_smg/gravity_models_smg.c", closure_anchor, closure_new)

print("SDMC NKp-v3 native-exact No-Slip patch complete.")
