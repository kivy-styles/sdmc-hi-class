#!/usr/bin/env python3
"""Add a localized transition-bridge EFT completion to the SDMC split sector.

This patch is applied after the tracker, full-background, NKp-v2 B0, smooth
No-Slip, split-sector, and split-bookkeeping patches.

The split-sector calculation showed that the early tracker fluid can be
separated cleanly from the late Horndeski scalar, but the frozen handoff then
passes through a short interval with a negative native hi_class scalar
sound-speed numerator.  A global change of alpha_B/alpha_M can remove that
instability, but earlier scans showed that it also produces excessive growth
and lensing.

The completion introduced here is deliberately local in e-fold time:

    alpha_B = alpha_B,late + A_B P(N),

where P is a compact-like smooth pulse made from two tanh edges,

    P(N) = 1/2 [tanh((N-N_on)/w) - tanh((N-N_off)/w)].

The pulse is zero before the tracker handoff and zero again at low redshift.
It changes neither M_*^2 nor alpha_M, so the frozen late Planck-mass trajectory
is preserved exactly.  In particular, Kp keeps its frozen late No-Slip release
and NKp-v2 keeps its recovered B0 F(a); exact No-Slip is violated only inside
the bridge window when A_B != 0.

The kinetic closure uses the exact native hi_class Horndeski sound-speed
numerator.  Whenever that numerator is positive, D is chosen to realize the
requested scalar sound speed; otherwise D is held at a positive spectator
floor and native hi_class stability rejection is allowed to fail the model.
Thus the bridge scan cannot hide a gradient instability by choosing alpha_K.

This is a numerical completion test.  Bridge parameters are to be fixed by
stability and minimum disturbance, not by a Planck likelihood fit.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: transition bridge already present")
        return
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{path}: expected one bridge anchor, found {n}\n{old}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: transition bridge patched")


# Register two bridge models without changing the existing split-sector controls.
replace_once(
    "include/background.h",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip,\n",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v2_transition_bridge, sdmc_kp_transition_bridge,\n",
)

# Parse the bridge parameters immediately before the stock EFT parser.
replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''  if (strcmp(string1,"eft_alphas_power_law") == 0) {''',
    '''  if (strcmp(string1,"sdmc_v2_transition_bridge") == 0) {
     pba->gravity_model_smg = sdmc_v2_transition_bridge;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 9;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||
                pba->parameters_2_smg[2] <= 0. || pba->parameters_2_smg[3] <= 0. ||
                pba->parameters_2_smg[4] <= 0. || pba->parameters_2_smg[6] <= 0. ||
                pba->parameters_2_smg[7] <= 0. || pba->parameters_2_smg[8] <= 0. ||
                pba->parameters_2_smg[6] <= pba->parameters_2_smg[7], errmsg,
                "sdmc_v2_transition_bridge requires positive Ns, cs2, zact, activation width, D floor, z_on, z_off, bridge width and z_on>z_off");
   }

  if (strcmp(string1,"sdmc_kp_transition_bridge") == 0) {
     pba->gravity_model_smg = sdmc_kp_transition_bridge;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 6;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||
                pba->parameters_2_smg[3] <= 0. || pba->parameters_2_smg[4] <= 0. ||
                pba->parameters_2_smg[5] <= 0. || pba->parameters_2_smg[3] <= pba->parameters_2_smg[4], errmsg,
                "sdmc_kp_transition_bridge requires positive cs2, D floor, z_on, z_off, bridge width and z_on>z_off");
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
)

# Add both bridge closures immediately before the stock EFT alpha model.
replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''  else if (pba->gravity_model_smg == eft_alphas_power_law) {''',
    r'''  else if (pba->gravity_model_smg == sdmc_v2_transition_bridge) {

    /* parameters_smg = Ns, cs2, zact, dNact, Dfloor,
     *                  A_bridge, z_on, z_off, bridge_width */
    double Ns_target = pba->parameters_2_smg[0];
    double cs2_target = pba->parameters_2_smg[1];
    double zact = pba->parameters_2_smg[2];
    double dNact = pba->parameters_2_smg[3];
    double Dfloor = pba->parameters_2_smg[4];
    double Abridge = pba->parameters_2_smg[5];
    double zon = pba->parameters_2_smg[6];
    double zoff = pba->parameters_2_smg[7];
    double wbri = pba->parameters_2_smg[8];
    double N = log(a);

    double F, am, amp;
    sdmc_v2_b0_profile(pba,N,Ns_target,zact,dNact,&F,&am,&amp);

    double Non = -log(1.+zon);
    double Noff = -log(1.+zoff);
    double ton = tanh((N-Non)/wbri);
    double toff = tanh((N-Noff)/wbri);
    double pulse = .5*(ton-toff);
    double pulse_N = .5*((1.-ton*ton)-(1.-toff*toff))/wbri;

    double bra = -2.*am + Abridge*pulse;
    double bra_N = -2.*amp + Abridge*pulse_N;
    double H = pvecback[pba->index_bg_H];
    double rho_o = pvecback[pba->index_bg_rho_tot_wo_smg];
    double p_o = pvecback[pba->index_bg_p_tot_wo_smg];
    double rho_s = pvecback[pba->index_bg_rho_smg];
    double p_s = pvecback[pba->index_bg_p_smg];

    /* Exact native hi_class Horndeski c_s^2 numerator for alpha_H=alpha_T=0. */
    double cs2num =
        .5*(2.-bra)*(bra+2.*am)
      + 1.5*(2.-bra)*(rho_s+p_s)/(H*H)
      - 1.5*(2.-2.*F+bra*F)/F*(rho_o+p_o)/(H*H)
      + bra_N;

    /* alpha_K is a designer kinetic completion.  It may normalize a positive
     * numerator, but it is never allowed to mask a negative one. */
    double D = Dfloor;
    if (cs2num > 0.) {
      double Dt = cs2num/cs2_target;
      if (Dt > D) D = Dt;
    }
    double ak = D-1.5*bra*bra;

    pvecback[pba->index_bg_kineticity_smg] = ak;
    pvecback[pba->index_bg_braiding_smg] = bra;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = am;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = F-1.;
    pvecback[pba->index_bg_M2_smg] = F;
  }

  else if (pba->gravity_model_smg == sdmc_kp_transition_bridge) {

    /* parameters_smg = cs2, Dfloor, A_bridge, z_on, z_off, bridge_width */
    const double AM = 0.02586;
    const double zc = 4.38;
    const double width = 0.4467;
    double cs2_target = pba->parameters_2_smg[0];
    double Dfloor = pba->parameters_2_smg[1];
    double Abridge = pba->parameters_2_smg[2];
    double zon = pba->parameters_2_smg[3];
    double zoff = pba->parameters_2_smg[4];
    double wbri = pba->parameters_2_smg[5];
    double N = log(a);

    double Nc = -log(1.+zc);
    double x = (N-Nc)/(2.*width);
    double th = tanh(x);
    double sech2 = 1./(cosh(x)*cosh(x));
    double S = .5*(1.+th);
    double F = exp(AM*S);
    double am = AM/(4.*width)*sech2;
    double amp = -AM/(4.*width*width)*sech2*th;

    double Non = -log(1.+zon);
    double Noff = -log(1.+zoff);
    double ton = tanh((N-Non)/wbri);
    double toff = tanh((N-Noff)/wbri);
    double pulse = .5*(ton-toff);
    double pulse_N = .5*((1.-ton*ton)-(1.-toff*toff))/wbri;

    double bra = -2.*am + Abridge*pulse;
    double bra_N = -2.*amp + Abridge*pulse_N;
    double H = pvecback[pba->index_bg_H];
    double rho_o = pvecback[pba->index_bg_rho_tot_wo_smg];
    double p_o = pvecback[pba->index_bg_p_tot_wo_smg];
    double rho_s = pvecback[pba->index_bg_rho_smg];
    double p_s = pvecback[pba->index_bg_p_smg];

    double cs2num =
        .5*(2.-bra)*(bra+2.*am)
      + 1.5*(2.-bra)*(rho_s+p_s)/(H*H)
      - 1.5*(2.-2.*F+bra*F)/F*(rho_o+p_o)/(H*H)
      + bra_N;

    double D = Dfloor;
    if (cs2num > 0.) {
      double Dt = cs2num/cs2_target;
      if (Dt > D) D = Dt;
    }
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
)

print("SDMC localized transition-bridge patch complete.")
