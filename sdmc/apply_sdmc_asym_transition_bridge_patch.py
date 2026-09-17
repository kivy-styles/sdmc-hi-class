#!/usr/bin/env python3
"""Refine the Kp transition bridge with three diagnostic modes.

Apply this after apply_sdmc_transition_bridge_patch.py.

Modes selected by z_off:

  z_off > 0 : asymmetric alpha_B pulse with independent turn-on/off widths.
  z_off = 0 : one-way alpha_B step with no finite-redshift return edge.
  z_off < 0 : localized Planck-mass bridge.  A_bridge is interpreted as the
              amplitude of a sech^2 bump in ln F centered at z_on with width
              width_on.  In this mode alpha_B=-2 alpha_M exactly, so No-Slip is
              preserved while F, alpha_M and alpha_M' are deformed only inside
              a compact transition region and return to the frozen Kp branch.

The third mode is motivated by two native-runtime null results: neither a
finite alpha_B-only pulse nor a one-way residual alpha_B offset produced a
stable Kp completion.  It tests whether the obstruction belongs to the frozen
Planck-mass trajectory itself rather than to the exact No-Slip relation.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    n = text.count(old)
    if n == 0 and new in text:
        print(f"{path}: extended Kp bridge already present")
        return
    if n != 1:
        raise RuntimeError(f"{path}: expected one bridge anchor, found {n}\n{old}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: extended Kp bridge patched")


replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''     pba->parameters_2_size_smg = 6;\
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||\
                pba->parameters_2_smg[3] <= 0. || pba->parameters_2_smg[4] <= 0. ||\
                pba->parameters_2_smg[5] <= 0. || pba->parameters_2_smg[3] <= pba->parameters_2_smg[4], errmsg,\
                "sdmc_kp_transition_bridge requires positive cs2, D floor, z_on, z_off, bridge width and z_on>z_off");''',
    '''     pba->parameters_2_size_smg = 7;\
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||\
                pba->parameters_2_smg[3] <= 0. ||\
                pba->parameters_2_smg[5] <= 0. || pba->parameters_2_smg[6] <= 0. ||\
                (pba->parameters_2_smg[4] > 0. && pba->parameters_2_smg[3] <= pba->parameters_2_smg[4]), errmsg,\
                "sdmc_kp_transition_bridge requires positive cs2, D floor, z_on, widths; z_on>z_off when z_off>0");'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    /* parameters_smg = cs2, Dfloor, A_bridge, z_on, z_off, bridge_width */\
    const double AM = 0.02586;\
    const double zc = 4.38;\
    const double width = 0.4467;\
    double cs2_target = pba->parameters_2_smg[0];\
    double Dfloor = pba->parameters_2_smg[1];\
    double Abridge = pba->parameters_2_smg[2];\
    double zon = pba->parameters_2_smg[3];\
    double zoff = pba->parameters_2_smg[4];\
    double wbri = pba->parameters_2_smg[5];\
    double N = log(a);''',
    '''    /* parameters_smg = cs2, Dfloor, A_bridge, z_on, z_off, width_on, width_off\
     * z_off = 0  : one-way alpha_B step.\
     * z_off < 0  : localized ln F sech^2 bridge, exact alpha_B=-2 alpha_M. */\
    const double AM = 0.02586;\
    const double zc = 4.38;\
    const double width = 0.4467;\
    double cs2_target = pba->parameters_2_smg[0];\
    double Dfloor = pba->parameters_2_smg[1];\
    double Abridge = pba->parameters_2_smg[2];\
    double zon = pba->parameters_2_smg[3];\
    double zoff = pba->parameters_2_smg[4];\
    double won = pba->parameters_2_smg[5];\
    double woff = pba->parameters_2_smg[6];\
    double N = log(a);'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    double won = pba->parameters_2_smg[5];\
    double woff = pba->parameters_2_smg[6];\
    double N = log(a);\
\
    double Nc = -log(1.+zc);\
    double x = (N-Nc)/(2.*width);\
    double th = tanh(x);\
    double sech2 = 1./(cosh(x)*cosh(x));\
    double S = .5*(1.+th);\
    double F = exp(AM*S);\
    double am = AM/(4.*width)*sech2;\
    double amp = -AM/(4.*width*width)*sech2*th;\
\
    double Non = -log(1.+zon);\
    double Noff = -log(1.+zoff);\
    double ton = tanh((N-Non)/wbri);\
    double toff = tanh((N-Noff)/wbri);\
    double pulse = .5*(ton-toff);\
    double pulse_N = .5*((1.-ton*ton)-(1.-toff*toff))/wbri;\
\
    double bra = -2.*am + Abridge*pulse;''',
    '''    double won = pba->parameters_2_smg[5];\
    double woff = pba->parameters_2_smg[6];\
    double N = log(a);\
\
    double Nc = -log(1.+zc);\
    double x = (N-Nc)/(2.*width);\
    double th = tanh(x);\
    double sech2 = 1./(cosh(x)*cosh(x));\
    double S = .5*(1.+th);\
    double lnF = AM*S;\
    double am = AM/(4.*width)*sech2;\
    double amp = -AM/(4.*width*width)*sech2*th;\
\
    double Non = -log(1.+zon);\
    double bra, bra_N;\
    if (zoff < 0.) {\
      /* Localized Planck-mass bridge in ln F.  q=sech^2(u) vanishes on\
       * both sides, so the original Kp normalization is recovered exactly. */\
      double u = (N-Non)/won;\
      double tb = tanh(u);\
      double qb = 1./(cosh(u)*cosh(u));\
      lnF += Abridge*qb;\
      am += -2.*Abridge*qb*tb/won;\
      amp += 2.*Abridge*qb*(3.*tb*tb-1.)/(won*won);\
      bra = -2.*am;\
      bra_N = -2.*amp;\
    }\
    else {\
      double ton = tanh((N-Non)/won);\
      double pulse, pulse_N;\
      if (zoff == 0.) {\
        pulse = .5*(1.+ton);\
        pulse_N = .5*(1.-ton*ton)/won;\
      }\
      else {\
        double Noff = -log(1.+zoff);\
        double toff = tanh((N-Noff)/woff);\
        pulse = .5*(ton-toff);\
        pulse_N = .5*((1.-ton*ton)/won-(1.-toff*toff)/woff);\
      }\
      bra = -2.*am + Abridge*pulse;\
      bra_N = -2.*amp + Abridge*pulse_N;\
    }\
    double F = exp(lnF);'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    double bra_N = -2.*amp + Abridge*pulse_N;\
    double H = pvecback[pba->index_bg_H];''',
    '''    double H = pvecback[pba->index_bg_H];'''
)

print("SDMC extended Kp transition-bridge runtime patch complete.")
