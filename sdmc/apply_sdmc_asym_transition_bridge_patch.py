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
    if new in text:
        print(f"{path}: extended Kp bridge already present")
        return
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{path}: expected one bridge anchor, found {n}\n{old}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: extended Kp bridge patched")


replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''     pba->parameters_2_size_smg = 6;\n     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\n     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||\n                pba->parameters_2_smg[3] <= 0. || pba->parameters_2_smg[4] <= 0. ||\n                pba->parameters_2_smg[5] <= 0. || pba->parameters_2_smg[3] <= pba->parameters_2_smg[4], errmsg,\n                "sdmc_kp_transition_bridge requires positive cs2, D floor, z_on, z_off, bridge width and z_on>z_off");''',
    '''     pba->parameters_2_size_smg = 7;\n     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\n     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||\n                pba->parameters_2_smg[3] <= 0. ||\n                pba->parameters_2_smg[5] <= 0. || pba->parameters_2_smg[6] <= 0. ||\n                (pba->parameters_2_smg[4] > 0. && pba->parameters_2_smg[3] <= pba->parameters_2_smg[4]), errmsg,\n                "sdmc_kp_transition_bridge requires positive cs2, D floor, z_on, widths; z_on>z_off when z_off>0");'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    /* parameters_smg = cs2, Dfloor, A_bridge, z_on, z_off, bridge_width */\n    const double AM = 0.02586;\n    const double zc = 4.38;\n    const double width = 0.4467;\n    double cs2_target = pba->parameters_2_smg[0];\n    double Dfloor = pba->parameters_2_smg[1];\n    double Abridge = pba->parameters_2_smg[2];\n    double zon = pba->parameters_2_smg[3];\n    double zoff = pba->parameters_2_smg[4];\n    double wbri = pba->parameters_2_smg[5];\n    double N = log(a);''',
    '''    /* parameters_smg = cs2, Dfloor, A_bridge, z_on, z_off, width_on, width_off\n     * z_off = 0  : one-way alpha_B step.\n     * z_off < 0  : localized ln F sech^2 bridge, exact alpha_B=-2 alpha_M. */\n    const double AM = 0.02586;\n    const double zc = 4.38;\n    const double width = 0.4467;\n    double cs2_target = pba->parameters_2_smg[0];\n    double Dfloor = pba->parameters_2_smg[1];\n    double Abridge = pba->parameters_2_smg[2];\n    double zon = pba->parameters_2_smg[3];\n    double zoff = pba->parameters_2_smg[4];\n    double won = pba->parameters_2_smg[5];\n    double woff = pba->parameters_2_smg[6];\n    double N = log(a);'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    double won = pba->parameters_2_smg[5];\n    double woff = pba->parameters_2_smg[6];\n    double N = log(a);\n\n    double Nc = -log(1.+zc);\n    double x = (N-Nc)/(2.*width);\n    double th = tanh(x);\n    double sech2 = 1./(cosh(x)*cosh(x));\n    double S = .5*(1.+th);\n    double F = exp(AM*S);\n    double am = AM/(4.*width)*sech2;\n    double amp = -AM/(4.*width*width)*sech2*th;\n\n    double Non = -log(1.+zon);\n    double Noff = -log(1.+zoff);\n    double ton = tanh((N-Non)/wbri);\n    double toff = tanh((N-Noff)/wbri);\n    double pulse = .5*(ton-toff);\n    double pulse_N = .5*((1.-ton*ton)-(1.-toff*toff))/wbri;\n\n    double bra = -2.*am + Abridge*pulse;''',
    '''    double won = pba->parameters_2_smg[5];\n    double woff = pba->parameters_2_smg[6];\n    double N = log(a);\n\n    double Nc = -log(1.+zc);\n    double x = (N-Nc)/(2.*width);\n    double th = tanh(x);\n    double sech2 = 1./(cosh(x)*cosh(x));\n    double S = .5*(1.+th);\n    double lnF = AM*S;\n    double am = AM/(4.*width)*sech2;\n    double amp = -AM/(4.*width*width)*sech2*th;\n\n    double Non = -log(1.+zon);\n    double bra, bra_N;\n    if (zoff < 0.) {\n      /* Localized Planck-mass bridge in ln F.  q=sech^2(u) vanishes on\n       * both sides, so the original Kp normalization is recovered exactly. */\n      double u = (N-Non)/won;\n      double tb = tanh(u);\n      double qb = 1./(cosh(u)*cosh(u));\n      lnF += Abridge*qb;\n      am += -2.*Abridge*qb*tb/won;\n      amp += 2.*Abridge*qb*(3.*tb*tb-1.)/(won*won);\n      bra = -2.*am;\n      bra_N = -2.*amp;\n    }\n    else {\n      double ton = tanh((N-Non)/won);\n      double pulse, pulse_N;\n      if (zoff == 0.) {\n        pulse = .5*(1.+ton);\n        pulse_N = .5*(1.-ton*ton)/won;\n      }\n      else {\n        double Noff = -log(1.+zoff);\n        double toff = tanh((N-Noff)/woff);\n        pulse = .5*(ton-toff);\n        pulse_N = .5*((1.-ton*ton)/won-(1.-toff*toff)/woff);\n      }\n      bra = -2.*am + Abridge*pulse;\n      bra_N = -2.*amp + Abridge*pulse_N;\n    }\n    double F = exp(lnF);'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    double bra_N = -2.*amp + Abridge*pulse_N;\n    double H = pvecback[pba->index_bg_H];''',
    '''    double H = pvecback[pba->index_bg_H];'''
)

print("SDMC extended Kp transition-bridge runtime patch complete.")
