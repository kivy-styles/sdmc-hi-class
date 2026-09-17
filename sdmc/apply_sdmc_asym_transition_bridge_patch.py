#!/usr/bin/env python3
"""Refine the Kp transition bridge to use independent turn-on/off widths.

Apply this *after* apply_sdmc_transition_bridge_patch.py.

The signed symmetric-pulse scan showed a characteristic edge migration:
positive A_bridge repairs the early edge but transfers the native c_s^2
failure to the turn-off edge, while negative A_bridge fails at the turn-on
edge. A compact bridge therefore needs asymmetric derivatives: a relatively
sharp positive turn-on correction near the high-z instability and a broad,
late negative return where the unmodified late sector has more stability
margin.

This patch changes only the Kp bridge parameterization from

  cs2, Dfloor, A_bridge, z_on, z_off, width

to

  cs2, Dfloor, A_bridge, z_on, z_off, width_on, width_off.

It leaves F(a), alpha_M(a), the frozen Kp expansion history, and the NKp-v2
bridge implementation untouched.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: asymmetric Kp bridge already present")
        return
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{path}: expected one asymmetric-bridge anchor, found {n}\n{old}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: asymmetric Kp bridge patched")


replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''     pba->parameters_2_size_smg = 6;\n     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\n     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||\n                pba->parameters_2_smg[3] <= 0. || pba->parameters_2_smg[4] <= 0. ||\n                pba->parameters_2_smg[5] <= 0. || pba->parameters_2_smg[3] <= pba->parameters_2_smg[4], errmsg,\n                "sdmc_kp_transition_bridge requires positive cs2, D floor, z_on, z_off, bridge width and z_on>z_off");''',
    '''     pba->parameters_2_size_smg = 7;\n     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\n     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||\n                pba->parameters_2_smg[3] <= 0. || pba->parameters_2_smg[4] <= 0. ||\n                pba->parameters_2_smg[5] <= 0. || pba->parameters_2_smg[6] <= 0. ||\n                pba->parameters_2_smg[3] <= pba->parameters_2_smg[4], errmsg,\n                "sdmc_kp_transition_bridge requires positive cs2, D floor, z_on, z_off, turn-on/off widths and z_on>z_off");'''
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    /* parameters_smg = cs2, Dfloor, A_bridge, z_on, z_off, bridge_width */\n    const double AM = 0.02586;\n    const double zc = 4.38;\n    const double width = 0.4467;\n    double cs2_target = pba->parameters_2_smg[0];\n    double Dfloor = pba->parameters_2_smg[1];\n    double Abridge = pba->parameters_2_smg[2];\n    double zon = pba->parameters_2_smg[3];\n    double zoff = pba->parameters_2_smg[4];\n    double wbri = pba->parameters_2_smg[5];\n    double N = log(a);''',
    '''    /* parameters_smg = cs2, Dfloor, A_bridge, z_on, z_off, width_on, width_off */\n    const double AM = 0.02586;\n    const double zc = 4.38;\n    const double width = 0.4467;\n    double cs2_target = pba->parameters_2_smg[0];\n    double Dfloor = pba->parameters_2_smg[1];\n    double Abridge = pba->parameters_2_smg[2];\n    double zon = pba->parameters_2_smg[3];\n    double zoff = pba->parameters_2_smg[4];\n    double won = pba->parameters_2_smg[5];\n    double woff = pba->parameters_2_smg[6];\n    double N = log(a);'''
)

# Anchor the pulse replacement to the Kp-only width_on/width_off declarations.
# This avoids matching the separate NKp-v2 symmetric bridge, which intentionally
# retains its single wbri parameter.
replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''    double won = pba->parameters_2_smg[5];\n    double woff = pba->parameters_2_smg[6];\n    double N = log(a);\n\n    double Nc = -log(1.+zc);\n    double x = (N-Nc)/(2.*width);\n    double th = tanh(x);\n    double sech2 = 1./(cosh(x)*cosh(x));\n    double S = .5*(1.+th);\n    double F = exp(AM*S);\n    double am = AM/(4.*width)*sech2;\n    double amp = -AM/(4.*width*width)*sech2*th;\n\n    double Non = -log(1.+zon);\n    double Noff = -log(1.+zoff);\n    double ton = tanh((N-Non)/wbri);\n    double toff = tanh((N-Noff)/wbri);\n    double pulse = .5*(ton-toff);\n    double pulse_N = .5*((1.-ton*ton)-(1.-toff*toff))/wbri;\n\n    double bra = -2.*am + Abridge*pulse;''',
    '''    double won = pba->parameters_2_smg[5];\n    double woff = pba->parameters_2_smg[6];\n    double N = log(a);\n\n    double Nc = -log(1.+zc);\n    double x = (N-Nc)/(2.*width);\n    double th = tanh(x);\n    double sech2 = 1./(cosh(x)*cosh(x));\n    double S = .5*(1.+th);\n    double F = exp(AM*S);\n    double am = AM/(4.*width)*sech2;\n    double amp = -AM/(4.*width*width)*sech2*th;\n\n    double Non = -log(1.+zon);\n    double Noff = -log(1.+zoff);\n    double ton = tanh((N-Non)/won);\n    double toff = tanh((N-Noff)/woff);\n    double pulse = .5*(ton-toff);\n    double pulse_N = .5*((1.-ton*ton)/won-(1.-toff*toff)/woff);\n\n    double bra = -2.*am + Abridge*pulse;'''
)

print("SDMC asymmetric Kp transition-bridge runtime patch complete.")
