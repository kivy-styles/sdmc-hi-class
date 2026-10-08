#!/usr/bin/env python3
from pathlib import Path
p=Path("gravity_smg/gravity_models_smg.c")
s=p.read_text()
old='''    double S=.5*(1.+th);
    double F=exp(AF*S);
    double am=AF/(4.*width)*sech2;
    double bra=-2.*am;
    double D=Dfloor + D0*pow(S,power);'''
new='''    double S=.5*(1.+th);
    double F=exp(AF*S);
    double am=AF/(4.*width)*sech2;
    double bra=-2.*am;
    /* SDMC_F_ONLY_FORCE_CENTER:
       keep the accepted independent-kinetic response on its original
       structural window while allowing the Planck-mass response F to use
       an independently derived force-balance center/width. */
    const double sdmc_D_zc=3.927876388467848;
    const double sdmc_D_width=0.33114133956842123;
    double sdmc_D_Nc=-log(1.+sdmc_D_zc);
    double sdmc_D_x=(N-sdmc_D_Nc)/(2.*sdmc_D_width);
    double sdmc_D_S=.5*(1.+tanh(sdmc_D_x));
    double D=Dfloor + D0*pow(sdmc_D_S,power);'''
if old not in s:
    raise SystemExit("independent kinetic anchor not found")
s=s.replace(old,new,1)
p.write_text(s)
print("F_ONLY_FORCE_CENTER_D_WINDOW_DECOUPLED")
