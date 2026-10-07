#!/usr/bin/env python3
from pathlib import Path

p=Path("gravity_smg/gravity_models_smg.c")
s=p.read_text()

marker='pba->gravity_model_smg = sdmc_v3_covariant_linear_audit;'
i=s.index(marker)
j=s.index('if (strcmp(string1,"eft_alphas_power_law")',i)
block=s[i:j]
if "pba->parameters_size_smg = 3;" not in block:
    raise RuntimeError("expected basin-audit parser size=3 before null stabilizer patch")
block=block.replace("pba->parameters_size_smg = 3;","pba->parameters_size_smg = 4;",1)
s=s[:i]+block+s[j:]

anchor='    pgf->G2_Xphiphi = k1pp+2.*k2pp*X;\n'
if anchor not in s:
    raise RuntimeError("G2 derivative anchor missing")

insert=r'''    pgf->G2_Xphiphi = k1pp+2.*k2pp*X;

    /* Exploratory off-trajectory null deformation:
       dG2 = eps * (X-Xs(phi))^4 / Xs(phi)^3,
       Xs(phi)=H_target(phi)^2/2.
       dG2 and all derivatives used by hi_class through total order 3
       vanish identically on the reconstructed late300 trajectory X=Xs.
       This therefore tests covariant-completion non-uniqueness without
       changing the target background/linear action on trajectory. */
    {
      double eps_null = pba->parameters_smg[3];
      double lnhs,lnh1,lnh2,lnh3;
      sdmc_lin_eval(sdmc_lin_lnH,phi,&lnhs,&lnh1,&lnh2,&lnh3);
      double Xs = 0.5*exp(2.*lnhs);
      double Xsp = 2.*lnh1*Xs;
      double Xspp = (2.*lnh2 + 4.*lnh1*lnh1)*Xs;
      double Y = X-Xs;
      double Q2 = Xs*Xs;
      double Q3 = Q2*Xs;
      double Q4 = Q3*Xs;
      double Q5 = Q4*Xs;
      double Y2 = Y*Y;
      double Y3 = Y2*Y;
      double Y4 = Y2*Y2;

      double d0 = eps_null*Y4/Q3;
      double dX = 4.*eps_null*Y3/Q3;
      double dXX = 12.*eps_null*Y2/Q3;
      double dXXX = 24.*eps_null*Y/Q3;
      double dphi = -eps_null*Xsp*Y3*(3.*X+Xs)/Q4;
      double dXphi = -12.*eps_null*X*Xsp*Y2/Q4;
      double dXXphi = 12.*eps_null*Y*(Xs-3.*X)*Xsp/Q4;
      double dphiphi = eps_null*Y2*
        (-3.*X*X*Xs*Xspp + 12.*X*X*Xsp*Xsp
         +2.*X*Xs*Xs*Xspp + Xs*Xs*Xs*Xspp)/Q5;
      double dXphiphi = 12.*X*eps_null*Y*
        (-Xs*Xspp*Y + 2.*Xsp*Xsp*(2.*X-Xs))/Q5;

      pgf->G2 += d0;
      pgf->G2_X += dX;
      pgf->G2_XX += dXX;
      pgf->G2_XXX += dXXX;
      pgf->G2_phi += dphi;
      pgf->G2_Xphi += dXphi;
      pgf->G2_XXphi += dXXphi;
      pgf->G2_phiphi += dphiphi;
      pgf->G2_Xphiphi += dXphiphi;
    }
'''
s=s.replace(anchor,insert,1)
p.write_text(s)
print("NULL_G2_STABILIZER_PATCHED parameters_size=4 epsilon_index=3")
