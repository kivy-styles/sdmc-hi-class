#!/usr/bin/env python3
from pathlib import Path

# Freeze the split-tracker helper to the combined candidate's established
# background constants when the true covariant model uses parameters_smg for
# its Lagrangian rather than for expansion_smg.
p=Path("source/background.c")
s=p.read_text()
old='''  double Ox = pba->parameters_smg[0];
  double lambda_e = pba->parameters_smg[1];
  double zt = pba->parameters_smg[2];
  double dNt = pba->parameters_smg[3];
  double A = pba->parameters_smg[4];
  double tauA = pba->parameters_smg[5];
  double B = pba->parameters_smg[6];
  double tauB = pba->parameters_smg[7];'''
new='''  double Ox = 0.6953970160452582;
  double lambda_e = 17.925;
  double zt = 17.775;
  double dNt = 0.5;
  double A = 0.01105624999;
  double tauA = 0.25;
  double B = 0.01951933685;
  double tauB = 1.5;'''
if old not in s:
    raise RuntimeError("split tracker constant anchor missing")
p.write_text(s.replace(old,new,1))

# Add enum entry after the already-installed EFT diagnostic models.
p=Path("include/background.h")
s=p.read_text()
anchor="sdmc_v3_lens_braiding,"
if anchor not in s:
    raise RuntimeError("lens-braiding enum anchor missing")
s=s.replace(anchor,anchor+" sdmc_combined_covariant_replay,",1)
p.write_text(s)

p=Path("gravity_smg/gravity_models_smg.c")
s=p.read_text()
if '#include "sdmc_combined_covariant_table.h"' not in s:
    s=s.replace('#include "gravity_models_smg.h"\n',
                '#include "gravity_models_smg.h"\n#include "sdmc_combined_covariant_table.h"\n',1)

anchor='''  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
parser='''  if (strcmp(string1,"sdmc_combined_covariant_replay") == 0) {
     pba->gravity_model_smg = sdmc_combined_covariant_replay;
     pba->field_evolution_smg = _TRUE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_size_smg = 1;
     class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);
     if (has_tuning_index_smg == _FALSE_) {
       pba->tuning_index_smg = 0;
       pba->tuning_dxdy_guess_smg = 1.;
     }
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
if s.count(anchor)!=1:
    raise RuntimeError("parser anchor count "+str(s.count(anchor)))
s=s.replace(anchor,parser,1)

ganchor='''  else if(pba->gravity_model_smg == galileon){'''
gcode=r'''  else if (pba->gravity_model_smg == sdmc_combined_covariant_replay) {
    double Nc=-log(1.+SDMC_COMB_COV_ZC);
    double xx=(phi-Nc)/(2.*SDMC_COMB_COV_WIDTH);
    double tt=tanh(xx);
    double uu=1.-tt*tt;
    double w=SDMC_COMB_COV_WIDTH;
    double S=.5*(1.+tt);
    double S1=uu/(4.*w);
    double S2=-tt*uu/(4.*w*w);
    double S3=-uu*(1.-3.*tt*tt)/(8.*w*w*w);
    double A=SDMC_COMB_COV_AF;
    double F=exp(A*S);
    double F1=F*(A*S1);
    double F2=F*(A*A*S1*S1+A*S2);
    double F3=F*(A*A*A*S1*S1*S1+3.*A*A*S1*S2+A*S3);

    double g,g1,g2,g3;
    double k1,k1p,k1pp,k1ppp;
    double k2,k2p,k2pp,k2ppp;
    double V,Vp,Vpp,Vppp;
    sdmc_comb_cov_eval(sdmc_comb_cov_g,phi,&g,&g1,&g2,&g3);
    sdmc_comb_cov_eval(sdmc_comb_cov_k1,phi,&k1,&k1p,&k1pp,&k1ppp);
    sdmc_comb_cov_eval(sdmc_comb_cov_k2,phi,&k2,&k2p,&k2pp,&k2ppp);
    sdmc_comb_cov_eval(sdmc_comb_cov_V,phi,&V,&Vp,&Vpp,&Vppp);
    V += 3.*pba->H0*pba->H0*pba->parameters_smg[0];

    pgf->G2 = k1*X+k2*X*X-V;
    pgf->G2_X = k1+2.*k2*X;
    pgf->G2_XX = 2.*k2;
    pgf->G2_phi = k1p*X+k2p*X*X-Vp;
    pgf->G2_Xphi = k1p+2.*k2p*X;
    pgf->G2_XXphi = 2.*k2p;
    pgf->G2_phiphi = k1pp*X+k2pp*X*X-Vpp;
    pgf->G2_Xphiphi = k1pp+2.*k2pp*X;

    pgf->G3_X = g;
    pgf->G3_phi = g1*X;
    pgf->G3_Xphi = g1;
    pgf->G3_phiphi = g2*X;
    pgf->G3_Xphiphi = g2;
    pgf->G3_phiphiphi = g3*X;

    pgf->DG4 = (F-1.)/2.;
    pgf->G4 = F/2.;
    pgf->G4_phi = F1/2.;
    pgf->G4_phiphi = F2/2.;
    pgf->G4_phiphiphi = F3/2.;
  }

  else if(pba->gravity_model_smg == galileon){'''
if s.count(ganchor)!=1:
    raise RuntimeError("G anchor count "+str(s.count(ganchor)))
s=s.replace(ganchor,gcode,1)

start=s.index("int gravity_models_initial_conditions_smg(")
pos=s.index("case propto_omega:",start)
ic='''case sdmc_combined_covariant_replay:
    {
      double ph=log(a);
      double lnh,lnhp,lnhpp,lnhppp;
      sdmc_comb_cov_eval(sdmc_comb_cov_lnH,ph,&lnh,&lnhp,&lnhpp,&lnhppp);
      pvecback_integration[pba->index_bi_phi_smg]=ph;
      pvecback_integration[pba->index_bi_phi_prime_smg]=a*exp(lnh);
    }
    break;

  '''
s=s[:pos]+ic+s[pos:]
p.write_text(s)
print("COMBINED_COVARIANT_MODEL_INSTALLED")
