#!/usr/bin/env python3
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument("--OmegaX",type=float,required=True)
ap.add_argument("--lambda-e",type=float,default=17.925)
ap.add_argument("--z-t",type=float,default=17.775)
ap.add_argument("--dNt",type=float,default=0.5)
ap.add_argument("--A",type=float,default=0.01105624999)
ap.add_argument("--tauA",type=float,default=0.25)
ap.add_argument("--B",type=float,default=0.01951933685)
ap.add_argument("--tauB",type=float,default=1.5)
args=ap.parse_args()

# The split tracker helper normally reads expansion_smg. The covariant replay
# instead consumes parameters_smg for the Lagrangian, so freeze the exact same
# background constants that were used to generate this candidate's target.
p=Path("source/background.c"); s=p.read_text()
old='''  double Ox = pba->parameters_smg[0];
  double lambda_e = pba->parameters_smg[1];
  double zt = pba->parameters_smg[2];
  double dNt = pba->parameters_smg[3];
  double A = pba->parameters_smg[4];
  double tauA = pba->parameters_smg[5];
  double B = pba->parameters_smg[6];
  double tauB = pba->parameters_smg[7];'''
# tolerate the historical indentation used by the patcher
if old not in s:
    old='''  double Ox = pba->parameters_smg[0];
    double lambda_e = pba->parameters_smg[1];
    double zt = pba->parameters_smg[2];
    double dNt = pba->parameters_smg[3];
    double A = pba->parameters_smg[4];
    double tauA = pba->parameters_smg[5];
    double B = pba->parameters_smg[6];
    double tauB = pba->parameters_smg[7];'''
new=f'''  double Ox = {args.OmegaX:.17g};
    double lambda_e = {args.lambda_e:.17g};
    double zt = {args.z_t:.17g};
    double dNt = {args.dNt:.17g};
    double A = {args.A:.17g};
    double tauA = {args.tauA:.17g};
    double B = {args.B:.17g};
    double tauB = {args.tauB:.17g};'''
if old not in s:
    raise RuntimeError("split tracker constant anchor missing")
p.write_text(s.replace(old,new,1))

p=Path("include/background.h"); s=p.read_text()
anchors=[
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip, sdmc_v3_independent_kinetic,\n",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip,\n",
]
done=False
for old_enum in anchors:
    if old_enum in s:
        new_enum=old_enum.rstrip("\n").rstrip(",")+" , sdmc_covariant_replay,\n"
        # repair spacing deterministically
        new_enum=new_enum.replace("noslip ,","noslip,").replace("kinetic ,","kinetic,")
        s=s.replace(old_enum,new_enum,1)
        done=True
        break
if not done:
    raise RuntimeError("enum anchor missing")
p.write_text(s)

p=Path("gravity_smg/gravity_models_smg.c"); s=p.read_text()
if '#include "sdmc_covariant_replay_table.h"' not in s:
    inc='#include "gravity_models_smg.h"\n'
    if inc not in s: raise RuntimeError("gravity include anchor missing")
    s=s.replace(inc,inc+'#include "sdmc_covariant_replay_table.h"\n',1)

anchor='''  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
parser='''  if (strcmp(string1,"sdmc_covariant_replay") == 0) {
     pba->gravity_model_smg = sdmc_covariant_replay;
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
if s.count(anchor)!=1: raise RuntimeError(f"parser anchor count {s.count(anchor)}")
s=s.replace(anchor,parser,1)

ganchor='''  else if(pba->gravity_model_smg == galileon){'''
gcode=r'''  else if (pba->gravity_model_smg == sdmc_covariant_replay) {
    double Nc=-log(1.+SDMC_COV_ZC);
    double xx=(phi-Nc)/(2.*SDMC_COV_WIDTH);
    double tt=tanh(xx);
    double uu=1.-tt*tt;
    double w=SDMC_COV_WIDTH;
    double S=.5*(1.+tt);
    double S1=uu/(4.*w);
    double S2=-tt*uu/(4.*w*w);
    double S3=-uu*(1.-3.*tt*tt)/(8.*w*w*w);
    double A=SDMC_COV_AF;
    double F=exp(A*S);
    double F1=F*(A*S1);
    double F2=F*(A*A*S1*S1+A*S2);
    double F3=F*(A*A*A*S1*S1*S1+3.*A*A*S1*S2+A*S3);

    double g,g1,g2,g3;
    double k1,k1p,k1pp,k1ppp;
    double k2,k2p,k2pp,k2ppp;
    double V,Vp,Vpp,Vppp;
    sdmc_cov_eval(sdmc_cov_g,phi,&g,&g1,&g2,&g3);
    sdmc_cov_eval(sdmc_cov_k1,phi,&k1,&k1p,&k1pp,&k1ppp);
    sdmc_cov_eval(sdmc_cov_k2,phi,&k2,&k2p,&k2pp,&k2ppp);
    sdmc_cov_eval(sdmc_cov_V,phi,&V,&Vp,&Vpp,&Vppp);
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
if s.count(ganchor)!=1: raise RuntimeError(f"G anchor count {s.count(ganchor)}")
s=s.replace(ganchor,gcode,1)

start=s.index("int gravity_models_initial_conditions_smg(")
pos=s.index("case propto_omega:",start)
ic='''case sdmc_covariant_replay:
    {
      double ph=log(a);
      double lnh,lnhp,lnhpp,lnhppp;
      sdmc_cov_eval(sdmc_cov_lnH,ph,&lnh,&lnhp,&lnhpp,&lnhppp);
      pvecback_integration[pba->index_bi_phi_smg]=ph;
      pvecback_integration[pba->index_bi_phi_prime_smg]=a*exp(lnh);
    }
    break;

  '''
s=s[:pos]+ic+s[pos:]
p.write_text(s)
print("COVARIANT_REPLAY_INSTALLED",
      {"OmegaX":args.OmegaX,"lambda_e":args.lambda_e,"z_t":args.z_t},flush=True)
