#!/usr/bin/env python3
from pathlib import Path

def replace_once(path, old, new):
    p=Path(path)
    s=p.read_text()
    if new in s:
        return
    if s.count(old)!=1:
        raise RuntimeError(f"{path}: expected one anchor, found {s.count(old)}")
    p.write_text(s.replace(old,new,1))

# This installer is meant to run after install_sdmc_v3_independent_kinetic.py.
replace_once(
    "include/background.h",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip, sdmc_v3_independent_kinetic,\n",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip, sdmc_v3_independent_kinetic, sdmc_v3_lens_braiding,\n"
)

anchor='''  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
parser='''  if (strcmp(string1,"sdmc_v3_lens_braiding") == 0) {
     pba->gravity_model_smg = sdmc_v3_lens_braiding;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 9;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||
                pba->parameters_2_smg[2] <= 0. || pba->parameters_2_smg[3] <= 0. ||
                pba->parameters_2_smg[4] < 0. || pba->parameters_2_smg[5] <= 0. ||
                pba->parameters_2_smg[7] < 0. || pba->parameters_2_smg[8] <= 0., errmsg,
                "sdmc_v3_lens_braiding requires AF,zc,width,D0,Dfloor,z_lens,sigma_lens positive and p non-negative");
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
replace_once("gravity_smg/gravity_models_smg.c",anchor,parser)

canchor='''  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
closure=r'''  else if (pba->gravity_model_smg == sdmc_v3_lens_braiding) {
    double AF=pba->parameters_2_smg[0];
    double zc=pba->parameters_2_smg[1];
    double width=pba->parameters_2_smg[2];
    double D0=pba->parameters_2_smg[3];
    double power=pba->parameters_2_smg[4];
    double Dfloor=pba->parameters_2_smg[5];
    double Alens=pba->parameters_2_smg[6];
    double zlens=pba->parameters_2_smg[7];
    double siglens=pba->parameters_2_smg[8];

    double N=log(a);
    double Nc=-log(1.+zc);
    double x=(N-Nc)/(2.*width);
    double th=tanh(x);
    double ch=cosh(x);
    double sech2=1./(ch*ch);
    double S=.5*(1.+th);
    double F=exp(AF*S);
    double am=AF/(4.*width)*sech2;

    double z=1./a-1.;
    double dz=(z-zlens)/siglens;
    double dlens=Alens*exp(-0.5*dz*dz);
    double bra=-2.*am + dlens;

    /* Hold the independent kinetic normalization D fixed while changing the
       braiding relation, so this diagnostic isolates the perturbation/slip
       effect rather than conflating it with a new kinetic normalization. */
    double D=Dfloor + D0*pow(S,power);
    double ak=D-1.5*bra*bra;

    pvecback[pba->index_bg_kineticity_smg]=ak;
    pvecback[pba->index_bg_braiding_smg]=bra;
    pvecback[pba->index_bg_tensor_excess_smg]=0.;
    pvecback[pba->index_bg_M2_running_smg]=am;
    pvecback[pba->index_bg_beyond_horndeski_smg]=0.;
    pvecback[pba->index_bg_delta_M2_smg]=F-1.;
    pvecback[pba->index_bg_M2_smg]=F;
  }

  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
replace_once("gravity_smg/gravity_models_smg.c",canchor,closure)
print("LENS_BRAIDING_MODEL_INSTALLED")
