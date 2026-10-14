#!/usr/bin/env python3
"""Add a generalized NKp-v2 EFT-completion diagnostic to hi_class.

This model deliberately keeps the frozen NKp-v2 expansion history unchanged.
It uses the smooth B0 profile only as a *shape basis* and scans

    F = F_B0**s,
    alpha_M = s alpha_M,B0,
    alpha_B = beta_B alpha_M,
    D = alpha_K + 3 alpha_B^2/2 = D0 > 0.

Hence alpha_K = D0 - 3 alpha_B^2/2.  Native hi_class computes the actual
scalar sound speed and rejects gradient-unstable choices.  This is a
completion-space diagnostic, not a refit and not a new frozen model.

parameters_smg = s, beta_B, D0, Ns_B0, z_activation, DeltaN_activation
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: generalized completion patch already present")
        return
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{path}: expected one patch anchor, found {n}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: patched")


# The script is applied after apply_sdmc_v2_noslip_patch.py.
replace_once(
    "include/background.h",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip,\n",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_completion,\n",
)

p = Path("gravity_smg/gravity_models_smg.c")
text = p.read_text(encoding="utf-8")
marker = "int gravity_models_get_alphas_par_smg("
if text.count(marker) != 1:
    raise RuntimeError("could not uniquely locate gravity_models_get_alphas_par_smg")
pre, post = text.split(marker, 1)

parser_anchor = '  if (strcmp(string1,"eft_alphas_power_law") == 0) {'
parser_block = '''  if (strcmp(string1,"sdmc_v2_completion") == 0) {
     pba->gravity_model_smg = sdmc_v2_completion;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 6;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] < 0. || pba->parameters_2_smg[2] <= 0., errmsg,
                "sdmc_v2_completion requires s>=0 and D0>0");
     class_test(pba->parameters_2_smg[3] <= 0. || pba->parameters_2_smg[4] <= 0. || pba->parameters_2_smg[5] <= 0., errmsg,
                "sdmc_v2_completion requires positive B0 seed/activation parameters");
   }

'''
if parser_block not in pre:
    if pre.count(parser_anchor) != 1:
        raise RuntimeError(f"expected one parser anchor before alpha getter, found {pre.count(parser_anchor)}")
    pre = pre.replace(parser_anchor, parser_block + parser_anchor)

model_anchor = "  else if (pba->gravity_model_smg == eft_alphas_power_law) {"
model_block = '''  else if (pba->gravity_model_smg == sdmc_v2_completion) {

    double scale = pba->parameters_2_smg[0];
    double betaB = pba->parameters_2_smg[1];
    double D0 = pba->parameters_2_smg[2];
    double Ns_seed = pba->parameters_2_smg[3];
    double zact = pba->parameters_2_smg[4];
    double dNact = pba->parameters_2_smg[5];
    double N = log(a);
    double F0, am0, amp0;

    /* B0 supplies only the smooth time-profile basis. */
    sdmc_v2_b0_profile(pba,N,Ns_seed,zact,dNact,&F0,&am0,&amp0);

    double F = pow(F0,scale);
    double am = scale*am0;
    double ab = betaB*am;
    double ak = D0-1.5*ab*ab;

    pvecback[pba->index_bg_kineticity_smg] = ak;
    pvecback[pba->index_bg_braiding_smg] = ab;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = am;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = F-1.;
    pvecback[pba->index_bg_M2_smg] = F;
  }

'''
if model_block not in post:
    if post.count(model_anchor) != 1:
        raise RuntimeError(f"expected one alpha-model anchor after getter, found {post.count(model_anchor)}")
    post = post.replace(model_anchor, model_block + model_anchor)

p.write_text(pre + marker + post, encoding="utf-8")
print("SDMC NKp-v2 generalized EFT-completion patch complete.")
