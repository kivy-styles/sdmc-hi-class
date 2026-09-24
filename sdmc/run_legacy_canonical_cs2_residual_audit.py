#!/usr/bin/env python3
"""
Decompose the tiny c_s^2-1 residual of the best fine-search legacy-canonical
future stitch.

The mature r=0 normal form has c_s^2==1 exactly for every Z, so any excess
above unity before the asymptotic regime must come from the finite transition
sector (running F, residual G3, or incomplete G2 coefficient release).

This audit evaluates the exact Bellini-Sawicki numerator term by term for the
best dense fine-search point (mu_k1,mu_k2,mu_V)=(46,7,100).
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/legacy_canonical_cs2_residual_audit.json")
SPEC={
  "name":"canonical_cs2_residual",
  "chi_inf":aud.FINF,
  "r":0.0,
  "muQ":46.0,
  "mu_k1":46.0,
  "mu_k2":7.0,
  "mu_V":100.0,
}
NMAX=4.0
NPTS=8001

model=aud.build_refined_noslip_candidate(SPEC,iterations=2,blend_rate=500.)
NN,yy,diag=aud.release_trajectory(model,Nmax=NMAX,npts=NPTS)

rows=[]
for ne,(sig,v),d in zip(NN,yy.T,diag):
    F=float(d["F"]); H=float(d["H"])
    b=float(d["alphaB"]); m=float(d["alphaM"])
    D=float(d["D_full"])
    xp=float(d["rho_smg_class"]+d["p_smg_class"])
    mp=float((d["rho_m_action"]+4.*d["rho_r_action"]/3.)/3.)
    t_noslip=((2.-b)*(b+2.*m)/2.)
    t_scalar=1.5*(2.-b)*xp/(H*H)
    t_matter=-1.5*(2.-2.*F+b*F)*mp/(H*H*F)
    t_deriv=float(d["alphaB_N"])
    num=t_noslip+t_scalar+t_matter+t_deriv
    cs2=num/D

    aq=model["action"](float(sig))
    Z=.5*float(v)*float(v)
    Zstar=float(model["Z_inf"])
    pure_r0_k1=2.*aud.FINF
    pure_r0_V=4.*aud.FINF*Zstar
    shape={
      "F_minus_Finf":float(aq["F"]-aud.FINF),
      "alphaM":m,
      "alphaB":b,
      "alphaB_plus_2alphaM":float(b+2.*m),
      "g":float(aq["g"]),
      "g_times_Z":float(aq["g"]*Z),
      "k1_over_target_minus_1":float(aq["k1"]/pure_r0_k1-1.),
      "k2_Zstar_over_Finf":float(aq["k2"]*Zstar/aud.FINF),
      "V_over_target_minus_1":float(aq["V"]/pure_r0_V-1.),
      "Z_over_Zstar_minus_1":float(Z/Zstar-1.),
    }
    rows.append({
      "ln_a":float(ne),
      "sigma":float(sig),
      "D":D,
      "cs2":float(cs2),
      "cs2_minus_1":float(cs2-1.),
      "numerator":float(num),
      "numerator_minus_D":float(num-D),
      "terms":{
        "noslip_combo":float(t_noslip),
        "scalar_effective_rho_plus_p":float(t_scalar),
        "matter":float(t_matter),
        "alphaB_N":float(t_deriv),
      },
      "kessence_proxy":float(d["cs2_kessence_proxy"]),
      "transition_shape":shape,
    })

imax=max(range(len(rows)),key=lambda i:rows[i]["cs2"])
maxrow=rows[imax]

def nearest(target):
    i=min(range(len(rows)),key=lambda j:abs(rows[j]["ln_a"]-target))
    return rows[i]

samples=[nearest(x) for x in [0.,.5,1.,2.,2.5,2.835,3.,3.5,4.]]

out={
  "status":(
    "exact term decomposition of the residual sound-speed excess on the "
    "best fine-search r=0 canonical transition"
  ),
  "spec":SPEC,
  "analytic_endpoint":{
    "r":0.0,
    "pure_mature_cs2":"1 exactly for every Z",
    "interpretation":(
      "therefore any finite cs2-1 in this audit is a transition residual, "
      "not an intrinsic property of the mature canonical action"
    )
  },
  "max_cs2":maxrow,
  "samples":samples,
  "diagnostic":{
    "max_excess":float(maxrow["cs2_minus_1"]),
    "max_at_ln_a":float(maxrow["ln_a"]),
    "max_numerator_minus_D":float(maxrow["numerator_minus_D"]),
    "dominant_absolute_numerator_term":
      max(maxrow["terms"],key=lambda k:abs(maxrow["terms"][k])),
  }
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("LEGACY_CANONICAL_CS2_RESIDUAL_AUDIT")
print("MAX",json.dumps(maxrow,sort_keys=True))
for row in samples:
    print("SAMPLE",json.dumps(row,sort_keys=True))
