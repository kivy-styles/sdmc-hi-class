#!/usr/bin/env python3
"""
Search a split-rate future stitch for the legacy canonical r=0 endpoint.

The one-rate canonical audit found a tradeoff:
  * slow common G2 release keeps c_s^2 <= 1 but produces a large transient
    Gamma=chi/F overshoot;
  * fast common release removes the Gamma overshoot but gives a small
    transient c_s^2 > 1.

This script asks whether that tradeoff is an artifact of forcing k1, k2 and V
to share one decay rate.  It keeps the accepted z>=0 action unchanged and the
mature endpoint fixed at
    r=0, chi_inf=F_inf, N_inf=sqrt(8*pi/3),
while allowing independent future-only rates.

The first scan uses mu_k1=mu_k2=muK and independent muV.  If no point satisfies
both Gamma<=1.0001 and c_s^2<=1.0001 while retaining D>0 and c_s^2>0, a small
second scan splits mu_k1 and mu_k2 around the best coarse point.

The best candidate is then re-run with two C3 No-Slip refinements.
"""
from pathlib import Path
import json,math
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/legacy_canonical_split_rate_audit.json")

R=0.0
MU_K=[3.753897593021392,5.0,10.0,20.0]
MU_V=[20.0,40.0,60.0,80.0,120.0]
TOL=1e-4

def evaluate(mu1,mu2,muv,nmax=6.,npts=1201,refine=False):
    spec={
      "name":f"canon_split_{mu1:g}_{mu2:g}_{muv:g}",
      "chi_inf":aud.FINF,
      "r":R,
      "muQ":mu1,
      "mu_k1":mu1,
      "mu_k2":mu2,
      "mu_V":muv,
    }
    try:
        model=(aud.build_refined_noslip_candidate(spec,iterations=2,blend_rate=500.)
               if refine else aud.build_candidate(spec))
        res=aud.evolve(model,Nmax=nmax,npts=npts)
        h=res["health_diagnostics"]
        gamma_max=float(h["Gamma_action_max"])
        csmax=float(h["max_cs2_horndeski"])
        csmin=float(h["min_cs2_horndeski"])
        Dmin=float(h["min_D_horndeski"])
        gamma_excess=max(0.,gamma_max-(1.+TOL))
        cs_super=max(0.,csmax-(1.+TOL))
        health_deficit=max(0.,-csmin)+max(0.,-Dmin)
        noslip=float(h["max_abs_alphaB_plus_2alphaM"])
        penalty=gamma_excess+10.*cs_super+100.*health_deficit+noslip
        feasible=(Dmin>0. and csmin>0.
                  and gamma_max<=1.+TOL and csmax<=1.+TOL)
        return {
          "status":"ok","mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "refined_noslip":refine,
          "Gamma_max":gamma_max,
          "Gamma_max_at_ln_a":float(h["Gamma_action_max_at_ln_a"]),
          "Gamma_final":float(h["Gamma_action_final"]),
          "min_D":Dmin,
          "min_cs2":csmin,
          "max_cs2":csmax,
          "max_abs_alphaB_plus_2alphaM":noslip,
          "q_max":float(res["future_kinematics"]["q_max"]),
          "N_final":float(res["final"]["N_struct"]),
          "p_final":float(res["final"]["p"]),
          "q_final":float(res["final"]["q"]),
          "penalty":penalty,
          "strict_feasible":feasible,
        }
    except Exception as e:
        return {
          "status":"failure","mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "refined_noslip":refine,"strict_feasible":False,
          "penalty":1e99,"error":str(e)
        }

coarse=[]
for muk in MU_K:
    for muv in MU_V:
        coarse.append(evaluate(muk,muk,muv))

best_coarse=min(coarse,key=lambda r:r["penalty"])
strict=[r for r in coarse if r.get("strict_feasible")]

split=[]
if not strict:
    base=float(best_coarse["mu_k1"])
    muv0=float(best_coarse["mu_V"])
    kgrid=sorted(set([max(.5,base/2.),base,min(120.,base*2.)]))
    vgrid=sorted(set([max(5.,muv0/1.5),muv0,min(160.,muv0*1.5)]))
    for mu1 in kgrid:
        for mu2 in kgrid:
            for muv in vgrid:
                split.append(evaluate(mu1,mu2,muv))
    strict=[r for r in split if r.get("strict_feasible")]

pool=[r for r in coarse+split if r["status"]=="ok"]
best=min(strict,key=lambda r:r["penalty"]) if strict else min(pool,key=lambda r:r["penalty"])
refined=evaluate(best["mu_k1"],best["mu_k2"],best["mu_V"],
                 nmax=10.,npts=2001,refine=True)

out={
  "status":(
    "future-only split-rate search for the legacy canonical r=0 endpoint. "
    "No accepted z>=0 coefficient is changed."
  ),
  "strict_limits":{
    "Gamma_max":1.+TOL,
    "cs2_max":1.+TOL,
    "D_min":0.0,
    "cs2_min":0.0,
  },
  "coarse_common_kinetic_scan":coarse,
  "split_followup_scan":split,
  "best_unrefined":best,
  "best_refined":refined,
  "strict_solution_found_before_refinement":bool(strict),
  "strict_solution_survives_refinement":bool(refined.get("strict_feasible",False)),
  "interpretation":(
    "If a strict solution survives refinement, the earlier one-rate tradeoff "
    "was interpolation-induced rather than a property of the r=0 endpoint. "
    "If not, the current C2/C3 tail family still cannot simultaneously enforce "
    "bounded covariant activation and subluminal scalar propagation."
  )
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("LEGACY_CANONICAL_SPLIT_RATE_AUDIT")
print("BEST_COARSE",json.dumps(best_coarse,sort_keys=True))
print("BEST_UNREFINED",json.dumps(best,sort_keys=True))
print("BEST_REFINED",json.dumps(refined,sort_keys=True))
print("STRICT_FOUND",bool(strict))
print("STRICT_REFINED",bool(refined.get("strict_feasible",False)))
