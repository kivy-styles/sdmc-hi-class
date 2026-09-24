#!/usr/bin/env python3
"""
Focused refinement of the split-rate r=0 canonical future stitch.

The first split-rate audit found an unrefined strict point near
(mu_k1,mu_k2,mu_V)=(40,10,120), while the two-step C3 No-Slip refinement
missed the cs2<=1.0001 screen by only ~9e-6.

This audit searches the local neighborhood with a tighter target:
  Gamma_max <= 1 + 1e-6
  c_s^2_max <= 1 + 1e-6
  D > 0
  c_s^2 > 0

Unrefined candidates are ranked first.  The best few are then re-run with two
C3 No-Slip refinements and full ln a=10 evolution.  The accepted z>=0 action
is unchanged.
"""
from pathlib import Path
import json
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/legacy_canonical_split_rate_refinement.json")

K1_GRID=[35.0,40.0,45.0,50.0]
K2_GRID=[5.0,7.5,10.0]
V_GRID=[120.0,140.0,160.0]
TOL=1e-6
NREFINE=8

def evaluate(mu1,mu2,muv,refine=False,nmax=6.,npts=1201):
    spec={
      "name":f"canon_local_{mu1:g}_{mu2:g}_{muv:g}",
      "chi_inf":aud.FINF,
      "r":0.0,
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
        gm=float(h["Gamma_action_max"])
        cmax=float(h["max_cs2_horndeski"])
        cmin=float(h["min_cs2_horndeski"])
        dmin=float(h["min_D_horndeski"])
        ns=float(h["max_abs_alphaB_plus_2alphaM"])
        gex=max(0.,gm-(1.+TOL))
        cex=max(0.,cmax-(1.+TOL))
        deficit=max(0.,-cmin)+max(0.,-dmin)
        # No-Slip residual is a secondary smoothness discriminator, not a hard
        # endpoint condition here because the iterative refinement already
        # drives the physical combination very small.
        penalty=gex+10.*cex+100.*deficit+0.05*ns
        strict=(gm<=1.+TOL and cmax<=1.+TOL and cmin>0. and dmin>0.)
        return {
          "status":"ok","mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "refined_noslip":refine,
          "Gamma_max":gm,
          "Gamma_final":float(h["Gamma_action_final"]),
          "max_cs2":cmax,"min_cs2":cmin,"min_D":dmin,
          "max_abs_alphaB_plus_2alphaM":ns,
          "q_max":float(res["future_kinematics"]["q_max"]),
          "N_final":float(res["final"]["N_struct"]),
          "p_final":float(res["final"]["p"]),
          "q_final":float(res["final"]["q"]),
          "strict_feasible":strict,
          "penalty":penalty,
        }
    except Exception as e:
        return {
          "status":"failure","mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "refined_noslip":refine,"strict_feasible":False,
          "penalty":1e99,"error":str(e)
        }

coarse=[]
for mu1 in K1_GRID:
    for mu2 in K2_GRID:
        for muv in V_GRID:
            coarse.append(evaluate(mu1,mu2,muv))

ranked=sorted([r for r in coarse if r["status"]=="ok"],
              key=lambda r:r["penalty"])
selected=ranked[:NREFINE]
refined=[
  evaluate(r["mu_k1"],r["mu_k2"],r["mu_V"],
           refine=True,nmax=10.,npts=2001)
  for r in selected
]
refined_ok=[r for r in refined if r["status"]=="ok"]
strict=[r for r in refined_ok if r["strict_feasible"]]
best=(min(strict,key=lambda r:r["penalty"])
      if strict else min(refined_ok,key=lambda r:r["penalty"]))

out={
  "status":(
    "focused local split-rate refinement for the r=0 canonical endpoint; "
    "all modifications are future-only."
  ),
  "strict_tolerance":TOL,
  "coarse_grid":coarse,
  "selected_for_refinement":selected,
  "refined_results":refined,
  "strict_refined_solutions":strict,
  "strict_refined_solution_found":bool(strict),
  "best_refined":best,
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("LEGACY_CANONICAL_SPLIT_RATE_REFINEMENT")
print("SELECTED",json.dumps(selected,sort_keys=True))
for r in refined:
    print("REFINED_LOCAL",json.dumps(r,sort_keys=True))
print("STRICT_REFINED_FOUND",bool(strict))
print("BEST_REFINED_LOCAL",json.dumps(best,sort_keys=True))
