#!/usr/bin/env python3
"""
Fine search for a genuinely subluminal refined legacy-canonical future stitch.

The focused split-rate search found six healthy refined r=0 solutions with
Gamma<=1.0001 and max(c_s^2)<=1.0001.  The best point,
(mu_k1,mu_k2,mu_V)=(44,8,110), missed exact c_s^2<=1 only by ~6.5e-6.

The max c_s^2 response is monotonic enough in the local mu_k2 direction that
a narrow scan below mu_k2=8 can test whether the residual is only a stitching
artifact.  We impose the stronger numerical ceilings

    Gamma_max <= 1 + 1e-8
    max(c_s^2) <= 1 + 1e-8

after two C3 No-Slip refinement iterations.

Accepted z>=0 coefficients are unchanged.  The mature target remains
r=0, chi_inf=F_inf, N_inf=sqrt(8*pi/3).
"""
from pathlib import Path
import json,itertools
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/legacy_canonical_subluminal_fine_search.json")

MU1=[42.0,44.0,46.0]
MU2=[7.0,7.2,7.4,7.6,7.8,8.0]
MUV=[100.0,110.0,120.0]
TOL=1e-8

def evaluate(mu1,mu2,muv,nmax=10.,npts=3001):
    spec={
      "name":f"canon_sublum_{mu1:g}_{mu2:g}_{muv:g}",
      "chi_inf":aud.FINF,
      "r":0.0,
      "muQ":mu1,
      "mu_k1":mu1,
      "mu_k2":mu2,
      "mu_V":muv,
    }
    try:
        model=aud.build_refined_noslip_candidate(
            spec,iterations=2,blend_rate=500.)
        res=aud.evolve(model,Nmax=nmax,npts=npts)
        h=res["health_diagnostics"]
        gm=float(h["Gamma_action_max"])
        cmin=float(h["min_cs2_horndeski"])
        cmax=float(h["max_cs2_horndeski"])
        dmin=float(h["min_D_horndeski"])
        ns=float(h["max_abs_alphaB_plus_2alphaM"])
        strict=(gm<=1.+TOL and cmax<=1.+TOL and cmin>0. and dmin>0.)
        # Rank by exact ceiling violations first, then No-Slip and q transient.
        gex=max(0.,gm-(1.+TOL))
        cex=max(0.,cmax-(1.+TOL))
        score=1e5*gex+1e6*cex+ns+1e-4*abs(float(res["future_kinematics"]["q_max"]))
        return {
          "status":"ok",
          "mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "Gamma_max":gm,
          "Gamma_max_at_ln_a":float(h["Gamma_action_max_at_ln_a"]),
          "Gamma_final":float(h["Gamma_action_final"]),
          "min_D":dmin,
          "min_cs2":cmin,
          "min_cs2_at_ln_a":float(h["min_cs2_horndeski_at_ln_a"]),
          "max_cs2":cmax,
          "max_cs2_at_ln_a":float(h["max_cs2_horndeski_at_ln_a"]),
          "max_abs_alphaB_plus_2alphaM":ns,
          "q_max":float(res["future_kinematics"]["q_max"]),
          "N_final":float(res["final"]["N_struct"]),
          "p_final":float(res["final"]["p"]),
          "q_final":float(res["final"]["q"]),
          "strict_subluminal":strict,
          "score":score,
        }
    except Exception as e:
        return {
          "status":"failure","mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "strict_subluminal":False,"score":1e99,"error":str(e)
        }

rows=[evaluate(*x) for x in itertools.product(MU1,MU2,MUV)]
ok=[r for r in rows if r["status"]=="ok"]
strict=[r for r in ok if r["strict_subluminal"]]
best=(min(strict,key=lambda r:r["score"]) if strict
      else min(ok,key=lambda r:r["score"]))

# Re-evaluate the selected point on a denser future grid.
validated=evaluate(best["mu_k1"],best["mu_k2"],best["mu_V"],
                   nmax=10.,npts=8001)

out={
  "status":(
    "fine two-iteration No-Slip search for exact-sub-luminal legacy-canonical "
    "future stitch; accepted z>=0 action unchanged"
  ),
  "limits":{
    "Gamma_max":1.+TOL,
    "cs2_max":1.+TOL,
    "D_min":0.0,
    "cs2_min":0.0
  },
  "grid":{"mu_k1":MU1,"mu_k2":MU2,"mu_V":MUV},
  "rows":rows,
  "strict_count":len(strict),
  "best_scan":best,
  "best_dense_validation":validated,
  "strict_dense_validation":bool(validated.get("strict_subluminal",False)),
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("LEGACY_CANONICAL_SUBLUMINAL_FINE_SEARCH")
print("STRICT_COUNT",len(strict))
print("BEST_SCAN",json.dumps(best,sort_keys=True))
print("DENSE_VALIDATION",json.dumps(validated,sort_keys=True))
print("STRICT_DENSE",bool(validated.get("strict_subluminal",False)))
