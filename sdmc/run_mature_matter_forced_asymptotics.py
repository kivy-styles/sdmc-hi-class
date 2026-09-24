#!/usr/bin/env python3
"""
Matter-forced asymptotics of the mature inverse-square SDMC scalar.

The fluid-level diagnostic
    N/N_inf = sqrt(1+rho_m/rho_X+...)
is not the exact action-level response of a healthy inverse-square scalar.
For
    G2 = sigma^-2 f(Z),   G3 -> 0,   F -> const,
the scalar equation and Friedmann constraint can be linearized around the
pure-scalar coasting fixed point.

Define
    kappa = K_*/(2 F) = 1+2r = 1/c_s,*^2,
    z = Z/Z_* - 1,
    u = p-1,
    m = rho_m sigma^2 / R_*,
    R_* = 6 F Z_*.

To first order in late matter loading,
    z_N + 2 z = -(3/kappa) m,   m_N=-m,
so
    z = -(3/kappa)m + C a^-2,
    u = -(3/(2 kappa))m + O(a^-2),
    N/N_inf - 1 = -(3/(2 kappa))m + O(a^-2),
    q = -(3/(2 kappa))m + O(a^-2).

Thus the intrinsic scalar perturbation still decays as a^-2, but the forced
matter response decays only as a^-1 and eventually dominates.  The mature
action approaches coasting from the very weakly accelerating side, q->0^-,
not from the decelerating side of the frozen-fluid heuristic.

This script verifies these asymptotic relations against the freely evolved
C3 future tails.
"""
from pathlib import Path
import json, math
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/mature_matter_forced_asymptotics.json")
CHECK_N=[8.0,10.0]

out={
  "status":(
    "linearized action-level asymptotic audit around the mature inverse-square "
    "fixed point; it supersedes the sign-of-q interpretation of the earlier "
    "frozen-fluid matter-loading diagnostic, while retaining its 2.868923 "
    "value only as a diagnostic interior lapse reference."
  ),
  "analytic":{
    "kappa":"Kstar/(2F)=1+2r=1/cs2_star",
    "forced_equation":"dz/dln(a)+2z=-(3/kappa)m, dm/dln(a)=-m",
    "late_solution":"z=-(3/kappa)m+C*a^-2",
    "p_minus_1":"-(3/(2kappa))*m",
    "N_over_Ninf_minus_1":"-(3/(2kappa))*m",
    "q":"-(3/(2kappa))*m",
    "rhoX_rescaled_minus_1":"-m",
    "interpretation":(
      "matter forces the scalar rescaled density downward at first order; "
      "the scalar response cancels the direct matter contribution to H^2 "
      "at O(m), leaving H*sigma fixed to first order and q->0 from below."
    )
  },
  "candidates":{}
}

for spec in aud.CANDIDATES:
    model=aud.build_candidate(spec)
    result=aud.evolve(model,Nmax=10.,npts=2001)
    r=float(result["tail_parameters"]["r"])
    kappa=1.+2.*r
    Zinf=aud.Z0*float(result["asymptotic_target"]["Z_inf_over_Z0"])
    Rstar=6.*aud.FINF*Zinf

    rows=[]
    for target in CHECK_N:
        s=min(result["samples"],key=lambda x:abs(x["ln_a"]-target))
        m=float(s["rho_m_action"])*float(s["sigma"])**2/Rstar
        z=float(s["Z_ratio"])/float(result["asymptotic_target"]["Z_inf_over_Z0"])-1.
        u=float(s["p"])-1.
        nd=float(s["N_struct"])/float(result["asymptotic_target"]["N_inf"])-1.
        q=float(s["q"])
        pred_z=-3.*m/kappa
        pred_common=-1.5*m/kappa

        def rel(obs,pred):
            return (obs-pred)/pred if abs(pred)>1e-300 else None

        rows.append({
          "ln_a":float(s["ln_a"]),
          "a":float(s["a"]),
          "matter_loading_m":m,
          "z_observed":z,
          "z_predicted":pred_z,
          "z_fractional_residual":rel(z,pred_z),
          "p_minus_1_observed":u,
          "N_over_Ninf_minus_1_observed":nd,
          "q_observed":q,
          "common_predicted":pred_common,
          "p_fractional_residual":rel(u,pred_common),
          "N_fractional_residual":rel(nd,pred_common),
          "q_fractional_residual":rel(q,pred_common),
        })

    out["candidates"][spec["name"]]={
      "r":r,
      "kappa":kappa,
      "cs2_star":1./kappa,
      "N_inf":float(result["asymptotic_target"]["N_inf"]),
      "rows":rows,
    }

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("MATURE_MATTER_FORCED_ASYMPTOTICS")
for name,data in out["candidates"].items():
    print("CANDIDATE",name,"kappa",data["kappa"])
    for row in data["rows"]:
        print("FORCED_ASYMPTOTIC",name,json.dumps(row,sort_keys=True))
