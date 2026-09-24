#!/usr/bin/env python3
"""
Stress-test bounded *covariant* activation for the mature matched endpoint.

The raw source activation is chi=rho_X/rho_v.  In the mature non-minimal
No-Slip theory the gravitationally weighted activation is

    Gamma = (G_eff rho_X)/(G rho_v) = chi/F.

If the Unified Balanced Identity is covariantized as a statement about
gravitational strength, the natural mature saturation condition is Gamma->1,
not chi->1.  This audit scans faster G2 transition rates on the
chi_inf=F_inf branch and asks whether the transient Gamma>1 overshoot can be
removed while D>0 and c_s^2>0 remain satisfied.

The accepted z>=0 action is unchanged.
"""
from pathlib import Path
import json,numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/future_covariant_activation_bound_stress_test.json")

R_FIXED=3.756118836402892
F_INF=aud.FINF
MU_GRID=[
    3.753897593021392,
    10.0,20.0,40.0,60.0,70.0,76.0,80.0,90.0,100.0,120.0
]

rows=[]
for mu in MU_GRID:
    spec=dict(
      name=f"matched_gamma_mu_{mu:g}",
      chi_inf=F_INF,
      r=R_FIXED,
      muQ=mu
    )
    model=aud.build_candidate(spec)
    NN,yy,diag=aud.release_trajectory(model,Nmax=10.,npts=2001)
    result=aud.evolve(model,Nmax=10.,npts=2001)

    # release_trajectory returns raw rhs diagnostics; reconstruct the positive
    # structural activation exactly as evolve() does:
    #   chi/chi0=(rhoX_struct/rhoX_struct,0)*sigma^2.
    rhoXa=np.asarray([x["rhoX_structural"] for x in diag])
    chi=aud.CHI0*(rhoXa/rhoXa[0])*yy[0]*yy[0]
    gamma=chi/np.asarray([x["F"] for x in diag])
    imax=int(np.argmax(gamma))
    h=result["health_diagnostics"]
    f=result["future_kinematics"]
    row={
      "muQ":mu,
      "Gamma_max":float(gamma[imax]),
      "Gamma_max_at_ln_a":float(NN[imax]),
      "Gamma_final":float(gamma[-1]),
      "transient_Gamma_overshoot_over_final":
        float(gamma[imax]/max(gamma[-1],1e-300)-1.),
      "transient_Gamma_overshoot_over_unit":
        float(gamma[imax]-1.),
      "min_D":h["min_D_horndeski"],
      "min_cs2":h["min_cs2_horndeski"],
      "max_cs2":h["max_cs2_horndeski"],
      "max_abs_alphaB_plus_2alphaM":
        h["max_abs_alphaB_plus_2alphaM"],
      "mapper_ratio_final":h["mapper_ratio_final"],
      "q_max":f["q_max"],
      "q_zero_crossings":f["q_zero_crossings"],
      "N_struct_final":result["final"]["N_struct"],
      "p_final":result["final"]["p"],
      "q_final":result["final"]["q"],
    }
    rows.append(row)

first_no_transient=None
for r in rows:
    # Compare against the late reconstructed Gamma floor rather than demanding
    # exact unity at finite ln a=10.
    if r["transient_Gamma_overshoot_over_final"] <= 1e-4:
        first_no_transient=r
        break

out={
  "status":(
    "covariant-activation stress test on the chi_inf=F_inf branch; "
    "Gamma=chi/F is the gravitationally weighted activation and tends to 1."
  ),
  "fixed_r":R_FIXED,
  "F_inf":F_INF,
  "rows":rows,
  "first_scanned_without_resolved_transient_overshoot":first_no_transient,
  "interpretation":(
    "A slow future G2 handoff can transiently drive Gamma above unity. "
    "The scan tests whether a faster future-only coefficient flow removes that "
    "overshoot without violating the exact Horndeski kinetic and sound-speed "
    "health conditions. This is diagnostic and does not promote a tail."
  )
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_COVARIANT_ACTIVATION_BOUND_STRESS_TEST")
for r in rows:
    print("GAMMA_SCAN",json.dumps(r,sort_keys=True))
print("FIRST_NO_TRANSIENT",json.dumps(first_no_transient,sort_keys=True))
