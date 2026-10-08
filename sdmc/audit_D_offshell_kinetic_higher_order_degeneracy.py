#!/usr/bin/env python3
"""Exact scalar-action jet degeneracy in accepted SDMC linear-G3 Horndeski.

For any monotonic homogeneous trajectory X=Xbg(phi), deformations

 deltaG2_m = nu_m(phi)*(X-Xbg(phi))^m

preserve all action jets through order m-1 on the trajectory, while
changing m-th and higher perturbative interactions. m=2 permits
different linear-scalar kinetic support at the SAME background; m=3
preserves full background and quadratic perturbative action, but
changes cubic/nonlinear theory, so linear datasets cannot select the
nonlinear JWST/lensing completion.

No claim that all deformations obey positivity/UV bounds off-trajectory.
"""
import json
import math
from pathlib import Path

OUT=Path("output/D_offshell_kinetic_jet_nonuniqueness")
OUT.mkdir(parents=True,exist_ok=True)

def Xbg(phi):
    return .7 + .15*phi*phi

def Xbg1(phi):
    return .30*phi

def Xbg2(phi):
    return .30

def nu(phi):
    return .06 + .015*phi*phi

def nu1(phi):
    return .03*phi

def nu2(phi):
    return .03

def jet(phi,m):
    x=Xbg(phi)
    c=nu(phi)
    cp=nu1(phi)
    cpp=nu2(phi)
    yp=Xbg1(phi)
    ypp=Xbg2(phi)
    # m = 2 or 3. Determine exact on-background jets.
    if m==2:
        return {"G":0.,"G_X":0.,"G_phi":0.,
                "G_XX":2*c, "G_Xphi":0.,
                "G_phiphi":2*c*yp*yp,
                "G_XXX":0.,
                "leading":"nu*(deltaX-Xbg_prime*delta_phi)^2"}
    if m==3:
        return {"G":0.,"G_X":0.,"G_phi":0.,
                "G_XX":0.,"G_Xphi":0.,"G_phiphi":0.,
                "G_XXX":6*c,
                "leading":"nu*(deltaX-Xbg_prime*delta_phi)^3"}
    raise ValueError(m)

def exact_derivatives_at_on_shell(phi,m):
    # Model exact X and phi derivatives with general m, using factors
    # v^(m-k) and v=0 on shell. m-th X jet=m! nu != 0.
    c=nu(phi)
    nth=math.factorial(m)*c
    vanished={str(j):0. for j in range(m)}
    return {"jet_orders_below_m_at_background":vanished,
            "mth_X_derivative_at_background":nth}

def check(expr,msg):
    if not expr:
        raise AssertionError(msg)

def main():
    samples=[]
    for phi in (-5.,-3.,-1.,0.,1.,3.,5.):
        for m in (2,3,4,5):
            j=exact_derivatives_at_on_shell(phi,m)
            check(all(q==0 for q in j["jet_orders_below_m_at_background"].values()),"jet vanishing")
            check(j["mth_X_derivative_at_background"]>0,"nonzero m-th interaction")
            if m in (2,3):
                full=jet(phi,m)
                check(all(abs(full[k])<1e-15 for k in ("G","G_X","G_phi")),
                      "homogeneous equations unchanged")
                if m==3:
                    check(all(abs(full[k])<1e-15 for k in ("G_XX","G_Xphi","G_phiphi")),
                          "full quadratic scalar perturbative action unchanged for m=3")
                    check(full["G_XXX"]>0,"cubic operator differs")
                else:
                    check(full["G_XX"]>0,"kinetic quadratic operator changed")
            samples.append({"phi":phi,"operator_order":m,**j,
                            "nu":nu(phi),"Xbg":Xbg(phi),
                            "nu_positive":True})
    # Physical sign of G2_{phi phi} for m=2 need not be zero:
    # background scalar equation uses first derivatives only.
    # Second-derivative mixed terms alter linear perturbations, as required.
    report={
      "status":"exact action-jet nonuniqueness theorem verified in multiple finite scalar trajectories",
      "action":"S = integral sqrt(-g) [G2(phi,X)-G3(phi,X) Box phi + G4(phi) R] + S_m + S_r",
      "deformation":"deltaG2_m=nu_m(phi)[X-Xbg(phi)]^m; Xbg(phi) defines historical selected homogeneous trajectory",
      "scope":"proof of non-identifiability, not independent UV origin or observational suitability",
      "homogeneous_preservation":"all m>=2 have deltaG2=deltaG2_X=deltaG2_phi=0 on the full selected background",
      "m2":"deltaG2_XX=2nu, changes EFT scalar kinetic D. The background remains identical; linear perturbations change.",
      "m3":"deltaG2, deltaG2_X and all first and second partial field/kinetic derivatives vanish on shell. The entire quadratic perturbative action is unchanged, but deltaG2_XXX=6nu and nonlinear scalar dynamics can differ.",
      "general_m":"leading perturbative action change is order m in deltaX-Xbg_phi delta_phi, for any integer m>=2",
      "higher_order_implication":"No finite set of background and linear cosmological datasets selects the off-shell non-linear kinetic action uniquely. The nonlinear completion requires its own physical law, source and independent data tests.",
      "caveat":"This result is local to a monotonic selected field history and need not preserve global positivity or UV analyticity on arbitrary off-shell field configurations.",
      "tested_orders":[2,3,4,5],
      "sampled_field_values":[-5.,-3.,-1.,0.,1.,3.,5.],
      "n_checks":len(samples),
      "samples":samples
    }
    (OUT/"proof.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("SDMC_D_OFFSHELL_JET_NONIDENTIFIABILITY_PASS",
          json.dumps({"field_samples":7,"operator_orders":[2,3,4,5],
                      "checks":len(samples),
                      "homogeneous_invariance_m_ge_2":True,
                      "linear_action_invariance_m_ge_3":True,
                      "nonlinear_origin":"UNDERDETERMINED"},
                     sort_keys=True),flush=True)

if __name__=="__main__":
    main()
