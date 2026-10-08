#!/usr/bin/env python3
"""Simultaneous mathematical-closure and non-identifiability audit for late300.

Only the accepted SDMC ansatz, reconstructed action algebra, already reported
GitHub results, and standard definitions are tested here. An algebraic
identity or inverse-reconstruction check is *not* a microscopic derivation.
No likelihood optimization is performed. This does not alter late300.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

PARAMS = json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]
ORPH = 4.17998772e-5
TAU_A, TAU_B = 0.25, 1.5
EARLY_HANDOFF_WIDTH = 0.5
DERIVED = {
    "AF": 0.020520,
    "zc": 4.034077502899133,
    "wF": 0.32925377073762596,
    "Dfloor": 0.05181542627513409,
    "D0": 0.3419932122356156,
    "pD": 1.0,
    "lambda_e": 18.40625,
    "zt": 16.742289660253434,
    "A_late": 0.024822032167576252,
    "B_late": 0.01264704344701022,
}
HSC_DELTA_CHI2 = -2.6668948
HSC_N = 60
PLANCK_DLOGZ, PLANCK_ERR = -2.12549, 0.69167
HSC_DLOGZ, HSC_ERR = 2.0390446846788706, 0.6148357940108409


def need(condition: bool, msg: str):
    if not condition:
        raise AssertionError(msg)


def close(x: float, y: float, tol: float = 1e-11):
    return math.isclose(x, y, abs_tol=tol, rel_tol=tol)


def s_step(N: float, zc: float, w: float) -> float:
    Nc = -math.log1p(zc)
    x = (N-Nc)/w
    return 1./(1.+math.exp(-x)) if x < 700. else 1.


def f_of_N(N: float, af: float, zc: float, width: float) -> float:
    return math.exp(af*s_step(N,zc,width))


def alpha_M(N: float, af: float, zc: float, width: float) -> float:
    s=s_step(N,zc,width)
    return af*s*(1.-s)/width


def d_support(N: float, floor: float, amp: float, power: float,
              zc: float, width: float) -> float:
    return floor + amp*s_step(N,zc,width)**power


def tracker_fraction(lam: float, radiation: bool) -> float:
    return (4. if radiation else 3.)/lam**2


def omega_background():
    p=PARAMS
    return 1.-(p["omega_b"]+p["omega_cdm"]+ORPH)/(p["H0"]/100.)**2


def W_hand(z: float, zt: float, width: float) -> float:
    N=-math.log1p(z)
    Nt=-math.log1p(zt)
    return .5*(1.-math.tanh((N-Nt)/width))


def delta(z: float, A: float, B: float) -> float:
    return A*z*math.exp(-z/TAU_A)-B*z*z*math.exp(-z/TAU_B)


def solve_late_from_two_residuals(z1: float, z2: float, y1: float, y2: float):
    g1=z1*math.exp(-z1/TAU_A)
    h1=-z1*z1*math.exp(-z1/TAU_B)
    g2=z2*math.exp(-z2/TAU_A)
    h2=-z2*z2*math.exp(-z2/TAU_B)
    det=g1*h2-g2*h1
    need(abs(det)>1e-14,"redshift anchors do not independently constrain A,B")
    return ((y1*h2-y2*h1)/det,(g1*y2-g2*y1)/det)


def covariant_kinetic_inversion(C: float, F: float, alphaK: float,
                               gp: float, X: float, g: float, H2: float):
    need(F>0 and X>0,"reconstruction domain")
    k2=(F*alphaK-C+4*gp*X-6*g*H2)/(4*X)
    k1=C-2*k2*X
    predicted_alphaK=(k1+6*k2*X-4*gp*X+6*g*H2)/F
    return k1,k2,predicted_alphaK


def audit():
    p=PARAMS
    tests=[]
    d={}
    # Exact flat present closure, not a new adjustable SDMC parameter.
    o=omega_background()
    need(close(o,p["Omega_x0"],1e-13),"Omega_x0 flat closure failed")
    d["Omega_x0"]={"status":"exact algebraic closure given H0, physical densities and flatness",
        "derived":o,"accepted":p["Omega_x0"],"independent_sdmc_parameter":False}
    tests.append("Omega_x0 flat closure")

    # Planck-mass trajectory: derive coordinate *identifiability* exactly.
    AF,zc,w=DERIVED["AF"],DERIVED["zc"],DERIVED["wF"]
    Nc=-math.log1p(zc)
    sm=s_step(Nc,zc,w)
    F0=f_of_N(Nc,AF,zc,w)
    peak=alpha_M(Nc,AF,zc,w)
    need(close(sm,.5,1e-14),"F transition centre")
    need(close(math.log(F0),AF/2.,1e-12),"F centre normalization")
    need(close(peak,AF/(4*w),1e-13),"F alphaM peak")
    for n in [Nc-1.,Nc,Nc+1.]:
        e=1e-5
        derivative=(math.log(f_of_N(n+e,AF,zc,w))-math.log(f_of_N(n-e,AF,zc,w)))/(2*e)
        need(close(derivative,alpha_M(n,AF,zc,w),2e-10),"alphaM=dlogF/dN")
    need(close(math.log(f_of_N(50,AF,zc,w)/f_of_N(-50,AF,zc,w)),AF,1e-12),
         "Planck-mass integrated amplitude")
    d["F_coordinates"]={
        "status":"action-trajectory identities plus force/stability conditional candidates; not uniquely predicted",
        "candidate":{"AF":AF,"zc":zc,"wF":w},
        "identities":{"lnF_at_centre":math.log(F0),"alphaM_peak":peak,
          "Ncentre":Nc,"delta_lnF_full_transition":AF},
        "stability_result_at_derived_handoff":{
          "AF_zero_approx":0.020517949990009835,
          "zero_is_lower_bound_not_equality":True},
        "missing":"Independent microscopic F(S), forcing boundary conditions, and selection of amplitude above health threshold"}
    tests.append("F transition identities and stability-bound classification")

    # Kinetic trajectory has arbitrary positive amplitude, floor and exponent
    # until a genuinely independent k1,k2,g,F evolution is imposed.
    d0=DERIVED["D0"]; floor=DERIVED["Dfloor"]
    center_N=-math.log1p(p["z_c"])
    for po in [0.5,1.,1.5,2.]:
        vals=[d_support(n,floor,d0,po,p["z_c"],p["width"]) for n in [-20.,center_N,0.,20.]]
        need(all(v>0 for v in vals),f"healthy positive D at p={po}")
    e=1e-9
    accepted_D= d_support(0.,floor,d0,1.,p["z_c"],p["width"])
    alphaM_today=alpha_M(0.,AF,zc,w)
    alphaB_today=-2.*alphaM_today
    alphaK_today=accepted_D-1.5*alphaB_today**2
    need(close(accepted_D,alphaK_today+1.5*alphaB_today**2,1e-14),"kinetic algebra")
    d["D_coordinates"]={"status":"conditionally extractable from independently known D(N), not dynamically selected",
        "candidate":{"Dfloor":floor,"D0":d0,"pD":1},
        "positive_D_for_alternative_powers":[.5,1.,1.5,2.],
        "identity":"D=alphaK+3*alphaB^2/2=alphaK+6*alphaM^2 on No-Slip background",
        "alphaK_reconstructed_from_D_at_z0":alphaK_today,
        "D_z0":accepted_D,
        "missing":"Independent kinetic coefficient functions rather than back-solving them from D"}
    tests.append("positive kinetic family with non-unique exponent")

    # The accepted covariant reconstruction is an explicit inverse map:
    # for fixed F, H, etc. *any* D target gives a different (k1,k2)
    # that exactly matches alphaK. This algebraic test does not assert
    # stability or matching observations for all alternatives.
    C,F,g,gp,X,H2,alphaB=1.4,1.025,.08,.015,.6,.04,-.006
    inversion_rows=[]
    for Dt in [.04723,.051815,.2,.35]:
        ak=Dt-1.5*alphaB**2
        k1,k2,akr=covariant_kinetic_inversion(C,F,ak,gp,X,g,H2)
        need(close(ak,akr,1e-12),"covariant inverse kinetic closure")
        inversion_rows.append({"prescribed_D":Dt,"recovered_alphaK":akr,
                               "k1":k1,"k2":k2})
    need(len(set(round(a["k2"],10) for a in inversion_rows))==len(inversion_rows),
         "non-uniqueness of inverse kinetic coefficient choice")
    d["covariant_inverse_reconstruction"]={
       "status":"verified algebraically: reproducing selected D after reconstruction does not independently select it",
       "source_equations":"sdmc/run_late300_linear_covariant.py k2,k1,alphaK identities",
       "illustrative_dimensionless_control_variables":{"C":C,"F":F,"g":g,"gp":gp,"X":X,"H2":H2,"alphaB":alphaB},
       "alternative_exact_reconstructions":inversion_rows,
       "limitation":"Algebraic reconstruction, not 4 different healthy freely evolving cosmologies; global dynamics and data must be separately tested."}
    tests.append("covariant inverse-reconstruction non-identifiability certificate")

    # Tracker: the accepted early structural fraction is exactly related
    # to the supplied lambda; it is not a selector of lambda.
    lambda_in=p["lambda_e"]
    fr=tracker_fraction(lambda_in,True)
    fm=tracker_fraction(lambda_in,False)
    need(close(lambda_in,math.sqrt(4/fr),1e-12),"tracker tautological inverse")
    alternatives={str(x):tracker_fraction(x,True) for x in [14.,16.,18.,18.40625,20.,22.]}
    need(len(set(round(v,10) for v in alternatives.values()))==len(alternatives),
         "tracker slope remains unconstrained")
    dcan=4*fr
    rk=(floor/dcan-1)/(6-2*(floor/dcan))
    need(close((1+6*rk)/(1+2*rk),floor/dcan),"conditional kessence floor map")
    d["tracker_and_rK"]={
        "status":"fr=4/lambda^2 identity only; k-essence inferred rK requires accepted floor",
        "lambda_input":lambda_in,"fr_inferred":fr,"fm_inferred":fm,
        "lambda_alternatives_and_fr":alternatives,
        "D_radiation_canonical":dcan,
        "rK_conditional_inferred":rk,
        "reconstruction_observation":"k1 changes sign in late300 linear-G3 reconstructed action near z ~ 24.5",
        "missing":"Independent early fraction, UV kinetic coefficient ratio, activation/transfer evolution"}
    tests.append("lambda non-circularity and kessence rK qualification")

    # Handoff: a root of a handoff-dependent curvature map is a conditional
    # fixed-point *equation*, not a prediction without specified map.
    zt=DERIVED["zt"]
    w_at_handoff=W_hand(zt,zt,EARLY_HANDOFF_WIDTH)
    need(close(w_at_handoff,.5,1e-14),"handoff midpoint identity")
    d["zt"]={"status":"conditional self-consistent fixed-point candidate",
         "fixedpoint_from_prior_run":zt,
         "W_at_input_midpoint":w_at_handoff,
         "historical_fitted":p["z_t"],
         "note":"W(z_t)=1/2 follows from the definition for any z_t. The separate curvature-root map still requires independent physical inputs and uniqueness checks.",
         "missing":"UV/activation-curvature dynamics independent of the already chosen handoff ansatz"}
    tests.append("handoff midpoint tautology is not a selector")

    # Re-read completed 8-helper curvature-root job diagnostics from
    # github.com/kivy-styles/sdmc-hi-class/actions/runs/37763535074.
    # A given helper produces several curvature-derivative zeros.
    # A near-match must not be mistaken for a unique microscopic transition.
    handoff_scan=[
      (16.738,16.74298834420653),
      (16.74,16.741331979409132),
      (16.742,16.741843194390437),
      (16.7422,16.74217300780501),
      (16.7423,16.743144040965365),
      (16.7424,16.74275794986044),
      (16.744,16.74272969535322),
      (16.746,16.742235841754663),
    ]
    roots=[item[1] for item in handoff_scan]
    resids=[item[1]-item[0] for item in handoff_scan]
    sign_reversals=sum((resids[j]>=0)!=(resids[j-1]>=0) for j in range(1,len(resids)))
    need(len(handoff_scan)==8 and sign_reversals>=2, "curvature-root diagnostics")
    d["handoff_fixedpoint_root_robustness"]={
       "source_workflow_run":37763535074,
       "helper_root_pairs":handoff_scan,
       "selected_root_min_over_scan":min(roots),
       "selected_root_max_over_scan":max(roots),
       "root_spread":max(roots)-min(roots),
       "nearby_root_minus_helper_min":min(resids),
       "nearby_root_minus_helper_max":max(resids),
       "root_residual_sign_changes":sign_reversals,
       "other_curvature_derivative_roots_at_helper_16p7422":[29.61244254821877,11.30344270081078],
       "conclusion":"The observed curvature-derivative map has several event roots and numerical variation across nearby helper inputs. The 16.742... event is a reproducible candidate, but 8-figure unique physical selection is not established.",
       "clarification":"Other event roots at fixed helper are not proven additional fixed points; they are different zeros of the curvature diagnostic.",
    }
    tests.append("handoff root multiplicity and numerical sensitivity")


    # Late lobes: exactly invertible once the H/Hflat residual is *given*,
    # but those residual anchors encode observational or additional
    # action-dependent background information.
    A,B=p["A_late"],p["B_late"]
    z1,z2=.2,.7
    av,bv=solve_late_from_two_residuals(z1,z2,delta(z1,A,B),delta(z2,A,B))
    need(close(A,av,1e-13) and close(B,bv,1e-13),"two point lobe extraction")
    # delta'(0)=A and delta''(0)=-2A/tauA-2B
    h=1e-4
    deriv1=(delta(h,A,B)-delta(-h,A,B))/(2*h)
    deriv2=(delta(h,A,B)-2*delta(0,A,B)+delta(-h,A,B))/(h*h)
    need(close(deriv1,A,1e-7) and close(deriv2,-2*A/TAU_A-2*B,1e-7),
         "local lobe derivatives")
    d["late_lobes"]={"status":"conditionally extractable from H(z) or distances, not uniquely action-selected",
        "A_late":A,"B_late":B,"tauA_fixed_ansatz":TAU_A,"tauB_fixed_ansatz":TAU_B,
        "delta_prime_at_z0":A,
        "delta_second_derivative_at_z0":-2*A/TAU_A-2*B,
        "inversion_from_two_given_residual_points":{"z1":z1,"z2":z2,"A_recovered":av,"B_recovered":bv},
        "q0_relation":"A_late=q0-q0_control-0.5*f_prime0/(1-f0) if the H profile includes tracker f(z); reduces to q0-q0_control when f0~0 and f_prime0~0",
        "missing":"Independent action-calculated q0 and expansion curvature/distances (and microphysical tauA, tauB)"}
    tests.append("late-lobe unique algebraic inversion given data anchors")

    # Infocriteria: conditional HSC TATT results; no AIC/BIC family
    # dimension can be asserted solely from microscopic analogies.
    sensitivity=[]
    for dk in [0,1,2,3,6,10]:
        sensitivity.append({"delta_k_assumed":dk,
            "delta_AIC_HSC":HSC_DELTA_CHI2+2*dk,
            "delta_BIC_HSC":HSC_DELTA_CHI2+dk*math.log(HSC_N)})
    joint_d=PLANCK_DLOGZ+HSC_DLOGZ
    joint_err=math.hypot(PLANCK_ERR,HSC_ERR)
    d["model_selection"]={
        "status":"conditional HSC AIC/BIC bands and nuisance evidence; full-family evidence still open",
        "HSC_n":HSC_N,"HSC_TATT_delta_chi2":HSC_DELTA_CHI2,
        "sensitivity_by_assumed_extra_fitted_parameters":sensitivity,
        "combined_conditional_Planck_HSC":{"delta_logZ":joint_d,"error":joint_err,"central_Bayes_factor":math.exp(joint_d)},
        "warning":"Reconstruction-stage delta_k=0 is not a fair automatic model-family parameter count; data-selected ancestors and prior volumes matter.",
    }
    tests.append("conditional model-selection arithmetic")

    d["coordinate_ledger"]={
       "AF":"lower stability bound conditional on derived zc,width,zt; value not uniquely fixed",
       "z_c":"force-equality root conditional on reconstructed action/profile",
       "DeltaN_F":"force-profile descriptor conditional on reconstructed action/profile",
       "D_floor":"conditional kessence/tracker normalization; independent rK missing",
       "D0":"candidate structural normalization; not uniquely action-derived",
       "pD":"pD=1 leading analytic response; higher-power alternatives remain",
       "lambda_e":"inverts fr but fr originates from imposed lambda unless independently predicted",
       "z_t":"curvature-helper fixed point conditional on definition of map",
       "A_late":"determined by q0 if q0 and tracker derivative independent",
       "B_late":"determined by next background derivative or independent acoustic-distance condition",
    }
    d["theory_independence_verdict"]={
      "fully_derived_with_present_inputs":["Omega_x0 from flat closure"],
      "not_proven_as_independent_first_principles_numbers":list(d["coordinate_ledger"].keys()),
      "requirements_to_close":[
        "independent F(S) or Planck-mass response field equation with physical boundary data",
        "independent kinetic action functions and UV state selecting D_floor,D0,pD",
        "independently predicted primordial/early structural fraction selecting lambda_e",
        "unambiguous global activation/handoff solution selecting z_t",
        "independent H(z) or late expansion law selecting A_late and B_late",
      ],
      "no_fitted_coordinate_removed_for_family_evidence_by_this_audit":True,
    }
    d["tests_passed"]=tests
    d["test_count"]=len(tests)
    return d


if __name__=="__main__":
    output=audit()
    folder=Path("output/late300_all_coordinate_closure")
    folder.mkdir(parents=True,exist_ok=True)
    path=folder/"coordinate_ledger.json"
    path.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("ALL_COORDINATE_CLOSURE_TESTS_PASS",len(output["tests_passed"]))
    print("ALL_COORDINATE_CLOSURE_SUMMARY",json.dumps({
        "status":"algebraic derivations tested; none of the remaining ten is proven uniquely selected by these identities",
        "test_count":output["test_count"],
        "Omega_x0":output["Omega_x0"]["derived"],
        "fr":output["tracker_and_rK"]["fr_inferred"],
        "rK_conditional":output["tracker_and_rK"]["rK_conditional_inferred"],
        "delta_logZ_conditional":output["model_selection"]["combined_conditional_Planck_HSC"],
        "coordinates":output["coordinate_ledger"],
    },sort_keys=True))
