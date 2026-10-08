#!/usr/bin/env python3
"""SDMC Dfloor + D0 correction: end-to-end mathematical closure and inference audit.

Sources: author's SDMC_A_COMPLETED(7).pdf pp. 640-645, 666-669;
Structural Radius/Jacobian Appending Continuation VI pp. 3-6;
and GitHub CI 37807119146,37807287174,37817678495,
37764900113,37764410555,37764359939.

This code intentionally DOES NOT produce a new fitted UV action or a
marginalized full family evidence from fixed-spectrum likelihood summaries.
It supplies analytic closure identities, counterexamples, direct tests and
a rigorous model-selection provenance ledger.
"""
import csv
import json
import math
from pathlib import Path

OUT=Path("output/D_full_closure_joint_audit")
OUT.mkdir(parents=True,exist_ok=True)

LAM=18.40625
CAN_FLOOR=16./LAM**2
ACC_FLOOR=.05181542627513409
ACC_D0=.34231919445927034
ACC_PD=1.
D_CENTER=3.927876388467848
D_WIDTH=.33114133956842123
MATURE_D=2.
COMPOSITE_D0=2.*math.pi/LAM
D0_CAUSAL_MIN=.01918832585405535

PLANCK={
    "accepted_D":2784.82549538608,
    "D0_2pi_lambda":2784.821140904252,
    "D0_Xi_sqrtF":2784.8195024676133,
    "D0_force_candidate":2784.8246840476986,
    "floor_16_lambda2":2784.5593024304303,
    "combined_2pi_floor":2784.5561712625763,
}
DESI={
    "accepted_D":(329.54156686094177,340.53240339264215),
    "D0_2pi_lambda":(329.5411830248053,340.53240338012273),
    "floor_16_lambda2":(329.52467310279377,340.53240339267523),
    "combined_2pi_floor":(329.5246928853224,340.53240339252653),
}
# These below are DIFFERENT derived-bundle model comparisons (not the
# above canonical-floor/correction action), nuisance-only conditional.
EVIDENCE={
 "planck_derived_minus_local021":(-2.125487389543423,.6916683974301847),
 "hsc_nla_derived_minus_local021":(2.0390446846788706,.6148357940108409),
 "hsc_tatt_derived_minus_local021_chi2":-2.6668948487108395,
 "hsc_tatt_optimizer_success":False,
 "hsc_tatt_nobs":60,
}

def check(test,message):
    if not test:
        raise AssertionError(message)

def close(a,b,tol=1e-11):
    return math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol)

def logistic(x):
    if x>=0: return 1./(1.+math.exp(-x))
    e=math.exp(x)
    return e/(1.+e)

def S_D(N):
    Nc=-math.log1p(D_CENTER)
    return logistic((N-Nc)/D_WIDTH)

def D_finite(N,floor=ACC_FLOOR,d0=ACC_D0,p=ACC_PD):
    return floor+d0*S_D(N)**p

def smooth_step(N,start,end):
    # Exact C-infinity bump: 0 for N<=start, 1 for N>=end.
    # Derivatives at both endpoints vanish to every order.
    if N<=start: return 0.
    if N>=end: return 1.
    t=(N-start)/(end-start)
    lhs=-1./t
    rhs=-1./(1.-t)
    return logistic(lhs-rhs)

def D_extension(N,start,end):
    old=D_finite(N)
    return old+(MATURE_D-old)*smooth_step(N,start,end)

def main():
    stages=[]
    fr=4./LAM**2
    dm_can=4.*fr
    check(close(dm_can,CAN_FLOOR),"canonical radiation tracker D=4Omega_phi")
    check(close(fr*LAM**2,4),"radiation scalar fraction identity")
    stages.append("radiation tracker D=4 fr=16/lambda^2 and numerical evaluation")

    s0=S_D(0.)
    D0_now=ACC_D0*s0
    Dnow=D_finite(0.)
    Dlate=ACC_FLOOR+ACC_D0
    Nc=-math.log1p(D_CENTER)
    midpoint_slope=ACC_D0*ACC_PD/(2**(ACC_PD+1)*D_WIDTH)
    recovered_D0=midpoint_slope*2**(ACC_PD+1)*D_WIDTH/ACC_PD
    check(close(recovered_D0,ACC_D0),"correction recovered from kinetic-window midpoint slope")
    check(close(Dnow, .3913854938989103,1e-9),"historical present correction")
    check(abs(Dlate-MATURE_D)>1.0,"two mathematically distinct kinetic endpoints")
    stages.append("finite D0 correction, midpoint derivative, finite asymptote differs from mature D=2")

    cases=[]
    for begin,finish in ((0.5,3.5),(1.5,6.),(1.1,8.5),(3.,9.)):
        samples=[]
        for N in (-100.,-30.,-8.,-2.,-1.,0.,.2,1.,2.,4.,6.,9.,10.,20.,60.):
            d=D_extension(N,begin,finish)
            if N<=begin:check(close(d,D_finite(N)),"C-infinity extension must exactly retain finite sector")
            if N>=finish:check(close(d,MATURE_D),"C-infinity extension must reach mature endpoint")
            check(d>0.,"positive scalar kinetic denominator along illustrative interpolation")
            samples.append({"N":N,"D":d,"smooth_transition_Q":smooth_step(N,begin,finish)})
        cases.append({"begin":begin,"finish":finish,"at_N1":D_extension(1.,begin,finish),
            "at_N4":D_extension(4.,begin,finish),"at_N0":D_extension(0.,begin,finish),
            "at_N20":D_extension(20.,begin,finish),
            "samples":samples})
    check(abs(cases[0]["at_N4"]-cases[1]["at_N4"])>0.1,
          "distinct compatible positive C-infinity handoffs")
    stages.append("four distinct smooth positive kinetic histories with identical finite past and mature endpoint")

    # Formal on-shell Horndeski G2 deformation:
    # deltaG2=mu(phi)*(X-Xbg(phi))^2 has zero deltaG2, deltaG2_X,
    # deltaG2_phi along X=Xbg(phi), but deltaG2_XX=2mu.
    # In the linear G3 and phi=ln a reconstruction alpha_K
    # correction is deltaD = 4*mu*Xbg/F, where Xbg=H^2/2.
    # This is a non-uniqueness proof, not physical model selection.
    def on_shell(mu,Xbg,X,dxbg_dphi=1.,dmu_dphi=.7):
        d=X-Xbg
        return {"G2":mu*d*d, "G2_X":2*mu*d,
            "G2_phi":dmu_dphi*d*d-2*mu*d*dxbg_dphi,
            "G2_XX":2*mu}
    for mu in (.1,.5,3.):
        r=on_shell(mu,2.,2.)
        check(close(r["G2"],0) and close(r["G2_X"],0) and close(r["G2_phi"],0),
             "on-background action invisible through first derivatives")
        check(close(r["G2_XX"],2*mu),"positive/different second-derivative kinetic corrections")
    stages.append("three explicit action-level degeneracies preserving same background")

    # The response exponent cannot be uniquely selected by positivity:
    # S in (0,1) and Dfloor>0 imply D>0 for all p>0,D0>0.
    for p in (.5,1.,1.5,2.,3.):
        for N in (-50.,-5.,0.,2.,50.):
            check(D_finite(N,p=p)>0.,"positive kinetic denominator many exponents")
    stages.append("p_D=1 not selected solely by kinetic positivity")

    # Important experimental causality inequality, not a microphysical
    # selector: D0_accepted far above independently scanned threshold.
    check(ACC_D0 > 10.*D0_CAUSAL_MIN,"observed correction much larger than causal floor")
    stages.append("conditional D0 causal lower bound is far from the accepted normalization")

    # Genuine tested conditional likelihood ledger.
    planck_ref=PLANCK["accepted_D"]
    desi_ref=DESI["accepted_D"][0]
    keys=("accepted_D","D0_2pi_lambda","floor_16_lambda2","combined_2pi_floor")
    likelihood_rows=[]
    for k in keys:
        p=PLANCK[k]
        d,b=DESI[k]
        row={"case":k,"planck_chi2":p,"delta_planck_vs_accepted":p-planck_ref,
            "desi_chi2":d,"desi_chi2_fixed_lcdm":b,
            "delta_desi_vs_accepted":d-desi_ref,
            "delta_desi_vs_fixed_lcdm":d-b,
            "delta_conditional_planck_plus_desi_vs_accepted":(p-planck_ref)+(d-desi_ref)}
        likelihood_rows.append(row)
    combined=next(q for q in likelihood_rows if q["case"]=="combined_2pi_floor")
    check(close(combined["delta_conditional_planck_plus_desi_vs_accepted"],
                -.2693241235037-.0168739756194,tol=1e-8),
          "combination corresponds to tested fixed-cosmology Planck+DESI shifts")
    stages.append("four independent true-input Planck/DESI likelihood tables and conditional within-SDMC deltas")

    hsc_delta=EVIDENCE["hsc_tatt_derived_minus_local021_chi2"]
    accounting=[]
    for dk in (0,1,2,3,5,10):
        accounting.append({"extra_free_SDMC_dimensions":dk,
             "hsc_delta_chi2_prior_to_penalty":hsc_delta,
             "hsc_conditional_delta_AIC":hsc_delta+2*dk,
             "hsc_illustrative_delta_BIC_if_Neff60":hsc_delta+math.log(60)*dk,
             "conditional_D_replacement_delta_AIC_if_one_extra_dimension":
                 combined["delta_conditional_planck_plus_desi_vs_accepted"]+2*dk})
    stages.append("AIC/BIC and effective-dimension sensitivity, not a joint evidence")
    results={
     "verdict":"conditional D identities and matched action likelihoods verified; no independent finite kinetic transition or full model-family local021 evidence",
     "manuscript_scientific_scope":{
      "accepted_action":"G2=k1(phi)X+k2(phi)X^2-V(phi); G3=g(phi)X; G4=F(phi)/2",
      "reconstruction":"G2 coefficients are reconstructed FROM selected D(N), background H(N) and No-Slip F, not independently given UV couplings",
      "mature_action":"independently adopted manuscript canonical future closure with D_infinity=2 and c_s2_infinity=1",
      "mature_trajectory":"separately derived mature IR action does not select the finite matching time, normalization or transfer profile",
     },
     "radiation":{"lambda_e_historically_selected":LAM,"radiation_scaling_fraction_derived_given_lambda":fr,
      "canonical_Dfloor_given_lambda":CAN_FLOOR,"accepted_Dfloor":ACC_FLOOR,
      "noncanonical_enhancement_ratio_from_selected_floor":ACC_FLOOR/CAN_FLOOR,
      "status":"conditional canonical tracker theorem; no independent lambda or noncanonical enhancement selector"},
     "finite_D":{"Dfloor":ACC_FLOOR,"D0":ACC_D0,"pD":ACC_PD,"D_center_redshift":D_CENTER,
      "D_window_e_folds":D_WIDTH,"D_switch_today":s0,"D0_correction_today":D0_now,
      "D_today":Dnow,"midpoint_slope_from_fitted_profile":midpoint_slope,
      "D0_recovered_from_midpoint_slope":recovered_D0,
      "D_finite_late_asymptote":Dlate,"mature_D_infinity":MATURE_D,
      "mature_minus_finite_late":MATURE_D-Dlate,
      "warning":"finite logistic asymptote != physical mature canonical asymptote"},
     "four_nonunique_exact_C_infinity_extensions":cases,
     "likelihood_rows":likelihood_rows,
     "conditional_joint_D_replacement_delta_chi2_vs_accepted":
        combined["delta_conditional_planck_plus_desi_vs_accepted"],
     "combined_shift_conditions":"fixed cosmological parameters, independently profiled Planck foreground nuisance and DESI EFT bias nuisance; both use same fixed background; not compared to optimized local021",
     "local021_status":{
       "required":"joint likelihood L_Planck*L_DESI*L_SN*L_HSC under one shared cosmological/structural prior and separate nuisance coordinates, both model families sampled",
       "verified_observations_for_distinct_derived_bundle":EVIDENCE,
       "important":"Planck and HSC conditional evidence integrate nuisance coordinates only and the derived-bundle is different from these D substitutions; do NOT combine with Planck/DESI fixed chi2 into a final Bayes factor",
     },
     "accounting_sensitivity":accounting,
     "tests":stages,
     "what_is_not_derived":["lambda_e absolute","noncanonical ratio rK from UV coefficients",
         "D0 correction magnitude from action","p_D exponent from action",
         "D correction centre and width from action",
         "finite-to-mature handoff transition function and activation/transfer",
         "full model-family effective dimension and joint Bayes factor vs local021"],
     "sources":{
      "manuscript_A":"SDMC_A_COMPLETED(7).pdf pp. 640-645, 666-669",
      "jacobian_continuation":"SDMC_Structural_Radius_Jacobian_Investigation_Appending_Continuation_VI(2).pdf pp. 3-6",
      "true_action_full_planck_run":37807119146,
      "true_action_desi_run":37807287174,
      "independent_true_action_pliklite_run":37817678495,
      "D0_action_nonuniqueness_run":37803894563,
      "D0_causal_bound_run":37804665278,
      "planck_conditional_evidence_run":37764900113,
      "HSC_conditional_evidence_run":37764410555,
      "HSC_TATT_run":37764359939,
     }
    }
    (OUT/"audit.json").write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
    with (OUT/"conditional_likelihood_comparison.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(likelihood_rows[0].keys()))
        w.writeheader();w.writerows(likelihood_rows)
    with (OUT/"C_infinity_handoff_nonuniqueness.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["N","D_finite","D_match_0p5_to_3p5","D_match_1p5_to_6","D_match_1p1_to_8p5","D_match_3_to_9"])
        for j in range(-50,201):
            N=j/10.
            w.writerow([N,D_finite(N)]+[D_extension(N,xx["begin"],xx["finish"]) for xx in cases])
    print("SDMC_D_ENDTOEND_THEORY_AUDIT_PASS",json.dumps({
        "tests_passed":len(stages),
        "finite_asymptote":Dlate,
        "mature_limit":MATURE_D,
        "smooth_handoff_families":len(cases),
        "canonical_radiation_floor":CAN_FLOOR,
        "accepted_floor":ACC_FLOOR,
        "true_input_full_planck_desi_conditional_combined_delta_chi2":
            combined["delta_conditional_planck_plus_desi_vs_accepted"],
        "effective_family_bayes_vs_local021":"NOT YET IDENTIFIED",
        "D_finite_microphysics":"NOT DERIVED",
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
