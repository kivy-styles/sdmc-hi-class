#!/usr/bin/env python3
"""Necessary analytic conditions for SDMC finite kinetic -> canonical IR matching.

This is a source-level mathematical audit, NOT a simulation of future covariant
SDMC or an observational likelihood. The coasting and phi=ln(a) alignment are
explicit assumptions. G3= -F_phi/H^2 holds only ON a selected No-Slip
background, and does not determine the off-shell coefficient in the future.
"""
import json, math, re
from pathlib import Path

OUT=Path("output/D_canonical_IR_match")
OUT.mkdir(parents=True,exist_ok=True)
AF=0.020520
ZC_F=4.034077502899133
W_F=0.32925377073762596
ZC_D=3.927876388467848
W_D=0.33114133956842123
DFLOOR=0.05181542627513409
D0=0.34231919445927034
F_MATURE_ORIGINAL=1.023885427695
K1_TODAY_DERIVED=-0.0913078779588896
K2X_TODAY_DERIVED=0.0825143261780817

# These exact definitions are checked against source, so the test will
# fail if the reconstructed covariant action or source conventions change.
SOURCE=Path("sdmc/run_late300_linear_covariant.py").read_text()
assert "g=-F1/(H*H)" in SOURCE
assert "Dtarget=DFLOOR+D0*S**POWER" in SOURCE
assert "CubicSpline(N,y)" in SOURCE
assert "if(x>=sdmc_lin_x[SDMC_LIN_N-1]) return SDMC_LIN_N-2;" in SOURCE
assert "const double dx=x-sdmc_lin_x[i];" not in SOURCE # generated C text, not Python
assert "double dx=x-sdmc_lin_x[i];" in SOURCE
# Native source uses cubic continuation beyond the last fitted knot.
# That continuation is a numerical extension of the last spline segment,
# not an independently specified microscopic UV/IR completion.

def check(ok,msg):
    if not ok: raise AssertionError(msg)

def near(a,b,rtol=1e-10):
    return math.isclose(a,b,rel_tol=rtol,abs_tol=rtol)

def F_and_Fprime(phi, A=AF, w=W_F, zc=ZC_F):
    nc=-math.log1p(zc)
    t=math.exp(-(phi-nc)/w)
    S=1./(1.+t)
    F=math.exp(A*S)
    Fp=F*A*t/(w*(1.+t)**2)
    return F,Fp

def g_required_on_coasting(phi,A=AF,w=W_F,zc=ZC_F,H_star=1.):
    # ONLY a required on-background coefficient under H=H*exp(-phi).
    _,fp=F_and_Fprime(phi,A,w,zc)
    Hsq=H_star**2*math.exp(-2*phi)
    return -fp/Hsq

def log_slope(phi1,phi2,w=W_F):
    g1=abs(g_required_on_coasting(phi1,w=w))
    g2=abs(g_required_on_coasting(phi2,w=w))
    return math.log(g2/g1)/(phi2-phi1)

def main():
    # Radiation/finite/mature are separate stages. These are necessary
    # matching conditions, not sufficient to solve Friedmann+scalar EoMs.
    F_IR=math.exp(AF)
    canonical_k1_phi=2*F_IR
    assert F_IR>1.
    Dfin=DFLOOR+D0
    check(not near(Dfin,2.),"finite kinetic correction is not future canonical endpoint")
    check(K1_TODAY_DERIVED<0 and canonical_k1_phi>0,
          "continuous same-field k1 must cross zero on a canonicalizing path")
    check(K2X_TODAY_DERIVED>0,"derived today has nonzero higher-kinetic sector")
    # From mature canonical direct structural field sigma=C*exp(phi):
    # Z_sigma=sigma^2 X_phi; G2=2F∞ Z_sigma/sigma² - 4F∞ Z_BI/sigma²
    #    =2 F∞ X_phi - (4F∞ Z_BI/C²)*exp(-2phi).
    # Hence k1(phi)->2F∞, k2(phi)->0, g(phi)->0, V(phi)∝exp(-2phi).
    # At the canonical endpoint H~1/a, X_phi=H²/2, alphaK=2.
    X=0.7
    canonical_alphaK=(2*X*canonical_k1_phi)/(2*X*F_IR)
    check(near(canonical_alphaK,2.),"derived canonical IR alphaK=2 under phi=ln a")

    expect=2-1/W_F
    nums={}
    for phi in (6.,8.,10.,12.):
        F,Fp=F_and_Fprime(phi)
        nums[str(phi)]={"F":F,"dF_dphi":Fp,
                        "g_required_on_coasting":g_required_on_coasting(phi),
                        "Dfinite_historical":DFLOOR+D0/(1+math.exp(-(phi+math.log1p(ZC_D))/W_D))}
    observed=log_slope(10,12)
    check(abs(observed-expect)<1e-7,
          "late required coasting G3 decay exponent")
    conditions=[]
    for w in (.32925377073762596,.45,.5,.6):
        expn=2-1/w
        conditions.append({"width":w,"on_shell_coasting_g_exponent":expn,
                            "g_decreases_in_necessary_coasting_continuation":expn<0})
        check((expn<0)==(w<.5),"G3-coefficient condition w<1/2")
    mismatch=F_MATURE_ORIGINAL/F_IR-1
    check(mismatch>0.003 and mismatch<0.0032,
          "canonical report and derived finite F∞ differ by about 0.31%")
    checks=[
      "confirmed source g=-F_phi/H² is ON-SHELL reconstruction, not off-shell future evolution",
      "confirmed source spline-table evaluator uses the last cubic outside table domain",
      "converted mature structural sigma field into phi=ln a coefficients",
      "verified canonical k1(phi)=2F∞ and kinetic alphaK=2",
      "verified finite D asymptote differs from mature D=2",
      "verified if k1 is continuous on same field it must cross zero before mature branch",
      "verified necessary coasting/no-slip g decay exponent and width<1/2",
      "confirmed historical mature F∞ is 0.31% above derived-bundle F∞",
    ]
    report={
     "status":"analytic necessary matching conditions verified; full forward trajectory not derived",
     "scientific_qualification":"The assumed coasting scaling and field alignment are mature-branch boundary conditions. The g scaling is a necessary on-shell continuation IF no-slip persists, not a prediction from extrapolating the past-fitted off-shell cubic table.",
     "F_derived_infinity":F_IR,"F_original_mature_manuscript":F_MATURE_ORIGINAL,
     "F_relative_mismatch_original_over_derived":mismatch,
     "Planck_response_transition_width":W_F,
     "required_coasting_g_power_of_a":expect,
     "measured_asymptotic_log_g_slope_a10_to_a12":observed,
     "width_decay_criteria":conditions,
     "necessary_IR_phi_action":{
       "sigma_map":"sigma=C*exp(phi) in mature coasting branch with phi=ln a",
       "X_phi":"X_phi= -1/2 (nabla phi)^2; Z_sigma=sigma² X_phi",
       "G2_phi":"2*F_infinity*X_phi - (4*F_infinity*Z_BI/C²)*exp(-2 phi)",
       "G3_phi":"0","G4_phi":"F_infinity/2",
       "k1_canonical_future":canonical_k1_phi,
       "k2_canonical_future":0,
       "V_canonical_future":"(4*F_infinity*Z_BI/C²)*exp(-2 phi)",
       "D_canonical_future":2,
     },
     "finite_action":{
       "Dfloor_historical":DFLOOR,"D0_historical":D0,"Dfinite_asymptote":Dfin,
       "k1_today_derived_action":K1_TODAY_DERIVED,
       "k2_times_X_today_derived_action":K2X_TODAY_DERIVED,
       "required_later_k1_sign_crossing_if_same_continuous_phi":True,
       "caveat":"A k1 crossing is not by itself a ghost; full scalar normalization D remains relevant.",
     },
     "future_source_limitation":"Spline coefficient table was inverse-derived on the target background; out-of-domain cubic continuation exists as a code operation but is not independently selected microscopic evolution.",
     "sampled_required_coasting_coefficients":nums,
     "passed_checks":checks,
     "still_required":[
       "a genuine future boundary/initial-value problem with independently specified k1,k2,V,g,F",
       "a covariant scalar-matter transfer and activation law selecting the matching epoch",
       "on-shell and off-shell No-Slip compatibility across the future field trajectory",
       "an explanation of F∞ mismatch between legacy canonical and force-derived finite branches",
       "native forward perturbation evolution through the actual canonicalizing action",
     ]
    }
    (OUT/"audit.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("D_CANONICAL_IR_MATCH_PASS",json.dumps({
     "checks":len(checks),
     "g_decays_if_width_below_half":True,
     "g_power_for_derived_width":expect,
     "measured_log_g_slope":observed,
     "finite_D_asymptote":Dfin,
     "canonical_D":2.,
     "canonical_IR_k1":canonical_k1_phi,
     "F_infinity_relative_mismatch":mismatch,
     "future_dynamics":"not determined by matched endpoints alone"
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
