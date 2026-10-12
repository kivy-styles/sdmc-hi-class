#!/usr/bin/env python3
"""Source-constrained SDMC/Late300 analytic transition and No-Slip audit.

Mathematical tests of an illustrative early/late two-exponential potential,
not a new autonomous Horndeski solution or parameter-free prediction.

References: sdmc/late300_candidate.json; sdmc/run_late300_linear_covariant.py;
sdmc/reports/activation_handoff_field_equations_2026-10-10.md.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np

PARAM_PATH=Path(__file__).resolve().parent/'late300_candidate.json'

def audit(p):
    AF=p['A_F']; zc=p['z_c']; w=p['width']; Dfloor=p['D_floor']; D0=p['D0']
    lam=p['lambda_e']; zt=p['z_t']; power=p['kinetic_power']
    assert power==1, 'Direct affine D-F relationship assumes kinetic_power=1'
    lir=math.sqrt(2)
    nc=-math.log1p(zc); nt=-math.log1p(zt)
    Delta=lam-lir
    def U_of_N(n):
        return 0.5*(1.+np.tanh((n-nc)/(2*w)))
    def U_of_z(z):
        return U_of_N(-np.log1p(z))
    z_grid=np.expm1(np.linspace(0.,math.log1p(1e7),1601))
    u=U_of_z(z_grid)
    F=np.exp(AF*u)
    D=Dfloor+D0*u
    alphaM=AF/w*u*(1-u)
    alphaB=-2*alphaM
    alphaK=D-1.5*alphaB**2
    # Necessary No-Slip kinetic composite from full restricted Horndeski alpha_K:
    # Kcov=2X(k1+6k2X-4g_phi X)/(M^2 F H^2)
    # D=Kcov-6*alphaM+6*alphaM^2. This is a target, not independent source.
    Kcov=D+6*alphaM-6*alphaM**2
    C=D0/AF
    ut=float(U_of_z(zt))
    u_today=float(U_of_z(0.))
    inverseD=Dfloor+C*np.log(F)
    nmid=nc
    canonical_field_rate=3/lam
    required_field_rate=1/(Delta*w)
    canonical_width=1/(Delta*canonical_field_rate)
    Dmid=Dfloor+D0/2
    braidmid=6*(AF/(4*w))**2
    # Matter-era canonical scaling fraction versus Late300's prescribed
    # early-tracker suppression at the midpoint z_t (its window is 1/2).
    matter_fraction_early=3/lam**2
    lambda_eff_t=lam-Delta*ut
    matter_fraction_slope=3/lambda_eff_t**2
    target_f_tr_at_zt=.5*matter_fraction_early
    # Include radiation's small correction in the actual published tracker formula.
    tracker_exact_zt=None
    if all(k in p for k in ('H0','omega_b','omega_cdm')):
        h=p['H0']/100
        matter=(p['omega_b']+p['omega_cdm'])/h**2*(1+zt)**3
        radiation=(4.17998772e-5/h**2)*(1+zt)**4
        wbg=(radiation/(matter+radiation))/3
        tracker_exact_zt=.5*3*(1+wbg)/lam**2
    tests_z=[1090.,100.,30.,zt,10.,zc,2.,1.,0.]
    samples=[]
    for zz in tests_z:
        ui=float(U_of_z(zz))
        Fi=math.exp(AF*ui); ai=AF/w*ui*(1-ui); Di=Dfloor+D0*ui
        samples.append({'redshift':zz,'U':ui,'F':Fi,'alphaM':ai,'D':Di,
                        'alphaK':Di-6*ai**2,
                        'Kcov_required':Di+6*ai-6*ai**2})
    # Approximate mismatch for an equal-φ kinetic velocity perturbation in the
    # accepted LINEAR G3 on-shell tuned model; not a nonlinear simulation.
    velocity_perturbations={
        str(f): {'X_over_Xbg':f*f,
                 'NoSlip_residual_over_alphaM':1-f*f,
                 'NoSlip_residual_at_center':AF/(4*w)*(1-f*f)}
        for f in (0.9,1.0,1.1)
    }
    # Hypothetical smooth logarithm is only an illustrative off-shell G3 candidate;
    # at X->0 it restores analyticity by allowing a No-Slip residual.
    regularized={
        str(r): {'Xstar_over_X':r,
                 'NoSlip_residual_over_alphaM':r/(1+r),
                 'NoSlip_residual_at_center':AF/(4*w)*r/(1+r)}
        for r in (0.001,.01,.1,1.)
    }
    out={
      'scope':'Analytic necessary conditions and candidate tests, NOT a freely evolved unified SDMC action',
      'source_parameters':{k:p[k] for k in ['A_F','z_c','width','D_floor','D0','kinetic_power','lambda_e','z_t']},
      'exact_late300_transition':{
        'C_D_over_lnF':C,
        'max_abs_identity_error':float(np.max(np.abs(D-inverseD))),
        'max_abs_derivative_identity_error':float(np.max(np.abs((D0/w)*u*(1-u)-C*alphaM))),
        'z_early_tracker_handoff':zt,'z_gravitational_center':zc,
        'N_early_tracker_handoff':nt,'N_gravitational_center':nc,
        'e_fold_separation':nc-nt,
        'U_at_early_handoff':ut,'U_at_F_midpoint':.5,'U_today':u_today,
        'steep_potential_fraction_at_early_handoff_if_U_is_potential_weight':1-ut,
        'F_at_early_handoff':float(math.exp(AF*ut)),
        'D_at_early_handoff':Dfloor+D0*ut,
        'D_at_transition_midpoint':Dmid,
        'D_finite_formula_asymptote':Dfloor+D0,
        'minimum_alphaK_on_grid':float(np.min(alphaK)),
        'minimum_Kcov_on_grid':float(np.min(Kcov)),
        'NoSlip_kinetic_composite_identity':'Kcov=2X(k1+6k2X-4g_phi X)/(M^2 F H^2)=D+6alphaM-6alphaM^2',
        'Kcov_at_midpoint':Dmid+6*AF/(4*w)-6*(AF/(4*w))**2,
        'samples_by_redshift':samples,
        'max_alphaM_on_grid':float(np.max(alphaM)),
        'direct_braid_squared_at_midpoint':braidmid,
        'braid_squared_over_midpoint_D_correction':braidmid/(D0/2)
      },
      'two_exponential_candidate':{
        'lambda_early':lam,'lambda_IR':lir,'slope_difference':Delta,
        'matter_tracker_canonical_phi_N':canonical_field_rate,
        'predicted_width_using_exact_matter_tracker_speed':canonical_width,
        'accepted_width':w,'width_fractional_mismatch':canonical_width/w-1,
        'field_speed_needed_at_transition':required_field_rate,
        'field_speed_fractional_excess_over_canonical_matter_tracker':required_field_rate/canonical_field_rate-1,
        'required_phi_N_per_structural_p_if_p_matter_1p5':required_field_rate/1.5,
        'effective_potential_slope_at_early_handoff_if_same_U':lambda_eff_t,
        'canonical_matter_scaling_fraction_from_steep':matter_fraction_early,
        'canonical_matter_scaling_fraction_from_effective_slope_at_zt':matter_fraction_slope,
        'canonical_matter_fraction_relative_change':matter_fraction_slope/matter_fraction_early-1,
        'Late300_matter_limit_tracker_fraction_at_zt':target_f_tr_at_zt,
        'Late300_actual_tracker_fraction_at_zt_with_radiation':tracker_exact_zt,
        'slope_needed_for_50pct_fraction_if_canonical_matter_tracker':math.sqrt(2)*lam,
        'Gamma_minus_1_at_zt_for_two_exponential':Delta**2*ut*(1-ut)/lambda_eff_t**2,
        'future_center_N_shift_per_log_amplitude_ratio_if_constant_speed':w,
        'missing_theory':'absolute potential amplitude ratio, scalar field normalization and dynamical speed are not fixed by slopes; two exponential model does not itself satisfy exact finite No-Slip'
      },
      'NoSlip_offshell_test':{
        'identity':'X G3_X=-(M^2/2) F_phi for G4=M^2 F/2, G5=0',
        'linear_G3_target_reconstruction':'g(phi)=-M^2 F_phi/(2 X_bg(phi)); NoSlip residual/alphaM=1-X/Xbg',
        'analytic_no_go':'If F_phi !=0 and G3_X stays finite at X=0, exact off-shell NoSlip at every X>=0 is impossible because X G3_X -> 0',
        'logarithmic_form_for_X_positive':'G3=-(M^2/2) F_phi ln(X/Xstar) + h(phi); exact off-shell NoSlip but nonregular at X=0',
        'regular_log_model':'G3=-(M^2/2) F_phi ln(1+X/Xstar) + h(phi); residual/alphaM=Xstar/(X+Xstar)',
        'same_phi_velocity_deviations':velocity_perturbations,
        'regularized_ratio_test':regularized,
        'warnings':['Off-shell logarithmic model is NOT the accepted linear-G3 action.',
                    'A small No-Slip residual is not a full scalar stability or CMB lensing likelihood.',
                    'No independent microscopic Xstar or C coefficient predicted.']
      },
      'scientific_verdict':{
        'connected_at_target_level':['same_U_for_F_and_D','affine_D_logF','same scalar transition may underlie selected profiles'],
        'not_independently_derived':['D0/AF ratio','A_F','z_c','width','z_t','finite k1,k2,g,V','field trajectory phi(N)','UV initial conditions'],
        'candidate_naive_equal_potential_epoch_at_zt':'DISFAVORED: at zt the steep term still contributes 97.75% if both share the accepted U.',
        'adiabatic_tracker_shutoff_by_potential_flattening':'FAIL: adiabatic canonical scaling fraction rises by ~4.28% whereas prescribed early component halves at zt.',
        'action_level_required_kinetic_composite':'DERIVED AS NECESSARY ON-SHELL TARGET, not independently sourced from microscopic coefficients.',
        'fundamental_delta_k':'UNDETERMINED',
        'new_full_observational_likelihood':'NOT PERFORMED'
      }
    }
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',default='output/late300_unified_transition_audit.json')
    opts=ap.parse_args()
    params=json.loads(PARAM_PATH.read_text())['parameters']
    data=audit(params)
    dst=Path(opts.output)
    dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('UNIFIED_TRANSITION',json.dumps({
      'U_at_zt':data['exact_late300_transition']['U_at_early_handoff'],
      'C':data['exact_late300_transition']['C_D_over_lnF'],
      'width_prediction':data['two_exponential_candidate']['predicted_width_using_exact_matter_tracker_speed'],
      'width_fractional_mismatch':data['two_exponential_candidate']['width_fractional_mismatch'],
      'identity_error':data['exact_late300_transition']['max_abs_identity_error'],
      'fundamental_delta_k':data['scientific_verdict']['fundamental_delta_k']},sort_keys=True))
if __name__=='__main__':main()