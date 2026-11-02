"""SDMC near-flatness suppression audit (Manuscript B Part I §§28-33, 53-55).

This is a *kinematical counterfactual* built from the ARCHIVED FLAT
late300 background numbers and an independent hypothetical spatial
curvature, not a numerically integrated non-flat late300 FLRW solution.

Manuscript's diagnostic:
 Omega_kR,rad / Omega_kR,P = 1/Nrad² =0.0858219;
 Omega_kJ,0 / Omega_kR,P = 0.0917953 * Ck0;
 Omega_kJ,freeze / Omega_kR,P = 6.35e-17 * Ckfreeze;
 Omega_kJ,zstar / Omega_kR,P =2.10696e-4 * Ckstar.

Ck=kJ/kR. If kJ and kR are independently CONSTANT curvature labels
on exact FLRW spatial slices, Ck is necessarily constant in time.
If gamma_ij=M(t)² h_ij and a_s=M(t)*a_J with a common spatial
coordinate/curvature normalization, intrinsic scalar curvature gives
  6 kR/a_s²=M^-2 * 6 kJ/a_J², so kR=kJ and Ck=1.
A genuinely evolving Ck thus requires leaving at least one of those
assumptions (non-FLRW transition geometry, time-dependent effective
curvature in coarse-grained slices, or a nontrivial spatial mapping).
It CANNOT be obtained by allowing an arbitrary scalar M(t) alone.

For Ck=1 and target |Omega_k,J,0|<1e-3,
  |Omega_k,R,P|< 0.0108938039.
For |Omega_k,R,P|=1,
  Ck0 <0.0108938039.

The 435.6765 growth factor between z*=1090 and today is the ratio
of manuscript's conditional Omega values and is fixed entirely by
the reference a_J H_J history for CONSTANT kJ.

Physical cautions:
 - Omega at the exact prethermal Planck boundary is a provisional
   *structural* curvature diagnostic, not classical Jordan FLRW.
 - Early mapped small Omega does not require a NEW curvature-damping
   mechanism when kJ is fixed; a² H² determines it kinematically.
 - This does not calculate the dynamically realized kJ, nor show
   that any primordial distribution produces small present curvature.
 - For order-unity counterfactual curvature the archived flat
   background is not a self-consistent curved solution.
"""
from dataclasses import dataclass
import math
import json

N_PLANCK=1.
N_CROSS=2/3
N_RAD=3.41350846
CURVATURE_PRESENT_FACTOR=.0917953
CURVATURE_RECOMB_FACTOR=.000210696
CURVATURE_FREEZE_FACTOR=6.35e-17
RECOMB_Z=1090.

def structural_history():
    p1=(N_PLANCK/N_CROSS)**2
    p2=(N_CROSS/N_RAD)**2
    total=p1*p2
    return dict(phase1=p1,phase2=p2,total=total,
                percent_reduction=100*(1-total),
                final_from_archived_table=.0858219,
                archive_residual=total-.0858219)

def conditional_curvature(initial=1.,C0=1.,Cstar=1.,Cfreeze=1.):
    if initial<0 or min(C0,Cstar,Cfreeze)<0:
        raise ValueError("nonnegative curvature magnitudes and ratios")
    return dict(initial_structural_fraction=float(initial),
      present=float(initial*C0*CURVATURE_PRESENT_FACTOR),
      recombination=float(initial*Cstar*CURVATURE_RECOMB_FACTOR),
      freezeout=float(initial*Cfreeze*CURVATURE_FREEZE_FACTOR),
      C0=float(C0),Cstar=float(Cstar),Cfreeze=float(Cfreeze))

def required_extra_at_present(initial=1.,target=.001):
    if initial<=0 or target<=0:
        raise ValueError("positive initial curvature magnitude and target")
    max_C=target/(CURVATURE_PRESENT_FACTOR*initial)
    return dict(initial=float(initial),
       threshold=float(target),
       maximum_C0=float(max_C),
       suppression_needed_from_C1=float(max(1.,1/max_C)),
       adequate_with_C1=bool(max_C>1.))

def archival_growth_factors():
    growth=CURVATURE_PRESENT_FACTOR/CURVATURE_RECOMB_FACTOR
    return dict(recombination_to_present_growth=growth,
                freezeout_to_present_growth=
                    CURVATURE_PRESENT_FACTOR/CURVATURE_FREEZE_FACTOR,
                inferred_Hstar_over_H0=
                    (1+RECOMB_Z)*math.sqrt(growth))

def curvature_mapping_check(aJ,M,HJ,p,kJ):
    """Identity via explicitly mapped FLRW kinematics.

    Evaluate physical dimensionless curvature fractions in common
    curvature normalization. HR=p HJ, aS=M aJ.
    """
    if min(aJ,M,HJ,p)<=0:
        raise ValueError("positive scales and expansion rates")
    aS=M*aJ
    HR=p*HJ
    kR=kJ  # pure time-only conformal locking on same spatial slice
    omJ=-kJ/(aJ*aJ*HJ*HJ)
    omR=-kR/(aS*aS*HR*HR)
    ratio=omJ/omR if kJ!=0 else math.nan
    return dict(omegaJ=omJ,omegaR=omR,
       omegaJ_from_bridge=p*p*M*M*omR,
       algebraic_ratio=ratio,
       bridge_ratio=p*p*M*M,
       effective_Ck=1.)

def nonflat_FLRW_consistency(kR_initial,kR_final,kJ_initial,kJ_final):
    """Finite difference test of the fixed-FLRW-k condition.

    A nonzero evolving Ck is impossible if both kR and kJ are
    separately constant. This does NOT assume global spatial flatness.
    """
    fixedR=math.isclose(kR_initial,kR_final,rel_tol=1e-13,abs_tol=0.)
    fixedJ=math.isclose(kJ_initial,kJ_final,rel_tol=1e-13,abs_tol=0.)
    if kR_initial==0 or kR_final==0:
        raise ValueError("Ck nonzero-curvature diagnostic requires kR nonzero")
    Cfirst=kJ_initial/kR_initial
    Clast=kJ_final/kR_final
    return dict(fixed_structural_k=fixedR,fixed_jordan_k=fixedJ,
      Cfirst=Cfirst,Clast=Clast,
      C_changed=not math.isclose(Cfirst,Clast,rel_tol=1e-13,abs_tol=0.),
      both_FLRW_fixed_and_C_changed=bool(fixedR and fixedJ and
         not math.isclose(Cfirst,Clast,rel_tol=1e-13,abs_tol=0.)))

def diagnostic():
    return dict(scope="counterfactual kinematic curvature audit, NOT a solved nonflat SDMC background",
      history=structural_history(),
      conditions=[required_extra_at_present(x) for x in [1,.3,.1,.03,.01]],
      growth=archival_growth_factors(),
      snapshots=[conditional_curvature(x) for x in [1,.1,.01]],
      conformal_mapping=curvature_mapping_check(
        aJ=.001,M=.0001,HJ=25000.,p=2.,kJ=.0001),
      no_C_rise_under_two_fixed_curvatures=nonflat_FLRW_consistency(
        1.,1.,.5,.5),
      limitations=[
       "curvature perturbation feedback on late300 H_J not integrated",
       "physical C_k handoff map not derived",
       "global Planck boundary is not already classical Jordan FLRW",
       "observational target .001 is illustrative, not an absolute proof threshold"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
