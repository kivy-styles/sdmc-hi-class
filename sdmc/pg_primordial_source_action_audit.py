"""SDMC Manuscript B Part III: physical primordial-well source-channel closure.

ACCEPTED action audit in Jordan variables:
  S_total=S_J[g_J,sigma]+S_b[g_J,psi_b]+S_r[g_J,psi_r]
          +S_support[g_J,sigma,L,M].

S_b and S_r are minimally coupled to the physical Jordan metric, and
S_J (F/G2/G3) contains no explicit independent L. Consequently at
fixed physical g_J, sigma and material fields:
   delta S_b/delta L=delta S_r/delta L=delta S_J/delta L=0.

For the accepted supporting-space scalar trace mode (Manuscript B,
Part III Eq 106.3), the ONLY displayed linear link to u=m-Phi_J is
through curvature-selector coefficient W:
  E(u,l) = 1/2 Qmetric u² +
           (K/4) omega0² l² +
           (C/2)(u+l/3)² - j_metric u,
  K=B f_L²/3>0,
  omega0²=m_L²+c_L² k²/a_J²>0,
  C=4W eta_K k^4/a_s^4>=0.

Qmetric is an externally specified, POSITIVE Euclidean diagnostic
stiffness, NOT the constrained GR lapse/shift/metric Hessian. j_metric
denotes an arbitrary metric-source drive, NOT a CDM density. It is used
ONLY to check mixing and Schur algebra, never to predict mu or Sigma.

The field Hessian in (u,l) is
  [[Qmetric+C, C/3],
   [C/3, K omega0²/2+C/9]].
After static elimination of l, its metric stiffness is
  Qmetric+K_geom,
  K_geom=C omega0²/(omega0²+Lambda),
  Lambda=2C/(9K).

At W=0 (or k=0) this EXACTLY factorizes into independent
u and l blocks. A baryonic metric source cannot drive a
stationary independent l within this action at linear order.

The total constrained Hessian of the physical rolling Horndeski
gravity is NOT the positive Qmetric toy. Nothing here claims that
positive curvature stiffness creates attractive extra gravity.

Separate route that survives W=0:
  scalar's accepted nonminimal term F(sigma) R_J/2
  has delta S/delta sigma ⊃ F_sigma R_J/2,
so Jordan curvature may drive the scalar even when L has no
direct source. That is not a new independently sourced L well
and is subject to the accepted late300 No-Slip constraints.

A nonzero PRE-EXCITED L can mix quadratically with sigma through
coefficient K(sigma), since mixed Hessian ∂²E/(∂sigma∂L)
~K_sigma omega0² L/2; this is zero in the exactly locked
L=0 branch. L energy/stress ~L² may gravitate without a direct
linear matter->L force. It requires an independently normalized
physical excitation population, not a static primordial-well
memory theorem.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np

@dataclass(frozen=True)
class SourceAudit:
    K:float=1.
    omega0_sq:float=1.8
    w_eta:float=.04
    k:float=1.5
    a_s:float=1.
    Qmetric:float=2.
    j_metric:float=.08
    K_sigma:float=.12
    F_sigma:float=.2

    def __post_init__(self):
        arr=np.asarray(list(vars(self).values()),dtype=float)
        if not np.all(np.isfinite(arr)):raise ValueError("finite arguments")
        if min(self.K,self.omega0_sq,self.a_s,self.Qmetric)<=0:
            raise ValueError("positive kinetic, restoring and metric diagnostics")
        if self.w_eta<0 or self.k<0:
            raise ValueError("positive curvature selector sector")

def curvature_stiffness(e:SourceAudit):
    return 4*e.w_eta*e.k**4/e.a_s**4

def lambda_trace(e:SourceAudit):
    return 2*curvature_stiffness(e)/(9*e.K)

def primordial_support_energy(e,u,ell):
    """Static reduced trace+metric diagnostic energy. No added L source."""
    C=curvature_stiffness(e)
    return (.5*e.Qmetric*u*u
       +e.K*e.omega0_sq*ell*ell/4
       +.5*C*(u+ell/3)**2-e.j_metric*u)

def bare_jordan_matter_action(e,u,ell):
    """Minimal matter contributes ONLY via physical Jordan metric.

    'u' here is a formal metric source channel; no assumption that the
    physical matter ADM action is literally -j*u is made.
    """
    return -e.j_metric*u

def analytic_action_gradient(e,u,ell):
    C=curvature_stiffness(e)
    return np.array([e.Qmetric*u+C*(u+ell/3)-e.j_metric,
         e.K*e.omega0_sq*ell/2+C*(u+ell/3)/3])

def analytic_hessian(e):
    C=curvature_stiffness(e)
    return np.array([[e.Qmetric+C,C/3],
        [C/3,e.K*e.omega0_sq/2+C/9]],dtype=float)

def finite_difference_hessian(e,step=1e-4):
    if step<=0:raise ValueError("positive step")
    out=np.zeros((2,2))
    q=np.zeros(2)
    for i in range(2):
        for j in range(2):
            if i==j:
                a=q.copy();b=q.copy();a[i]=step;b[i]=-step
                out[i,i]=(primordial_support_energy(e,*a)-
                      2*primordial_support_energy(e,*q)+
                      primordial_support_energy(e,*b))/step**2
            else:
                a=q.copy();b=q.copy();c=q.copy();d=q.copy()
                a[i]+=step;a[j]+=step;b[i]+=step;b[j]-=step
                c[i]-=step;c[j]+=step;d[i]-=step;d[j]-=step
                out[i,j]=(primordial_support_energy(e,*a)-
                   primordial_support_energy(e,*b)-
                   primordial_support_energy(e,*c)+
                   primordial_support_energy(e,*d))/(4*step*step)
    return out

def stationary_solution(e):
    H=analytic_hessian(e)
    rhs=np.array([e.j_metric,0.])
    u,ell=np.linalg.solve(H,rhs)
    C=curvature_stiffness(e)
    lam=lambda_trace(e)
    Kgeom=C*e.omega0_sq/(e.omega0_sq+lam)
    u_schur=e.j_metric/(e.Qmetric+Kgeom)
    ell_schur=-3*lam*u_schur/(e.omega0_sq+lam)
    # This Euclidean effective source is not a physical gravitational mu.
    return dict(u=float(u),ell=float(ell),baseline_u=e.j_metric/e.Qmetric,
        offdiag=float(H[0,1]),stiffness=Kgeom,
        extra_metric_stiffness=Kgeom,
        schur_u=float(u_schur),schur_ell=float(ell_schur),
        hessian_det=float(np.linalg.det(H)),
        physical_mu_or_lensing_predicted=False)

def metric_source_grad_ell(e,step=1e-5):
    return ((bare_jordan_matter_action(e,.1,step)-
             bare_jordan_matter_action(e,.1,-step))/(2*step))

def curvature_loads_scalar(F0,F_sigma,R,dsigma=0.):
    """Off-shell F(sigma)R/2 scalar source; F=F0+F_sigma sigma."""
    if F0<=0 or not all(np.isfinite([F0,F_sigma,R,dsigma])):
        raise ValueError("finite curvature loading and F>0")
    return dict(sigma_equation_curvature_term=.5*F_sigma*R,
      action_density=.5*(F0+F_sigma*dsigma)*R,
      depends_on_structural_L=False,
      nonzero_if_curvature_and_Fsigma=bool(F_sigma*R!=0))

def scalar_kinetic_coefficient_mixing(e,ell):
    """K(sigma)*omega² L²/4 yields mixed Hessian K_sigma omega² L/2."""
    return .5*e.K_sigma*e.omega0_sq*ell

def independent_well_energy_at_W_zero(K,mass2,ell,ell_dot):
    """Scalar trace energy can gravitate at O(L²) even if linear mixing zero."""
    if K<=0 or mass2<=0:
        raise ValueError("healthy quadratic branch")
    return K*(ell_dot**2+mass2*ell**2)/4

def source_free_static_lump_identity(field,K,mass2,c2,a_j,length):
    """Positive integrated elliptic identity forbids a regular static lump.

    For W=0 and no L-linear source, the scalar-trace static field eq is
      K [m² ell -(c²/a_J²)Delta ell]=0,
    so on a periodic domain or a decaying, no-flux boundary,
      <ell * EL> =
         K[m² <ell²>+(c²/a_J²)<|grad ell|²>] >=0.
    The on-shell condition EL=0 forces ell=0 for m²>0 and K>0.
    A constant massless zero-mode is a separate degenerate branch.

    This is a positive-energy field equation on a fixed spatial metric;
    it does not prevent ordinary Jordan gravitational curvature.
    """
    from sdmc.pg_full_spatial_stress import _cube,_diff
    field=_cube(field)
    if K<=0 or mass2<=0 or c2<0 or a_j<=0 or length<=0:
        raise ValueError("positive massive healthy branch required")
    grad=[_diff(field,i,length) for i in range(3)]
    lap=sum(_diff(grad[i],i,length) for i in range(3))
    EL=K*(mass2*field-(c2/a_j**2)*lap)
    positive=K*(mass2*float(np.mean(field*field))
       +(c2/a_j**2)*sum(float(np.mean(g*g)) for g in grad))
    integrated=float(np.mean(field*EL))
    return dict(integrated_field_times_EL=integrated,
      positive_mass_gradient_norm=positive,
      identity_residual=integrated-positive,
      supports_nonzero_sourcefree_static_lump=False,
      Jordan_gravitational_well_forbidden=False)

def candidate_explicit_direct_source(e,u,ell,g_new):
    """ONLY illustrates an extra operator -g_new*j_metric*ell.

    g_new is NOT in the displayed accepted SDMC action. It is kept
    separate to prevent silently claiming new direct matter coupling.
    """
    return primordial_support_energy(e,u,ell)-g_new*e.j_metric*ell

def candidate_source_solution_at_W_zero(e,g_new):
    """Formal conditional operator demonstration, NOT SDMC prediction."""
    if curvature_stiffness(e)!=0:
        raise ValueError("this diagnostic requires selector off")
    return 2*g_new*e.j_metric/(e.K*e.omega0_sq)

def on_shell_charged_source_audit(e):
    C=curvature_stiffness(e)
    stat=stationary_solution(e)
    H=analytic_hessian(e)
    return dict(
      curvature_selector=C,mode_kinetic=e.K,
      omega0_sq=e.omega0_sq,
      minimal_material_direct_ell_source=metric_source_grad_ell(e),
      linear_metric_trace_cross_Hessian=float(H[0,1]),
      static_response=stat,
      scalar_curvature_load=curvature_loads_scalar(
                   F0=1.,F_sigma=e.F_sigma,R=.08),
      excited_scalar_trace_mixing_at_ell_0=scalar_kinetic_coefficient_mixing(e,0),
      excited_scalar_trace_mixing_at_ell_0p1=scalar_kinetic_coefficient_mixing(e,.1),
      L_energy_at_nonzero_initial_state=independent_well_energy_at_W_zero(
                     e.K,e.omega0_sq,.1,.15),
      physical_constrained_gravity_solved=False,
      physical_primordial_excitation_amplitude_solved=False)

def report():
    on=SourceAudit()
    off=SourceAudit(w_eta=0.)
    near=SourceAudit(w_eta=1e-6)
    zero_k=SourceAudit(k=0.)
    return dict(
      scope="actual accepted action direct source/mixing audit, not a fitted no-CDM cosmological theory",
      cases=dict(selector_on=on_shell_charged_source_audit(on),
       W_near_zero=on_shell_charged_source_audit(near),
       W_exactly_zero=on_shell_charged_source_audit(off),
       k_zero=on_shell_charged_source_audit(zero_k)),
      conditional_new_operator_example=dict(
        extra_operator_not_in_accepted_action="-g_new*j_metric*ell",
        g_new=1.,
        extra_equilibrium_at_W_zero=candidate_source_solution_at_W_zero(off,1.)),
      verdict="minimal Jordan matter has no direct L source; curvature selector only linearly mixes L with metric while W>0; no persistent independent static L well proved after W=0",
      remaining=[
          "physical complete 4D mapper action and any permitted direct matter/structural operator",
          "nonzero primordial well on-shell geometry from exact Planck boundary with S_A≈1",
          "background geometric excitation abundance from independent initial state",
          "complete combined F/G2/G3/ADM observable Poisson and Weyl constraints",
          "no-CDM Planck equality acoustic and lensing comparison"
      ])
if __name__=="__main__":
    print(json.dumps(report(),indent=2))
