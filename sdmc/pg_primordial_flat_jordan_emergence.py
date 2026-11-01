"""SDMC Part III: primordial-flat supporting space at emergent Jordan B=0.

Tests the user's clue against Manuscript B Parts I-II-III:
  gamma_ij = a_s² exp(2q) delta_ij,
  h_ij     = a_J² exp(-2Phi_J)delta_ij,
  q_k      = m_k - Phi_J(k) + ell_k/3
            = u_k + ell_k/3.
A nonzero k-mode of gamma is intrinsically flat to first order
iff q_k=0, or ell=-3u. This is a compatible INITIAL condition,
but NOT a consequence of B=0.

The manuscript illustrates B∈[0,1], W=4B(1-B).
In the reduced positive structural action,
  K=B f_L²/3,
  C=4 eta_K W k^4/a_s^4,
  E=K omega0² ell²/4+(C/2)(u+ell/3)².
  dE/dell=K omega0² ell/2 +C(u+ell/3)/3.
  Lambda=2C/(9K) for B>0.
When B=0, W=0, K=C=0, E=0 for every finite ell,
so the classical relative-strain equation loses kinetic rank.
Although Lambda(B→0+) can be finite, it is mathematically
invalid to use that normalized equation AT B=0.
The finite right-limit ell_eq does not select initial ell.
A regular prethermal parent is needed to specify shape data.

For B>0, the metric potential u can drive ell through W, if and only
if the Jordan material/metric constraints yield nonzero u.
At mature B=1, W=0, the static response ell_eq=0, and
gamma inherits the local Jordan curvature q=u.

Nothing here derives extra CDM-like gravity, independent abundance,
primordial curvature-zero state, or a Planck-era mode equation.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np

@dataclass(frozen=True)
class FlatHandoff:
    fL2:float=3.
    eta:float=.04
    k:float=1.5
    a_s:float=1.
    omega0_sq:float=1.8
    u:float=.08
    def __post_init__(self):
        v=np.array(list(vars(self).values()),dtype=float)
        if not np.isfinite(v).all() or min(self.fL2,self.a_s,self.omega0_sq)<=0:
            raise ValueError("finite positive kinetic and restoring parameters")
        if self.eta<0 or self.k<0:raise ValueError("nonnegative selector and k")

def coefficient(e,B):
    """B=0 is accepted as a DEGENERATE boundary, never divide by B."""
    if not np.isfinite(B) or not 0<=B<=1:
        raise ValueError("0<=B<=1 required")
    W=4*B*(1-B)
    K=B*e.fL2/3
    C=4*e.eta*W*e.k**4/e.a_s**4
    L=None if B==0 else 2*C/(9*K)
    return dict(B=B,W=W,K=K,C=C,Lambda=L,
       kinetic_rank_nonzero=bool(K>0))

def energy(e,B,ell,u=None):
    c=coefficient(e,B)
    v=e.u if u is None else u
    return c["K"]*e.omega0_sq*ell*ell/4+c["C"]*(v+ell/3)**2/2

def source_and_hessian(e,B,ell,u=None):
    c=coefficient(e,B)
    v=e.u if u is None else u
    return dict(
       strain_force=c["K"]*e.omega0_sq*ell/2+
                   c["C"]*(v+ell/3)/3,
       H_ell_ell=c["K"]*e.omega0_sq/2+c["C"]/9,
       H_u_ell=c["C"]/3,
       selector_source=-c["C"]*v/3,
       physical_matter_ell_source_direct=0.)
def gamma_curvature_mode(e,ell,u=None):
    v=e.u if u is None else u
    q=v+ell/3
    return dict(q=q,delta_R3=4*e.k**2*q/e.a_s**2,
       exact_linear_flat=bool(abs(q)<1e-13),
       delta_Jordan_phi_alone_not_extra_geometry=True)

def matched_flat_initial_shape(e,u=None):
    """Geometry alone: R3_gamma,k=0 for k>0 requires ell=-3u."""
    if e.k<=0:
        raise ValueError("k=0 is a homogeneous gauge/scale mode, no local R3 constraint")
    v=e.u if u is None else u
    return -3*v

def static_positive_B(e,B,u=None):
    """Valid only in the regular 0<B<=1 classical relative sector."""
    c=coefficient(e,B)
    if B<=0:raise ValueError("no unique B=0 equilibrium: kinetic rank zero")
    v=e.u if u is None else u
    lam=c["Lambda"]
    ell=-3*lam*v/(e.omega0_sq+lam)
    return dict(B=B,ell=ell,q=v+ell/3,
       gamma_R3=gamma_curvature_mode(e,ell,u=v)["delta_R3"],
       field_energy=energy(e,B,ell,u=v),
       strain_force=source_and_hessian(e,B,ell,u=v)["strain_force"],
       selector_energy=c["C"]/2*(v+ell/3)**2)

def right_limit_B_zero(e):
    """Algebraic B→0+ limit, *not* a B=0 equation or initial datum."""
    if e.fL2<=0:raise ValueError("positive fL2 required")
    Lambda_lim=(32*e.eta*e.k**4)/(3*e.fL2*e.a_s**4)
    ell=-3*Lambda_lim*e.u/(e.omega0_sq+Lambda_lim)
    return dict(Lambda_limit=Lambda_lim,
       conditional_stationary_strain_right_limit=ell,
       conditional_stationary_curvature_right_limit=(
          gamma_curvature_mode(e,ell)["delta_R3"]),
       initial_strain_determined_at_B_zero=False,
       B0_lagrangian_vanishes=True)

def activation_energy_derivative(e,B,ell,u=None):
    """dE/dB holding (ell,u) fixed.

    Shows how a coupling handoff can do work on an occupied geometric
    state; a parent covariant action must supply compensating energy.
    """
    coefficient(e,B)
    v=e.u if u is None else u
    D=e.fL2/3
    A=16*e.eta*e.k**4/e.a_s**4
    return D*e.omega0_sq*ell*ell/4 + A*(1-2*B)*(v+ell/3)**2/2

def constraint_rank_check(e,B):
    c=coefficient(e,B)
    return dict(B=B,kinetic_hessian_ell=c["K"],
       static_Hessian_ell=c["K"]*e.omega0_sq/2+c["C"]/9,
       canonical_momentum_for_fixed_velocity=c["K"]*.1,
       no_classical_differential_equation_at_B0=bool(B==0),
       W=c["W"])

def diagnostic():
    e=FlatHandoff()
    nums=[0,1e-5,.01,.1,.5,.9,1.]
    vals=[]
    for B in nums:
        c=coefficient(e,B)
        row=dict(B=B,W=c["W"],K=c["K"],C=c["C"],
          mixed_Hessian=c["C"]/3,kinetic_rank_nonzero=c["kinetic_rank_nonzero"])
        if B:
            s=static_positive_B(e,B)
            row.update(ell_equilibrium=s["ell"],
                       structural_R3=s["gamma_R3"],
                       geometric_energy=s["field_energy"])
        else:
            row.update(ell_equilibrium=None,structural_R3=None,
                       geometric_energy=None)
        vals.append(row)
    flat=matched_flat_initial_shape(e)
    return dict(
       scope="candidate SDMC B=0 Jordan emergence: analytic singular rank and geometrical flat-state matching",
       material_gravitates_through_Jordan_u_not_inserted_direct_L_source=True,
       structural_flat_initial=dict(u=e.u,ell=flat,
                  structural_ricci=gamma_curvature_mode(e,flat)["delta_R3"]),
       formal_Bzero=constraint_rank_check(e,0),
       formal_right_limit=right_limit_B_zero(e),
       activation_work_at_flat_initial=[
            dict(B=B,dE_dB=activation_energy_derivative(e,B,flat))
                     for B in (0.,.2,.5,1.)],
       activation_rows=vals,
       physical_initial_state_from_Bzero_undetermined=True,
       source_W_zero_both_Bzero_and_Bone=True,
       late_ordinary_Jordan_well_can_remain_curved=True,
       independent_extra_CDM_density_not_predicted=True)
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
