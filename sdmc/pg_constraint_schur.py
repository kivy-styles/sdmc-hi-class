"""SDMC supporting-space/Jordan metric constraint-block test.

Eliminates scalar structural trace ell from its restricted quadratic action
and exposes what this supplies to the Jordan potential Phi and mapping
m=delta ln(M). Other Jordan/active-scalar Hessian blocks are not fixed by
the supporting-space sector.

An optional positive-stiffness toy Jordan Hessian A,D tests the Schur
complement, without claiming a physical mu, Sigma or CDM replacement.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import sqrt
import json
import numpy as np

@dataclass(frozen=True)
class TraceAction:
    K:float=1.
    c2:float=1.
    mass2:float=.3
    weta:float=.05
    a_j:float=1.
    a_s:float=1.

    def __post_init__(self):
        v=(self.K,self.c2,self.mass2,self.weta,self.a_j,self.a_s)
        if not np.isfinite(v).all(): raise ValueError("nonfinite coefficients")
        if self.K<=0 or self.a_j<=0 or self.a_s<=0:
            raise ValueError("positive K and spatial scales required")
        if min(self.c2,self.mass2,self.weta)<0:
            raise ValueError("nonnegative reduced quadratic coefficients required")

def trace_stiffness(k,act):
    """Spatial-k in action units; physical k units are not yet fixed."""
    if k<0 or not np.isfinite(k): raise ValueError("finite k>=0 required")
    omega2=act.mass2+act.c2*k*k/act.a_j**2
    C=4*act.weta*k**4/act.a_s**4
    Lam=2*C/(9*act.K)
    return dict(omega2=omega2,C=C,Lambda=Lam,
                trace_frequency_squared=omega2+Lam)

def effective_kernel(k,act):
    """E_eff = Keff (m-Phi_J)^2; trace eliminated at fixed drivers."""
    v=trace_stiffness(k,act)
    if v["trace_frequency_squared"]==0:return 0.
    return v["C"]*v["omega2"]/v["trace_frequency_squared"]

def quasi_static_frequency_ratio(k,act,H_j):
    """Necessary ratio omega_eff/H_J for static trace elimination."""
    if H_j<=0 or not np.isfinite(H_j): raise ValueError("H_j>0")
    return sqrt(trace_stiffness(k,act)["trace_frequency_squared"])/H_j

def full_static_hessian(k,act,A,D):
    """Toy Hessian [Phi,m,ell] for E=A Phi²/2 +D m²/2+
       K omega² ell²/2 + C(m-Phi+ell/3)² -J Phi.
    A and D are not supplied by the reduced structural action.
    D stands for an existing mapper/scalar stiffness, not a new species.
    """
    if A<0 or D<0: raise ValueError("nonnegative toy Hessian stiffness")
    v=trace_stiffness(k,act)
    baseline=np.diag([A,D,act.K*v["omega2"]])
    relative=np.array([-1.,1.,1./3])
    return baseline+2*v["C"]*np.outer(relative,relative)

def schur_complement_trace(k,act,A,D):
    """Exact Hessian after eliminating ell."""
    H=full_static_hessian(k,act,A,D)
    if H[2,2]<=0:
        if effective_kernel(k,act)==0:return np.diag([A,D])
        raise ValueError("singular structural trace")
    return H[:2,:2]-np.outer(H[:2,2],H[2,:2])/H[2,2]

def analytic_effective_hessian(k,act,A,D):
    Ke=effective_kernel(k,act)
    return np.diag([A,D])+2*Ke*np.array([[1.,-1.],[-1.,1.]])

def generic_schur(H,retained):
    """Eliminate specified block of a real symmetric quadratic form.

    The full physical Hessian must be derived independently.
    Gauge/constraint degeneracies require special handling; do not
    replace an undefined block inverse by a pseudoinverse.
    """
    H=np.asarray(H,dtype=float)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or not np.allclose(
            H,H.T,rtol=1e-10,atol=1e-10):
        raise ValueError("real symmetric square Hessian required")
    if not np.isfinite(H).all():raise ValueError("nonfinite Hessian")
    keep=list(retained)
    if not keep or len(set(keep))!=len(keep) or any(
            x<0 or x>=len(H) for x in keep):
        raise ValueError("retained indices invalid")
    take=[i for i in range(len(H)) if i not in keep]
    if not take:return H[np.ix_(keep,keep)].copy()
    Hcc=H[np.ix_(take,take)]
    Hkc=H[np.ix_(keep,take)]
    return H[np.ix_(keep,keep)]-Hkc@np.linalg.solve(Hcc,Hkc.T)

def toy_jordan_response(k,act,A,D,J=1.):
    """Positive-stiffness toy Phi response; NOT cosmological mu(k,a).

    A physical gravitational constraint can have an indefinite Hessian
    signature and all the F/G2/G3 ADM blocks are absent from this toy.
    """
    if A<=0 or D<0: raise ValueError("A>0,D>=0 required")
    Ke=effective_kernel(k,act)
    deltaA=2*Ke*D/(D+2*Ke) if D+2*Ke>0 else 0.
    Aeff=A+deltaA
    phi=J/Aeff
    m=(2*Ke/(D+2*Ke)*phi) if D+2*Ke>0 else 0.
    v=trace_stiffness(k,act)
    ell=(-3*v["Lambda"]/v["trace_frequency_squared"]*(m-phi)
         if v["trace_frequency_squared"]>0 else 0.)
    return dict(k=k,baseline_response=J/A,phi=phi,m=m,trace_ell=ell,
       effective_kernel=Ke,delta_jordan_stiffness=deltaA,
       toy_screening_ratio=A/Aeff,
       structural_adjustment_fraction=2*Ke/(D+2*Ke)
           if D+2*Ke>0 else 0.)

def direct_full_hessian_solution(k,act,A,D,J=1.):
    H=full_static_hessian(k,act,A,D)
    return np.linalg.solve(H,np.array([J,0.,0.]))

def diagnostics():
    act=TraceAction(K=1.,c2=.8,mass2=.1,weta=.03)
    levels=[.01,.1,1.,10.,100.]
    return {
       "scope":"Reduced action Schur complement and positive-stiffness Jordan toy only",
       "trace_coefficients":vars(act),
       "k_sweep":[toy_jordan_response(k,act,A=k*k,D=.5) for k in levels],
       "limits":{
           "k_to_zero":"Keff ~ 4 W eta k^4/a_s^4 if omega2 nonzero",
           "k_to_infinity":"Keff ~ (9K/2)c2 k^2/a_j^2 if c2,Weta>0",
           "free_mapper_D_zero":"relative strain relaxes and induces no Phi stiffness",
           "stiff_mapper_D_large":"additional positive stiffness tends to 2 Keff"},
       "missing":["complete F/G2/G3 Hessian","ADM lapse and shift constraints",
           "physical gravitational mu and Sigma","ghost and sound-speed analysis",
           "independent rho_g background normalization","no-CDM Boltzmann likelihood"]}
if __name__=="__main__":
    print(json.dumps(diagnostics(),indent=2))
