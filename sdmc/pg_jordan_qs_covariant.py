"""SDMC Jordan quasi-static constraint checkpoint from a covariant ACTION subclass.

Physical basis: S = int sqrt(-g)[F(sigma)R/2 + k1(sigma)Z - V(sigma)]
              + S_m + the *restricted static trace-eliminated*
                SDMC spatial supporting-space action,
Z=-g^ab grad_a sigma grad_b sigma/2.

This is a frozen-coefficient, weak-field, quasi-static subclass with
NO background kinetic braiding. The accepted late300 Horndeski action
has G3=g(sigma) Z and a rolling structural scalar; those terms are not
contained in these gravitational constraints.

Newtonian-gauge reduced spatial-gradient action per Jordan volume:
 L2 = F Q(Phi^2-2 Phi Psi)
    + F_sigma Q s (Psi-2 Phi) - (k1 Q+mass2)s^2/2
    - rho_delta Psi - K_geom (s/sigma-Phi)^2,
 Q = k_comoving^2/a_J^2, s=delta sigma.

The last kernel is independently derived by eliminating the structural
trace at fixed driver; it is not an effective CDM matter density.
Observable metric strengths require a background satisfying the action
and a full Horndeski/ADM perturbation hierarchy.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import sqrt
import json
import numpy as np
from sdmc.pg_constraint_schur import TraceAction, effective_kernel, trace_stiffness

@dataclass(frozen=True)
class JordanSubclass:
    F: float = 1.
    F_sigma: float = .2
    k1: float = 1.
    mass2: float = .1
    sigma: float = 1.
    F_ref: float = 1.
    a_j: float = 1.

    def __post_init__(self):
        if any(not np.isfinite(v) for v in vars(self).values()):
            raise ValueError("nonfinite scalar-tensor coefficient")
        if self.F<=0 or self.F_ref<=0 or self.sigma<=0 or self.a_j<=0:
            raise ValueError("F,F_ref,sigma and a_j must be positive")
        if self.k1<0 or self.mass2<0:
            raise ValueError("restricted healthy canonical branch requires k1,mass2>=0")

def static_quadratic_action(psi,phi,s,rho,k,bg,structural=None):
    """Lagrangian per Jordan physical volume, no time-gradient terms.

    Sign of matter term is -rho Psi. Vary Psi, Phi and s to solve
    physical static Einstein/sigma constraints in *this subclass*.
    """
    if k<=0: raise ValueError("k>0 required")
    Q=(k/bg.a_j)**2
    Kg=0. if structural is None else effective_kernel(k,structural)
    return (bg.F*Q*(phi**2-2*phi*psi)
            +bg.F_sigma*Q*s*(psi-2*phi)
            -.5*(bg.k1*Q+bg.mass2)*s*s
            -rho*psi
            -Kg*(s/bg.sigma-phi)**2)

def constraint_hessian(k,bg,structural=None):
    """Symmetric Hessian of static L2 in order (Psi,Phi,delta_sigma).

    H x = (rho,0,0). The signed gravitational constraint Hessian is
    indefinite; positive-definite toy screening arguments do not apply.
    """
    if k<=0:raise ValueError("k>0 required")
    Q=(k/bg.a_j)**2
    KG=0. if structural is None else effective_kernel(k,structural)
    ff=bg.F_sigma
    return np.array([
        [0.,                 -2*bg.F*Q,              ff*Q],
        [-2*bg.F*Q,          2*bg.F*Q-2*KG,
                            -2*ff*Q+2*KG/bg.sigma],
        [ff*Q,              -2*ff*Q+2*KG/bg.sigma,
                            -(bg.k1*Q+bg.mass2)-2*KG/bg.sigma**2]
    ],dtype=float)

def static_perturbation_solution(k,bg,rho=1.,structural=None):
    """Quasi-static Jordan potentials from restricted action.

    mu=Psi/Psi_GR, Sigma=(Phi+Psi)/(2 Psi_GR), using
    Psi_GR=-rho/(2 F_ref (k/a_J)^2).
    Neither is late300's full no-CDM gravitational response.
    """
    if not np.isfinite(rho) or rho==0:raise ValueError("finite nonzero source")
    H=constraint_hessian(k,bg,structural)
    x=np.linalg.solve(H,np.array([rho,0.,0.]))
    baseline=-rho/(2*bg.F_ref*(k/bg.a_j)**2)
    return dict(k=k,psi_J=float(x[0]),phi_J=float(x[1]),
       delta_sigma=float(x[2]),
       delta_F=float(bg.F_sigma*x[2]),
       mu_subclass=float(x[0]/baseline),
       Sigma_subclass=float((x[1]+x[0])/(2*baseline)),
       slip=float(x[1]-x[0]),
       slip_ratio=float(x[1]/x[0]) if x[0] else None,
       structural_kernel=0. if structural is None else effective_kernel(k,structural),
       determinant=float(np.linalg.det(H)))

def brans_dicke_mu_analytic(k,bg):
    """Exact restricted unbraided scalar-tensor quasistatic response.

    mu = Fref/F [1+F_sigma² Q/(2F(k1 Q+mass2)+3F_sigma² Q)]
    Sigma=Fref/F. If F_sigma!=0 and k1=mass2=0, massless scalar
    Brans-Dicke / f(R) limit gives mu=4/3*Fref/F.
    """
    if k<=0:raise ValueError("k>0")
    Q=(k/bg.a_j)**2
    denom=2*bg.F*(bg.k1*Q+bg.mass2)+3*bg.F_sigma**2*Q
    if denom==0:
        if bg.F_sigma==0: return bg.F_ref/bg.F
        raise ValueError("singular physical scalar constraint")
    return bg.F_ref/bg.F*(1+bg.F_sigma**2*Q/denom)

def horndeski_alpha_coefficients(F,F_sigma,H,sigma_dot,Z,
                                  k1,k2,g,g_sigma):
    """Background EFT alpha functions of accepted G2/G3/G4 action.

    G2=k1 Z+k2 Z²-V, G3=g Z, G4=F/2, G5=0,
    S3=-G3 Box sigma. From the EFT-of-Horndeski definitions:
    alpha_M = sigma_dot F_sigma/(H F)
    alpha_B = sigma_dot(2 Z g - F_sigma)/(H F)
    alpha_K = 2 Z/(H²F)*(k1+6k2 Z-4g_sigma Z)
               +12 sigma_dot Z g/(H F).
    These coefficients alone do not solve the full Jordan constraints.
    """
    if F<=0 or H==0 or not np.isfinite(
       [F,F_sigma,H,sigma_dot,Z,k1,k2,g,g_sigma]).all():
        raise ValueError("positive F, nonzero H and finite coefficients required")
    am=sigma_dot*F_sigma/(H*F)
    ab=sigma_dot*(2*Z*g-F_sigma)/(H*F)
    ak=2*Z/(H*H*F)*(k1+6*k2*Z-4*g_sigma*Z)+12*sigma_dot*Z*g/(H*F)
    D=ak+1.5*ab**2
    return dict(alpha_M=am,alpha_B=ab,alpha_K=ak,D=D,
                no_slip_residual=ab+2*am,
                no_slip_covariant_residual=Z*g+F_sigma/2)

def sweep():
    """Finite coefficient test, not a late300 parameter reconstruction."""
    bg=JordanSubclass()
    st=TraceAction(K=1.,c2=.8,mass2=.1,weta=.03,
                   a_j=1.,a_s=1.)
    ks=(.01,.1,.5,1.,2.,10.)
    rows=[]
    for k in ks:
        b=static_perturbation_solution(k,bg)
        y=static_perturbation_solution(k,bg,structural=st)
        rows.append(dict(k=k,mu_F_only=b["mu_subclass"],
            Sigma_F_only=b["Sigma_subclass"],
            mu_with_structural=y["mu_subclass"],
            Sigma_with_structural=y["Sigma_subclass"],
            structural_kernel=y["structural_kernel"]))
    return {
        "scope":"Physically signed quasistatic F-R/canonical-Z SUBCLASS + reduced structural kernel",
        "not_theory":"Not the full rolling G3, k2, late300 No-Slip model",
        "background":vars(bg),"structural":vars(st),
        "model_rows":rows,
        "missing":["G3 kinetic braiding and time-dependent background",
                   "full lapse/shift and scalar Horndeski equations",
                   "structural branch physical homogeneous rho_g",
                   "cosmic initial conditions and radiation",
                   "physical lensing and no-CDM cosmological likelihood"]}
if __name__=="__main__":
    print(json.dumps(sweep(),indent=2))
