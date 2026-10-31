"""SDMC supporting-space trace: action-derived Jordan energy exchange test.

This file derives the on-shell reduced-mode balance including evolving
K=B f_L^2/3, c_L^2, m_L^2, W eta_K, a_J and a_s=S/S0.
It does NOT identify the mode energy with an independently conserved
cosmological rho_g, nor determine the physical Jordan lapse/shift constraints.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import exp
import json
import numpy as np

@dataclass(frozen=True)
class Epoch:
    """Reduced supporting-space coefficients and Jordan-time derivatives.

    Every dot is d/d t_c. a_s follows structural S/S0, not Jordan a_J.
    H_s = d ln a_s/dt_c = p_v H_J, where p_v is the structural tracker.
    """
    a_j: float
    a_s: float
    H_j: float
    H_s: float
    K: float
    c2: float
    mass2: float
    weta: float
    Kdot: float = 0.
    c2dot: float = 0.
    mass2dot: float = 0.
    wetadot: float = 0.

    def __post_init__(self):
        if self.a_j<=0 or self.a_s<=0 or self.K<=0:
            raise ValueError("a_j, a_s, K must be positive")
        if self.c2<0 or self.mass2<0 or self.weta<0:
            raise ValueError("requires nonnegative stability coefficients")
        for v in vars(self).values():
            if not np.isfinite(v): raise ValueError("nonfinite epoch quantity")
    @property
    def p_v(self):
        if self.H_j == 0: raise ValueError("p_v undefined for H_j=0")
        return self.H_s/self.H_j

def mode_parameters(k, e):
    """omega^2, C, derivatives of omega^2 and C in Jordan time."""
    if k<0 or not np.isfinite(k): raise ValueError("k>=0")
    r=k*k/e.a_j**2
    omega2=e.mass2+e.c2*r
    C=4*e.weta*k**4/e.a_s**4
    omega2dot=e.mass2dot+r*(e.c2dot-2*e.H_j*e.c2)
    Cdot=4*k**4/e.a_s**4*(e.wetadot-4*e.H_s*e.weta)
    return omega2,C,omega2dot,Cdot

def mode_kinematics(ell,ell_dot,u,u_dot,k,e):
    """On-shell acceleration; u=(delta ln M)-Phi_J is an external driver."""
    omega2,C,_,_=mode_parameters(k,e)
    q=ell/3+u
    acceleration=-(3*e.H_j+e.Kdot/e.K)*ell_dot-omega2*ell-2*C*q/(3*e.K)
    return dict(q=q,ell_ddot=acceleration,q_dot=ell_dot/3+u_dot)

def reduced_energy_and_pressure(ell,ell_dot,u,k,e):
    """Reduced mode energy and restricted isotropic Jordan scale pressure.

    P obtained by varying Jordan a at fixed structural S and field
    profiles, NOT from the full mixed ADM+Horndeski metric variation.
    """
    omega2,C,_,_=mode_parameters(k,e)
    q=ell/3+u
    kin=.5*e.K*ell_dot**2
    mass=.5*e.K*e.mass2*ell**2
    grad=.5*e.K*e.c2*(k/e.a_j)**2*ell**2
    curv=C*q*q
    E=kin+mass+grad+curv
    P=kin-mass-grad/3-curv
    return dict(E=E,P_restricted=P,kinetic=kin,mass=mass,
                gradient=grad,curvature=curv)

def energy_rate_on_shell(ell,ell_dot,u,u_dot,k,e):
    """Derivative of reduced E using Euler-Lagrange trace equation.

    EXACT for the quadratic reduced mode with prescribed external u and
    specified time-dependent coefficients; does not require dust closure.
    """
    omega2,C,omega2dot,Cdot=mode_parameters(k,e)
    q=ell/3+u
    return (-3*e.H_j*e.K*ell_dot**2
       + .5*e.Kdot*(omega2*ell**2-ell_dot**2)
       + .5*e.K*omega2dot*ell**2
       + Cdot*q*q+2*C*q*u_dot)

def jordan_energy_exchange_residual(ell,ell_dot,u,u_dot,k,e):
    """Jred=E_dot+3 H_J(E+P_red); Jred=0 is reduced conservation.

    The physical J_g can differ due to remaining field equations,
    ADM constraints and metric/scalar mixing.
    """
    quantities=reduced_energy_and_pressure(ell,ell_dot,u,k,e)
    energy_rate=energy_rate_on_shell(ell,ell_dot,u,u_dot,k,e)
    residual=energy_rate+3*e.H_j*(quantities["E"]+quantities["P_restricted"])
    return residual

def jordan_exchange_analytic(ell,ell_dot,u,u_dot,k,e):
    """Analytical simplified exchange law, exact in reduced model:

    J_red = .5 Kdot (omega2 ell² - ell_dot²)
          + .5 K (mass2dot + c2dot k²/a_J²) ell²
          + Cdot (ell/3+u)² + 2 C (ell/3+u) u_dot.
    """
    omega2,C,_,Cdot=mode_parameters(k,e)
    q=ell/3+u
    return (.5*e.Kdot*(omega2*ell**2-ell_dot**2)
          +.5*e.K*(e.mass2dot+e.c2dot*k*k/e.a_j**2)*ell**2
          +Cdot*q*q+2*C*q*u_dot)

def required_driver_rate_for_zero_exchange(ell,ell_dot,u,k,e):
    """Formal u_dot necessary for Jred=0 if C*q !=0.

    It is NOT an equation of motion for u; actual u_dot is determined
    jointly by Jordan metric constraints and the structural scalar.
    """
    omega2,C,_,Cdot=mode_parameters(k,e)
    q=ell/3+u
    if abs(C*q)<1e-14: raise ValueError("driver cannot control zero curvature coupling")
    numerator=(.5*e.Kdot*(omega2*ell**2-ell_dot**2)
       +.5*e.K*(e.mass2dot+e.c2dot*k*k/e.a_j**2)*ell**2
       +Cdot*q*q)
    return -numerator/(2*C*q)

def exchange_demo():
    """Compare curvature-only source and a formally balanced driver."""
    e=Epoch(a_j=1.,a_s=1.,H_j=.05,H_s=.075,K=1.,
            c2=1.,mass2=.2,weta=.04)
    ell=0.; vel=0.;u=.3; k=2.
    t=reduced_energy_and_pressure(ell,vel,u,k,e)
    J=jordan_energy_exchange_residual(ell,vel,u,0.,k,e)
    udot=required_driver_rate_for_zero_exchange(ell,vel,u,k,e)
    J_bal=jordan_energy_exchange_residual(ell,vel,u,udot,k,e)
    out={
        "scope":"reduced action-only trace energy; not global rho_g",
        "structural_tracker_p_v":e.p_v,
        "curvature_mode_energy":t["curvature"],
        "curvature_only_restricted_pressure":t["P_restricted"],
        "unforced_exchange_per_H_E":J/(e.H_j*t["E"]),
        "formal_required_driver_u_dot_for_Jred_zero":udot,
        "balanced_exchange_residual":J_bal,
        "not_yet_derived":["actual driver evolution","full Jordan F and X constraints",
            "physical geometric stress tensor","global cold-like pressureless fixed point",
            "independent source normalization","CMB and weak lensing"]
    }
    return out

if __name__=="__main__":
    print(json.dumps(exchange_demo(),indent=2))
