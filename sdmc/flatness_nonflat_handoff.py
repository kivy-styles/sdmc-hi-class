"""SDMC Planck-to-radiation flatness: non-flat classical constraint audit.

NO CDM sector is involved. Uses the original structural-to-Jordan
flatness question, not geometric-gravity abundance.

Part A: homogeneous spatial curvature in candidate restricted
UV + Jordan ADM constraint, *after* earliest classical start S_A≈1.
The exact Planck boundary is prethermal; no Jordan FLRW constraint
is imposed at S_P=1.

Manuscript B canonical ST lapse constraint (spatially flat)
  3F H² = rho_b+rho_r+ k1 v²/2+3k2 v^4/4+V
                -3H F_sigma v+3g H v³-.5g_sigma v^4 +rho_UV.
On non-flat FLRW (same action convention), gravitational term
  3F(H²+K/a²).
For exactly locked uniform structural gamma=M²h and an
EXPLICIT fixed-lapse, fixed-a_s positive operator
  U_K=6W eta_K K²/a_s^4,
the restricted lapse constraint becomes
  3F(H²+x)=E0+lambda*x²,
  x=K/a_J², lambda=6 W eta_K/M^4, a_s=M*a_J.
The term U_K and its simple lapse variation require W,eta,M
to be independent of lapse and no further UV off-shell terms.
A complete F/G2/G3+spatial parent action could alter this result.

Given independently solved E0, F,H, this quadratic can calculate
allowed x but cannot predict E0 and H from the same constraint,
and the constraint can have multiple roots. To be a physical
FLRW solution, K=a_J²*x must be fixed during the classical epoch.
Finite curvature selection from exact prethermal branch requires
additional transition matching and dynamical equations.

Part B: separately conserved transfer-only two-component toy
that covers a stiff to radiation-like local classical transition
and lets K vary across different *solutions*, never in time.
It deliberately includes no inflation or independently varying K.
X evolves from w≈1 to w≈1/3 and transfers energy to radiation:
 d rho_X/dln a = -3(1+w_X)rho_X-Gamma rho_X,
 d rho_r/dln a = -4rho_r+Gamma rho_X,
 Gamma=Q/(H rho_X)>=0.
Its Jordan toy Friedmann equation is
 H² = rho_total/(3F)-K/a²,
with F=1/3 in normalized units. H is solved self-consistently
for each K in this simple fixed-F toy.
This is NOT the SDMC full late300 scalar, UV, or thermalization
solution. Gamma(x) and w_X(x) are illustrative free functions.

In this classical decelerating toy, |Omega_K| increases;
nonzero Jordan K is not damped merely by entropy-producing
inter-component transfer. A possible prethermal spatial
handoff can still set a small curvature at Jordan emergence.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import json
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class ClassicalTransfer:
    x_start:float=0.
    x_end:float=4.
    initial_x_density:float=1.
    initial_rad_density:float=1e-8
    center_w:float=1.
    slope_w:float=6.
    gamma_max:float=5.
    gamma_center:float=1.2
    gamma_width:float=.5
    F:float=1/3

    def __post_init__(self):
        pars=np.array(list(asdict(self).values()),dtype=float)
        if not np.isfinite(pars).all() or self.x_end<=self.x_start or \
           self.initial_x_density<=0 or self.initial_rad_density<0 or \
           self.slope_w<=0 or self.gamma_width<=0 or self.gamma_max<0 or \
           self.F<=0: raise ValueError("positive stable toy transfer data required")

def wx(x,p:ClassicalTransfer):
    # (1+r)/(1+3r): manuscript early reduced G2 kinetic w.
    r=np.exp(np.clip(p.slope_w*(x-p.center_w),-100,100))
    return float((1+r)/(1+3*r))

def gamma(x,p:ClassicalTransfer):
    return float(p.gamma_max*np.exp(-.5*((x-p.gamma_center)/p.gamma_width)**2))

def transfer_rhs(x,y,p:ClassicalTransfer):
    X,R=y
    Qh=gamma(x,p)*X
    return np.array([-3*(1+wx(x,p))*X-Qh,-4*R+Qh])

def classical_transfer(p=ClassicalTransfer(),positions=None):
    if positions is None:
        positions=[0,.4,.8,1.2,2,3,4]
    xx=np.array(positions,dtype=float)
    if np.any(np.diff(xx)<=0) or min(xx)<p.x_start or max(xx)>p.x_end:
        raise ValueError("ordered sample inside classical domain required")
    sol=solve_ivp(lambda x,y:transfer_rhs(x,y,p),
         (p.x_start,p.x_end),
         [p.initial_x_density,p.initial_rad_density],
         dense_output=True,method="DOP853",rtol=2e-11,atol=2e-13,
         max_step=.03)
    if not sol.success:raise RuntimeError(sol.message)
    out=[]
    for x in xx:
        X,R=sol.sol(float(x))
        if min(X,R)<-1e-11:raise RuntimeError("negative transfer densities")
        E=float(X+R)
        out.append(dict(x=float(x),a=float(np.exp(x)),
            rho_X=float(X),rho_r=float(R),rho_total=E,
            effective_w=float((wx(x,p)*X+R/3)/E),
            scalar_w=wx(x,p),gamma_over_H=gamma(x,p)))
    return out

def curved_classical_history(K,p=ClassicalTransfer(),positions=None):
    """Cosmic expansion for *fixed* K and toy density evolution."""
    states=classical_transfer(p,positions)
    records=[]
    for s in states:
        a=s["a"]
        rho=s["rho_total"]
        curv=K/a**2
        H2=rho/(3*p.F)-curv
        if H2<=0:
            raise ValueError("closed model turns around before test endpoint")
        H=np.sqrt(H2)
        Omega=-curv/H2
        # In fixed-F GR toy, acceleration q=(rho+3P)/(6F H²).
        q=(rho*(1+3*s["effective_w"])/(6*p.F*H2))
        records.append(dict(**s,H=float(H),H2=float(H2),
                 K=float(K),K_reconstructed=float(
                       (rho/(3*p.F)-H2)*a*a),
                 Omega_k=float(Omega),deceleration=float(q)))
    return records

def nonflat_restricted_constraint(F,H,E0,W,eta,aJ,M):
    """Return positive-operator restricted lapse-constrained K roots.

    x=K/aJ² and positive UV energy lambda x², not a unique
    primordial curvature prediction.
    """
    vals=np.array([F,H,E0,W,eta,aJ,M],dtype=float)
    if not np.isfinite(vals).all() or min(F,aJ,M)<=0 or \
         min(H,E0,W,eta)<0:
        raise ValueError("regular positive action branch required")
    lam=6*W*eta/M**4
    D=E0-3*F*H*H
    A=3*F
    if lam==0:
        x=D/A
        roots=[x]
        disc=A*A
    else:
        disc=A*A-4*lam*D
        if disc<0: return dict(discriminant=float(disc),
                        physical_roots=[],real_solution=False)
        sqrt=np.sqrt(disc)
        xsmall=2*D/(A+sqrt)
        xhigh=(A+sqrt)/(2*lam)
        roots=sorted([float(xsmall),float(xhigh)])
    return dict(F=F,H=H,E0=E0,lambda_curv=float(lam),
                residual_flat=float(D),
                discriminant=float(disc),
                roots_x=roots,
                roots_K=[float(x*aJ*aJ) for x in roots],
                real_solution=True,
                curvature_dynamically_selected=False)

def constraint_at_x(F,H,E0,W,eta,aJ,M,x):
    lam=6*W*eta/M**4
    return 3*F*(H*H+x)-E0-lam*x*x

def archive_clock_check():
    """Independent direct kinematic structural results from the manuscript."""
    NP=1.
    NC=2/3
    NR=3.41350846
    return dict(crossover_growth=(NP/NC)**2,
        thermal_handoff_suppression=(NC/NR)**2,
        total_structural_suppression=(NP/NR)**2,
        jordan_freeze_ratio_if_C1=6.35e-17,
        jordan_recomb_ratio_if_C1=2.10696e-4,
        jordan_present_ratio_if_C1=.0917953,
        caveat="These are archived flat-reference kinematics, not calculated UV-sourced nonflat K")

def transition_curvature_residual(N,p,chi,f,GF,beta_b=0.,gamma_g=0.):
    """Dimensionless nonflat Jordan Omega_k from the SDMC master constraint.

    Uses post-thermal classical quantities, never exact S=1 boundary:
      rho_v=rho_P/S², rho_X=chi*rho_v,
      f=rho_X/(rho_X+rho_r),
      H_J=N/(p*t_P*S),
      rho_P*t_P²=1/G.

    Omega_k,J=1-p²/(3 G F N²)*[chi/f+beta_b+Gamma_g].
    Here GF=G*F; F=1/(8piG) -> GF=1/(8pi).
    beta_b=rho_b/rho_v, Gamma_g=rho_g/rho_v, with rho_g
    from action-derived nonflat spatial stress (not inserted CDM).
    Flat-reference endpoint chi was *fit to this same closure*;
    using it is algebraic consistency, not independent prediction.
    """
    vals=np.asarray([N,p,chi,f,GF,beta_b,gamma_g],dtype=float)
    if not np.isfinite(vals).all() or min(N,p,f,GF)<=0 or f>1 or \
          min(chi,beta_b,gamma_g)<0:
        raise ValueError("positive, regular postthermal branch required")
    balance=chi/f+beta_b+gamma_g
    expected=3*GF*N*N/(p*p)
    Omega=1-balance/expected
    return dict(omega_k_J=float(Omega),
           balance=float(balance),
           required_flat_balance=float(expected),
           absolute_source_residual=float(expected-balance),
           flat_closure_exact=bool(abs(Omega)<2e-14),
           independently_derived_curvature=False)

def endpoint_rounding_test():
    import math
    N=3.41350846
    p=2.
    f=.01
    GF=1/(8*math.pi)
    flat_chi=(3*GF*N*N/p**2)*f
    return dict(rounded=transition_curvature_residual(
            N,p,chi=.003477,f=f,GF=GF),
       exact_flat=transition_curvature_residual(
            N,p,chi=flat_chi,f=f,GF=GF),
       flat_chi_required=float(flat_chi),
       rounded_chi=.003477,
       note="rounded endpoint values of flat-reference branch do not independently predict curvature")

def diagnostic():
    p=ClassicalTransfer()
    histories={}
    for name,K in [("flat",0.),("open",-1e-6),
                   ("closed",1e-6),("open_100times",-1e-4)]:
        try:
            h=curved_classical_history(K,p)
            histories[name]=dict(K=K,
                   initial_omega=h[0]["Omega_k"],
                   final_omega=h[-1]["Omega_k"],
                   growth=abs(h[-1]["Omega_k"]/h[0]["Omega_k"])
                    if K else None,
                   sample=h)
        except ValueError as e:
            histories[name]=dict(K=K,exception=str(e))
    return dict(scope="restricted exact constraint algebra + fully conservative decelerating fixed-F two-fluid toy; not full nonflat SDMC",
       clock=archive_clock_check(),
       radiation_endpoint=endpoint_rounding_test(),
       flat_residual_selector_on=nonflat_restricted_constraint(
              F=1/3,H=1,E0=1,W=1,eta=.01,aJ=1,M=1),
       mismatched_constraint=nonflat_restricted_constraint(
              F=1/3,H=1,E0=1.002,W=1,eta=.01,aJ=1,M=1),
       selector_off=nonflat_restricted_constraint(
              F=1/3,H=1,E0=1.002,W=0,eta=.01,aJ=1,M=1),
       histories=histories,
       missing=["prethermal 4D boundary constraint at S_P=1",
                "microscopic source and initial 3-curvature distribution",
                "Jordan radiation bath Q from full action",
                "full nonflat UV ADM lapse, mapper and spatial equations",
                "nonflat scalar-tensor ODE with kinetic braiding and density backreaction",
                "nonflat CMB and BAO propagation"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
