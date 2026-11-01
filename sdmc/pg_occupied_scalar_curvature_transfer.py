"""SDMC Part III: occupied-well curvature -> scalar -> geometric trace response.

PURPOSE
-------
At W=0, minimally coupled Jordan matter has no direct *linear*
material -> independent structural trace ell source on the L=0 branch.
However, the EXISTING scalar-dependent structural coefficients may
cause a PRE-EXCITED relative-strain well to respond indirectly:

  delta R_J  -- F_sigma R_J/2 --> delta sigma
           -- K_sigma(sigma) --> delta ell (if ell0_dot !=0).

This is a TWO-STAGE CONDITIONAL reduced model, not the full rolling
G3 ADM + radiation/matter constraints or an SDMC no-CDM fit.

DERIVATION
----------
For a real cosine structural mode (Manuscript B Part III Sec. 106):

  S_ell = int dt_c a_J^3 K_ell(sigma)/4
                 [dot A^2 - omega_L^2(a_J) A^2],
  K_ell(sigma)=K0 exp(beta sigma) >0,
  omega_L^2=m_L^2+c_L^2 k^2/a_J^2; W=0.

With x=ln a_J, H_J=H0 exp(-h x), sigma_bar=v_sigma*x:

  A0_xx +(3-h+beta*v_sigma) A0_x
           +(omega_L^2/H_J^2)A0=0.

Introduce sigma=sigma_bar+epsilon*s(x) and
A=A0+epsilon*dA+O(epsilon^2). At fixed background Jordan metric:

  dA_xx +(3-h+beta*v_sigma)dA_x
           +(omega_L^2/H_J^2)dA
      = -beta*s_x*A0_x.

The scalar disturbance s(x) is driven by an EXTERNALLY PRESCRIBED
Jordan spacetime Ricci perturbation delta R_J(x) through an
illustrative canonical subcase

  s_xx+(3-h)s_x+(omega_sigma^2/H_J^2)s
        = F_sigma*deltaR_J/(2 K_sigma H_J^2).

The full accepted rolling Horndeski G2/G3 scalar equation and physical
ADM lapse/shift constraints are NOT substituted by this oscillator.
The scalar source is used only as a mathematically transparent probe.

The pre-existing structural mode's instantaneous reduced mode energy:

  E_L=K_ell/4 [H_J^2 A_x^2 + omega_L^2 A^2].

Its first variation at epsilon=0 is

  dE_L =
      beta*s*E_L0
      +K_ell/2[H_J^2 A0_x dA_x+omega_L^2 A0*dA].

If A0=A0_x=0 at the initial state, the response dA and dE both
vanish identically, even for nonzero F_sigma*delta R_J. This is the
precise background-occupation gate missing from the empty-branch
direct-source conclusion.

A nonzero A0 does NOT establish a homogeneous isotropic cosmic pg
background; it is a conditional excited-mode or local-patch test.
The full physical source amplitude and Weyl/Poisson response have
NOT been derived.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class OccupiedWell:
    H0:float=.5
    h:float=1.5
    beta:float=.8
    K0:float=1.
    sigma_bar_x:float=.15
    m_L2:float=1.
    c_L2:float=.3
    k_L:float=.8
    m_sigma2:float=.6
    c_sigma2:float=.2
    k_sigma:float=.3
    K_sigma:float=1.
    F_sigma:float=.25
    Ricci_pulse:float=.12
    pulse_center:float=0.
    pulse_width:float=.22
    initial_strain:float=.18
    initial_strain_x:float=.0
    x_start:float=-.5
    x_end:float=.85

    def __post_init__(self):
        vs=np.asarray(list(vars(self).values()),dtype=float)
        if not np.isfinite(vs).all():raise ValueError("finite coefficients required")
        if min(self.H0,self.K0,self.m_L2,self.m_sigma2,
               self.K_sigma,self.pulse_width)<=0:
            raise ValueError("positive H, kinetic, mass and width required")
        if min(self.c_L2,self.k_L,self.c_sigma2,self.k_sigma)<0:
            raise ValueError("nonnegative mode gradients and wavenumbers")
        if self.x_end<=self.x_start:
            raise ValueError("ordered Jordan e-fold range")

def Hubble(e,x):
    return e.H0*np.exp(-e.h*x)

def omega_L2(e,x):
    return e.m_L2+e.c_L2*e.k_L**2*np.exp(-2*x)

def omega_sigma2(e,x):
    return e.m_sigma2+e.c_sigma2*e.k_sigma**2*np.exp(-2*x)

def imposed_Ricci(e,x):
    """External curvature pulse, NOT the self-consistent physical Ricci."""
    z=(x-e.pulse_center)/e.pulse_width
    return e.Ricci_pulse*np.exp(-.5*z*z)

def scalar_rhs(e,x,s,s_x):
    H=Hubble(e,x)
    return (-(3-e.h)*s_x-(omega_sigma2(e,x)/H**2)*s
            +e.F_sigma*imposed_Ricci(e,x)/(2*e.K_sigma*H**2))

def l_friction(e,sigma_x=0.):
    return 3-e.h+e.beta*(e.sigma_bar_x+sigma_x)

def integrated_linear_response(e,n=501):
    """Simultaneously solve externally driven scalar and tangent L system.

    State y=(s,s_x,A0,A0_x,dA,dA_x), all with zero perturbation
    initial data but optional preexisting structural A0.
    """
    if n<25:raise ValueError("n>=25")
    x=np.linspace(e.x_start,e.x_end,n)
    y0=np.asarray([0.,0.,e.initial_strain,e.initial_strain_x,0.,0.])
    def rhs(x,y):
        s,sx,A,Ax,d,dx=y
        H=Hubble(e,x)
        omg=omega_L2(e,x)/H**2
        g=l_friction(e)
        return [sx,scalar_rhs(e,x,s,sx),
                Ax,-g*Ax-omg*A,
                dx,-g*dx-omg*d-e.beta*sx*Ax]
    sol=solve_ivp(rhs,(x[0],x[-1]),y0,
       t_eval=x,method="DOP853",rtol=2e-11,atol=2e-12,max_step=.005)
    if not sol.success:raise RuntimeError(sol.message)
    s,sx,A,Ax,d,dx=sol.y
    H=Hubble(e,x)
    K=e.K0*np.exp(e.beta*e.sigma_bar_x*x)
    om=omega_L2(e,x)
    E=K/4*(H**2*Ax**2+om*A**2)
    dE=e.beta*s*E+K/2*(H**2*Ax*dx+om*A*d)
    direct=e.beta*s*E
    induced=dE-direct
    return dict(
      x=x,s=s,s_x=sx,A0=A,A0_x=Ax,
      delta_A=d,delta_A_x=dx,E0=E,delta_E=dE,
      delta_E_direct_K=direct,delta_E_dynamical=induced,
      curvature_input=np.asarray(imposed_Ricci(e,x)),
      K=K,
      scope="externally driven canonical scalar and preoccupied geometric trace at fixed Jordan background")

def integrate_nonlinear_perturbed_strain(e,s,sx,epsilon,x):
    """Independently solve FULL exponential K(sigma) kinetic friction.

    sigma=sigma_bar+epsilon*s, exact K'/K=beta*(v+epsilon*s_x).
    Prescribed s is obtained from same linear scalar response;
    no direct material-to-L force is introduced.
    """
    xx=np.asarray(x)
    if xx.ndim!=1 or len(xx)<2 or np.any(np.diff(xx)<=0):
        raise ValueError("ordered x")
    def rhs(t,y):
        sigma_x=float(np.interp(t,xx,sx))
        return [y[1],-l_friction(e,epsilon*sigma_x)*y[1]
          -omega_L2(e,t)/Hubble(e,t)**2*y[0]]
    sol=solve_ivp(rhs,(xx[0],xx[-1]),
          [e.initial_strain,e.initial_strain_x],
          t_eval=xx,method="DOP853",rtol=2e-11,atol=2e-12,
          max_step=.002)
    if not sol.success:raise RuntimeError(sol.message)
    A,Ax=sol.y
    K=e.K0*np.exp(e.beta*(e.sigma_bar_x*xx+epsilon*s))
    E=K/4*(Hubble(e,xx)**2*Ax**2+omega_L2(e,xx)*A**2)
    return dict(A=A,Ax=Ax,E=E)

def centered_fd_variation(e,epsilon=1e-3,n=501):
    """Compare coupled tangent variation against independent ±eps runs."""
    if epsilon<=0:raise ValueError("positive epsilon")
    out=integrated_linear_response(e,n=n)
    p=integrate_nonlinear_perturbed_strain(
        e,out["s"],out["s_x"],epsilon,out["x"])
    m=integrate_nonlinear_perturbed_strain(
        e,out["s"],out["s_x"],-epsilon,out["x"])
    fdA=(p["A"]-m["A"])/(2*epsilon)
    fdE=(p["E"]-m["E"])/(2*epsilon)
    dA=out["delta_A"];dE=out["delta_E"]
    errA=float(np.max(np.abs(fdA-dA)))
    errE=float(np.max(np.abs(fdE-dE)))
    return dict(
      max_abs_strain_linearization_error=errA,
      max_abs_energy_linearization_error=errE,
      max_induced_deltaA=float(np.max(np.abs(dA))),
      max_energy_response=float(np.max(np.abs(dE))),
      max_scalar_perturbation=float(np.max(np.abs(out["s"]))),
      max_curvature_pulse=float(np.max(np.abs(out["curvature_input"]))),
      max_occupied_energy=float(np.max(out["E0"])),
      end_occupied_energy=float(out["E0"][-1]),
      initial_strain=e.initial_strain,
      K_sigma_dependence_beta=e.beta,
      requires_nonzero_background_occupation=True,
      source_is_prescribed_not_physical_constraint=True)

def occupied_mode_continuity(e,x,A,Ax,s=0.,sx=0.):
    """Exact e-fold energy continuity identity for variable K(sigma).

    The physical *isotropic gradient* pressure in this reduced cosine
    oscillator is
      P=K/4[H² Ax²-m² A²-(G/3) A²],
    whereas the scalar transfer pressure factor is
      Pstar=K/4[H² Ax²-(m²+G) A²].
    For gamma=3-h+Kx/K, the exact W=0 oscillator obeys:
      E_x+3(E+P) = -(K_x/K)*Pstar.

    This algebra holds at every x even if sigma is time-dependent.
    A WKB mass oscillator with <Pstar>≈0 consequently has negligible
    averaged exchange from slow running K, despite its instantaneous
    nonzero scalar coupling.

    Input A,Ax is the exact solution of the corresponding background
    with specified scalar perturbation s and sx; all fields and metrics
    use the same Jordan x-coordinate.
    """
    H=Hubble(e,x)
    sigma=e.sigma_bar_x*x+s
    kappa=e.beta*(e.sigma_bar_x+sx)
    K=e.K0*np.exp(e.beta*sigma)
    G=e.c_L2*e.k_L**2*np.exp(-2*x)
    omg2=e.m_L2+G
    Axx=-(3-e.h+kappa)*Ax-omg2/H**2*A
    E=K/4*(H**2*Ax**2+omg2*A**2)
    P=K/4*(H**2*Ax**2-e.m_L2*A**2-G*A**2/3)
    Pstar=K/4*(H**2*Ax**2-omg2*A**2)
    # Independent differentiation of E, not using continuity formula.
    Hprime=-e.h*H
    E_x=kappa*E+K/4*(2*H*Hprime*Ax**2+
                     2*H**2*Ax*Axx-2*G*A**2+2*omg2*A*Ax)
    J_k=-kappa*Pstar
    return dict(E=float(E),P=float(P),Pstar=float(Pstar),
      E_x=float(E_x),J_over_H=float(J_k),
      balance_residual=float(E_x+3*(E+P)-J_k),
      K_log_derivative=float(kappa))

def response_comparison():
    occupied=OccupiedWell()
    empty=OccupiedWell(initial_strain=0.,initial_strain_x=0.)
    no_coeff=OccupiedWell(beta=0.)
    no_curvature=OccupiedWell(Ricci_pulse=0.)
    no_scalar_coupling=OccupiedWell(F_sigma=0.)
    def stats(e):
        r=integrated_linear_response(e,n=261)
        return dict(
          max_scalar_perturbation=float(np.max(abs(r["s"]))),
          max_perturbed_trace=float(np.max(abs(r["delta_A"]))),
          max_perturbed_geometric_energy=float(np.max(abs(r["delta_E"]))),
          energy_nonzero=bool(np.max(r["E0"])>0))
    return dict(
      scope="conditional curvature-to-scalar-to-occupied-L response; full SDMC constraints remain open",
      occupied=stats(occupied),
      empty=stats(empty),
      beta_zero=stats(no_coeff),
      Ricci_zero=stats(no_curvature),
      Fsigma_zero=stats(no_scalar_coupling),
      finite_difference_crosscheck=centered_fd_variation(occupied,n=301),
      continuity_check=occupied_mode_continuity(occupied,x=.1,A=.08,Ax=-.2),
      key_result="F_sigma deltaR can modulate existing L through K_sigma; it cannot excite a zero-amplitude L through the homogeneous quadratic action",
      physical_constraints_not_yet_solved=[
        "derive deltaR_J and sigma from full rolling F/G2/G3 matter ADM constraints",
        "include sigma-L mutual backreaction and complete gauge-invariant perturbations",
        "derive primordial structural occupation spectrum and amplitude at S_A≈1",
        "check Jordan isotropic average and physical rho_g, P_g, J_g",
        "test observable Poisson/Weyl response, equality, CMB and weak lensing"
      ])
if __name__=="__main__":
    print(json.dumps(response_comparison(),indent=2))
