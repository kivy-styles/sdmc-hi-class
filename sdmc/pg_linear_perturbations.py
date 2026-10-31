"""SDMC pg: action trace and conditional Jordan perturbations.

Two strictly separate tracks:
(1) linearized quadratic supporting-space trace response from Manuscript B;
(2) an *assumed* pressureless effective pg stress as a GR/subhorizon
    code-control. The accepted spatial action has not derived this closure.
The tests are not a no-CDM Boltzmann/lensing validation.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import exp, log, sqrt
import json
import numpy as np
from scipy.integrate import solve_ivp

C_KMS=299792.458

@dataclass(frozen=True)
class Background:
    """Flat constant-F Jordan background; rho_X is smooth constant for control."""
    h: float
    omega_b: float
    omega_pg: float
    omega_r: float
    omega_x: float

    def __post_init__(self):
        if self.h<=0 or any(v<0 for v in
           (self.omega_b,self.omega_pg,self.omega_r,self.omega_x)):
            raise ValueError("positive h and nonnegative densities required")
        if abs(self.omega_b+self.omega_pg+self.omega_r+self.omega_x-1)>1e-8:
            raise ValueError("flat background requires sum Omega=1")

    def e2(self,a):
        if a<=0: raise ValueError("a>0")
        return (self.omega_b+self.omega_pg)/a**3+self.omega_r/a**4+self.omega_x

    def hconf(self,a):
        """a H_J/c in Mpc^-1, where H0=100*h km/s/Mpc."""
        return a*100*self.h/C_KMS*sqrt(self.e2(a))

    def dln_h_dln_a(self,a):
        return -0.5*(3*(self.omega_b+self.omega_pg)/a**3+4*self.omega_r/a**4)/self.e2(a)

    def frac_b(self,a): return self.omega_b/a**3/self.e2(a)
    def frac_pg(self,a): return self.omega_pg/a**3/self.e2(a)

def uniform_nested_sphere_mass_ratio(ru_over_r):
    """Uniform nested sphere ONLY, not a general Vu/Vs or rho_pg formula."""
    if not 0<=ru_over_r<1: raise ValueError("0 <= r_u/R < 1 required")
    return 1/(1-ru_over_r**3)

def gamma_for_dust(a,s_over_s0,gamma0):
    """Conditional Gamma_g target = Gamma_g0 (S/S0)^2/a^3.

    Algebraic reconstruction imposing a^-3, NOT structural action evolution.
    """
    if a<=0 or s_over_s0<=0 or gamma0<0: raise ValueError("invalid normalization")
    return gamma0*s_over_s0**2/a**3

def dust_conformal_derivatives(hconf,k_per_mpc,delta,theta,phi_prime,psi):
    """Exact separately conserved cold source in conformal Newtonian gauge.

    delta'= -theta+3Phi', theta'= -Hconf*theta+k^2Psi.
    Phi/Psi must be supplied by independent Jordan metric constraints.
    """
    if hconf<=0 or k_per_mpc<0: raise ValueError("invalid hconf or k")
    return (-theta+3*phi_prime, -hconf*theta+k_per_mpc**2*psi)


def dust_synchronous_derivatives(hconf,theta,htrace_prime):
    """Exact separately conserved dust in scalar synchronous gauge.

    delta_g' = -theta_g - h_sync'/2; theta_g' = -Hconf theta_g.
    htrace_prime is the synchronous spatial-metric trace derivative,
    NOT the Jordan Hubble rate; do NOT impose theta_g=0 without a
    justified pg-comoving gauge choice for the completed geometric action.
    """
    if hconf<=0: raise ValueError("hconf>0")
    return (-theta-0.5*htrace_prime, -hconf*theta)

def trace_stiffness(k,a_j,a_s,bf2,c_l2,m_l2,w_eta):
    """Quadratic structural trace, fixed coefficients and dimensionless k.

    omega0^2 = m_L^2 + c_L^2 k^2/a_J^2,
    Lambda = 8 W eta_K k^4/(3 B f_L^2 a_s^4).
    Does NOT fix physical metric constraints.
    """
    if min(k,a_j,a_s,bf2,c_l2,m_l2,w_eta)<0 or a_j==0 or a_s==0 or bf2==0:
        raise ValueError("require positive scales and B f_L^2, others >=0")
    omega2=m_l2+c_l2*k*k/a_j**2
    lam=8*w_eta*k**4/(3*bf2*a_s**4)
    return omega2,lam,omega2+lam

def trace_static_response(k,a_j,a_s,bf2,c_l2,m_l2,w_eta,phi_minus_m):
    """Reduced static solution driven by external Phi_J - delta ln M."""
    omega2,lam,total=trace_stiffness(k,a_j,a_s,bf2,c_l2,m_l2,w_eta)
    if total<=0: raise ValueError("static inverse undefined")
    ell=3*lam/total*phi_minus_m
    q=-phi_minus_m+ell/3
    return dict(ell=ell,q=q,q_over_original=omega2/total,
                omega2=omega2,lambda_k=lam)


def trace_dynamic_derivatives(ell,ell_dot,h_j,dln_kell_dt,
                              phi_minus_m,k,a_j,a_s,bf2,c_l2,m_l2,w_eta):
    """Reduced action-derived TRACE equation in Jordan proper time.

    ddot(ell) + [3 H_J + dot(K_ell)/K_ell] dot(ell)
       + [omega0^2 + Lambda_k] ell = 3 Lambda_k (Phi_J - m).

    NOT full delta-rho_g or Jordan Einstein constraints; the metric/mapper
    are prescribed external drivers in this restricted scalar-trace sector.
    """
    omega2,lam,total=trace_stiffness(k,a_j,a_s,bf2,c_l2,m_l2,w_eta)
    if h_j<0: raise ValueError("this control assumes an expanding Jordan epoch")
    return (ell_dot, -(3*h_j+dln_kell_dt)*ell_dot-total*ell+3*lam*phi_minus_m)

def quadratic_energy_variation(k,a_j,a_s,bf2,c_l2,m_l2,w_eta,
                               ell0,ell0_dot,q0,dell,dell_dot,dq):
    """First variation at fixed coefficients of quadratic trace energy.

    E=K/2[ell_dot^2+omega0^2 ell^2]+C q^2,
    K=B f_L^2/3, C=4 W eta_K k^4/a_s^4.
    At ell0=ell0_dot=q0=0, delta E vanishes identically.
    """
    omega2,_,_=trace_stiffness(k,a_j,a_s,bf2,c_l2,m_l2,w_eta)
    return (bf2/3)*(ell0_dot*dell_dot+omega2*ell0*dell)+(
        8*w_eta*k**4/a_s**4)*q0*dq

def _growth_rhs(x,y,bg,k_hmpc,cs2_pg,smooth_pg,mu):
    """Subhorizon Jordan Newtonian-limit closure, x=ln(a).

    y=(delta_b,u_b,delta_g,u_g); u=theta/(aH/c), delta_x=-u.
    u_x = -(2+H_x/H)u - 1.5 mu*(Om_b db + Om_g dg)
          + cs2_g*(k/aH)^2 dg   (last term only geometric source).
    With cs2_g nonzero, this is a Newtonian/subhorizon approximation,
    not the exact relativistic nonadiabatic fluid system.
    """
    a=exp(x); db,ub,dg,ug=y
    source=1.5*mu*(bg.frac_b(a)*db+(0 if smooth_pg else bg.frac_pg(a)*dg))
    damp=2+bg.dln_h_dln_a(a)
    u_g_x=(0 if smooth_pg else -damp*ug-source+
           cs2_pg*(k_hmpc*bg.h/bg.hconf(a))**2*dg)
    return np.array([-ub,-damp*ub-source,0 if smooth_pg else -ug,u_g_x])

def solve_subhorizon_growth(bg,k_hmpc=0.1,z_initial=100,cs2_pg=0,
                            smooth_pg=False,mu=1,n=240):
    """After recombination, background radiation only; no baryon drag/shear.

    Equal illustrative initial seeds are NOT primordial transfer functions.
    A dust pg is mathematically identical to an equal-density CDM source.
    """
    if k_hmpc<=0 or cs2_pg<0 or z_initial<=0 or mu<=0 or n<3:
        raise ValueError("invalid solver settings")
    a0=1/(1+z_initial)
    y0=np.array([a0,-a0,0 if smooth_pg else a0,
                 0 if smooth_pg else -a0])
    xs=np.linspace(log(a0),0.,n)
    sol=solve_ivp(_growth_rhs,(xs[0],xs[-1]),y0,t_eval=xs,
                  args=(bg,k_hmpc,cs2_pg,smooth_pg,mu),method='DOP853',
                  atol=1e-11,rtol=3e-9,max_step=0.035)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(sol.message)
    return dict(a=np.exp(sol.t),delta_b=sol.y[0],u_b=sol.y[1],
                delta_pg=sol.y[2],u_pg=sol.y[3],k_hmpc=k_hmpc)

def psi_subhorizon(bg,a,k_hmpc,delta_b,delta_pg,mu=1):
    """GR-control Poisson potential, NOT SDMC-derived Jordan lensing."""
    if k_hmpc<=0: raise ValueError("k>0")
    return -1.5*(bg.hconf(a)/(k_hmpc*bg.h))**2*mu*(
        bg.frac_b(a)*delta_b+bg.frac_pg(a)*delta_pg)

def archived_reference():
    h=.6971482083084993
    ob=.022083219194622913/h**2
    og=.12299536722293603/h**2
    rad=4.18395e-5/h**2
    return Background(h,ob,og,rad,1-ob-og-rad)

def demo():
    bg=archived_reference()
    dust=solve_subhorizon_growth(bg)
    smooth=solve_subhorizon_growth(bg,smooth_pg=True)
    press=solve_subhorizon_growth(bg,cs2_pg=1e-4)
    out={
      "scope":"conditional Jordan GR/subhorizon closure, NOT action-derived no-CDM",
      "omega_b":bg.omega_b,"omega_pg":bg.omega_pg,
      "delta_b_today_dust_pg":float(dust['delta_b'][-1]),
      "delta_pg_today_dust_pg":float(dust['delta_pg'][-1]),
      "delta_b_today_smooth_pg":float(smooth['delta_b'][-1]),
      "delta_b_today_cs2_pg_1e-4":float(press['delta_b'][-1]),
      "delta_pg_today_cs2_pg_1e-4":float(press['delta_pg'][-1]),
      "limits":"no photons/neutrino Boltzmann hierarchy, no sigma dynamics, no baryon drag, no weak lensing"
    }
    print(json.dumps(out,indent=2))
    return out

if __name__=='__main__': demo()
