"""SDMC Part III — canonical occupied-mode production from a finite K(t) handoff.

This is an explicitly CONDITIONAL quantum-oscillator diagnostic of the
already displayed structural trace action at W=0, not a derivation from
the exact quantum Planck boundary or the full gravitating EFT.

The classical mode A of Manuscript B Eq 106.3 has
  S_A=∫dt (a_J³ K(t)/4)[dot A²-omega_L²(t) A²].
Define b(t)=a_J³ K(t)/2, v=sqrt(b)*A, Gamma=dot b/b.
After removing a total derivative, its canonical action is
  S_v=∫dt 1/2 [dot v²-Omega_v² v²],
  Omega_v²=omega_L²-0.5 dot Gamma-0.25 Gamma².

At EXACT classical A=dot A=0 with no source, all classical solutions
remain zero. Nonadiabatic change in canonical Omega_v, if the field
admits quantization, can nevertheless mix positive and negative
frequency vacuum modes and generate |beta_k|²>0 occupation.

To avoid claiming an uncomputed expanding prethermal spacetime,
we quantify only a locally frozen-Jordan expansion benchmark:
  H_J=0, a_J=1, W=0, omega_L=sqrt(m_L²+c_L² k²)>0,
  K(t)=K_initial exp[ln(K_final/K_initial) sigmoid(t/tau)].
The initial and final positive-frequency bases are well-defined
because dot K/K->0 away from the finite transition.

Canonical Wronskian normalization:
  v_i=(2omega_L)^(-1/2),
  dot v_i=-i(omega_L/2)^(1/2).
At t_f:
  alpha=sqrt(omega/2)v+i dot v/sqrt(2omega),
  beta =sqrt(omega/2)v-i dot v/sqrt(2omega),
  |alpha|²-|beta|²=1, n_k=|beta|².

A fast but stable coefficient transition can produce a nonzero
canonical occupation; a slow adiabatic ramp yields tiny occupation.
This gives a possible EFT *source of initial occupation* even with
no explicit material-to-A force, but not a primordial cosmological
amplitude, density spectrum or CDM replacement.

CRITICAL LIMITS:
 * B_initial=0 makes the isolated kinetic rank degenerate and this
   canonical construction invalid. A positive B_floor is required.
 * H=0 is a local benchmark only; full SDMC H_J, p_v, structural
   clock, EFT cutoff, scalar/material constraints and mode spectrum
   are not included.
 * Quantum vacuum and choice of initial basis are assumptions.
 * Number density rho_g ~ int d^3k omega_k n_k/a_J³ still requires
   a physical cutoff, state, mode measure and gravitational stress.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class CanonicalHandoff:
    B_initial:float=.2
    B_final:float=1.
    tau:float=.4
    m_L:float=1.
    c2:float=.1
    k:float=.3
    t_window:float=12.

    def __post_init__(self):
        vals=np.asarray(list(vars(self).values()),dtype=float)
        if not np.isfinite(vals).all():raise ValueError("finite coefficients")
        if min(self.B_initial,self.B_final,self.tau,self.m_L,
               self.t_window)<=0:
            raise ValueError("strictly positive canonical kinetic floor, mass, width")
        if self.c2<0 or self.k<0:raise ValueError("nonnegative gradients")

def omega_out(e):
    return np.sqrt(e.m_L**2+e.c2*e.k**2)

def K_history(e,t):
    z=np.clip(t/e.tau,-100,100)
    u=1/(1+np.exp(-z))
    Delta=np.log(e.B_final/e.B_initial)
    K=e.B_initial*np.exp(Delta*u)
    rate=Delta*u*(1-u)/e.tau
    ratedot=Delta*u*(1-u)*(1-2*u)/e.tau**2
    return dict(K=K,rate=rate,ratedot=ratedot,u=u)

def canonical_frequency_squared(e,t):
    co=K_history(e,t)
    return omega_out(e)**2-.5*co["ratedot"]-.25*co["rate"]**2

def canonical_mode(e,samples=451):
    """Canonical complex in-vacuum mode through finite classical handoff."""
    if samples<20:raise ValueError("minimum resolution")
    t0=-e.t_window*e.tau;tf=e.t_window*e.tau
    times=np.linspace(t0,tf,samples)
    omega=omega_out(e)
    v0=1/np.sqrt(2*omega)
    vdot0=-1j*np.sqrt(omega/2)
    def rhs(t,y):
        return [y[1],-canonical_frequency_squared(e,t)*y[0]]
    sol=solve_ivp(rhs,(t0,tf),[v0,vdot0],
      t_eval=times,method="DOP853",rtol=3e-11,atol=3e-12,
      max_step=min(.04,e.tau/35))
    if not sol.success:raise RuntimeError(sol.message)
    v,vdot=sol.y
    alpha=np.sqrt(omega/2)*v[-1]+1j*vdot[-1]/np.sqrt(2*omega)
    beta=np.sqrt(omega/2)*v[-1]-1j*vdot[-1]/np.sqrt(2*omega)
    n=float(abs(beta)**2)
    Wronskian=1j*(np.conjugate(v)*vdot-v*np.conjugate(vdot))
    final_E=.5*(abs(vdot[-1])**2+omega**2*abs(v[-1])**2)
    occupation_energy=omega*n
    min_omega_sq=float(min(
       canonical_frequency_squared(e,float(t)) for t in times))
    return dict(times=times,v=v,vdot=vdot,
      n=n,alpha=alpha,beta=beta,
      alpha_beta_norm=float(abs(alpha)**2-abs(beta)**2),
      max_Wronskian_residual=float(max(abs(Wronskian-1))),
      final_energy=float(final_E),
      vacuum_energy=.5*omega,
      final_excess_energy=float(final_E-.5*omega),
      occupation_energy=float(occupation_energy),
      final_excess_energy_error=float(
                       (final_E-.5*omega)-occupation_energy),
      min_canonical_omega_squared=min_omega_sq,
      change_ratio=e.B_final/e.B_initial,
      initial_rate=float(K_history(e,t0)["rate"]),
      final_rate=float(K_history(e,tf)["rate"]),
      assumed_initial_vacuum=True)

def fast_slow_comparison():
    fast=CanonicalHandoff(tau=.4)
    medium=CanonicalHandoff(tau=.8)
    slow=CanonicalHandoff(tau=1.4)
    unchanged=CanonicalHandoff(B_initial=.2,B_final=.2,tau=.4)
    rows={}
    for tag,e in [("fast",fast),("medium",medium),
                  ("slow",slow),("unchanged",unchanged)]:
        x=canonical_mode(e)
        rows[tag]=dict(
          tau=e.tau,B_initial=e.B_initial,B_final=e.B_final,
          n_k=x["n"],
          nonadiabatic_amplification=bool(x["n"]>1e-6),
          positive_omega_squared=bool(x["min_canonical_omega_squared"]>0),
          min_canonical_omega_squared=x["min_canonical_omega_squared"],
          Wronskian_residual=x["max_Wronskian_residual"],
          norm_check=x["alpha_beta_norm"],
          vacuum_subtracted_energy=x["final_excess_energy"],
          energy_vs_nomega_error=x["final_excess_energy_error"])
    return dict(
      scope="finite positive kinetic-ramp Bogoliubov diagnostic in frozen Jordan local limit; not a primordial SDMC abundance",
      cases=rows,
      classically_exact_zero_field_remains_zero=True,
      physically_fixed_spectrum=False,
      source_curvature_selector_W=0,
      missing=[
        "positive-rank canonical birth and initial quantum state near S_A≈1",
        "actual post-Planck background a_J(t_c), K(sigma(t_c)), mapper and cutoffs",
        "3D tensor/trace gauge constraints, strong coupling, ghost/gradient checks",
        "k-dependent occupation across the physical SDMC spectrum",
        "renormalized rho_g, pressure, Jg, clustering and Jordan Weyl response",
        "matter-radiation equality, CMB, weak-lensing and late300 no-CDM likelihoods"
      ])
if __name__=="__main__":
    print(json.dumps(fast_slow_comparison(),indent=2))
