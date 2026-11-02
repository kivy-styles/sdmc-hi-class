"""SDMC Part III: conserved material redistribution -> reduced structural spectrum.

This is a minimal *diagnostic*, not a complete no-CDM cosmological solver.

Material inventory is conserved:
  Mu=M_remaining+M_inferred,
  delta(k)=[M_remaining * W_core(k) e^(-ik x_core)
             + sum_i m_i W_i(k) e^(-ik x_i)]/Mu
  for nonzero periodic k.
This is the matter density contrast mode of a finite reference box;
M_inferred denotes actual material distributed into the expanded region,
NOT newly created dark mass.

For independent, uniform positions and normalized weights w_i:
  E[delta(k)]=0,
  E|delta(k)|²=sum_i w_i² =1/N_eff.

For iid positions with distribution p(x)=(1+eps cos(kx))/L:
  E delta(k)=eps/2 (real phase convention),
  E|delta(k)|²=|mu|²+(1-|mu|²)/N_eff.

The reduced supporting curvature action for mode u=m-Phi_J:
  E(l,u)=K_l omega0² |l|²/2+C |u+l/3|²,
  K_l=B f_L²/3; C=4W eta_K k⁴/a_s⁴,
  Lambda=2C/(9K_l).
Its exact quasi-static minimum is
  T=omega0²/(omega0²+Lambda)
  ell=-3(1-T)u; q=T u;
  E_min=C T |u|²;
  E_curvature=C T² |u|²;
  E_relative=C T(1-T)|u|².

Poisson u=A_m delta/k², with A_m=4pi G_eff a_J² rho_m.
This approximation is meaningful only inside the horizon, for a
weakly perturbed nonrelativistic matter field. It MUST NOT be used
at the exact Planck boundary or as the full relativistic radiation
or Jordan modified-gravity constraint.

Then E_min=4W eta_K/a_s⁴ * A_m²*T*|delta(k)|²:
the apparent k⁴ in the curvature operator cancels the k^-4 of
the Newtonian potential squared, leaving T(k) and the matter power.
This energy is NOT yet a physical global Friedmann rho_g.

For conserved nonrelativistic rho_m~a_J^-3 and fixed G_eff:
 A_m~a_J^-1. The diagnostic ratio E_min/rho_v obeys
 Gamma_stat~ W a_s^-2 a_J^-2 T(k) P_delta(k).
For p_v=dln(a_s)/dln(a_J), frozen W,T:
 dlnGamma/dlna_J=-2*p_v-2+dln(P_delta)/dlna_J.
To match the manuscript's dust condition 2p_v-3 would demand
 dln(P_delta)/dlna_J=4p_v-1 (radiation 7, matter 5).
This is a CONDITIONAL stress-spectrum diagnostic, NOT a
physical source identity rho_g=E_min or a no-go theorem.
When W->0, E_min and the source of ell vanish on this forced
quasi-static branch; independent *homogeneous* modes may survive.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np

@dataclass(frozen=True)
class SelectorParameters:
    B: float=.5
    W: float=1.
    eta: float=.01
    fL: float=1.
    m2: float=1.
    c2: float=.1
    aJ: float=1.
    aS: float=1.

    def __post_init__(self):
        vals=np.asarray(list(vars(self).values()),dtype=float)
        if not np.isfinite(vals).all():
            raise ValueError("Finite parameters required")
        if self.B<=0 or self.B>1 or min(self.eta,self.W,self.fL,
                     self.m2,self.c2)<0 or self.aJ<=0 or self.aS<=0:
            raise ValueError("Regular stable B>0, nonnegative sectors and positive rulers")
        if self.fL==0 or self.m2+self.c2==0:
            raise ValueError("Nonzero kinetic and oscillator restoring normalization")

def response(k,p:SelectorParameters):
    """Static structural response of one positive k Fourier mode."""
    if k<=0 or not np.isfinite(k):
        raise ValueError("Only finite k>0, not global FLRW K")
    K=p.B*p.fL**2/3
    C=4*p.W*p.eta*k**4/p.aS**4
    omega02=p.m2+p.c2*k*k/p.aJ**2
    lam=2*C/(9*K)
    transmission=omega02/(omega02+lam)
    return dict(K=K,C=C,omega02=omega02,Lambda=lam,
                transmission=transmission,relative_fraction=1-transmission)

def source_energy(k,u,p:SelectorParameters):
    """Evaluate exactly the on-shell stationary quadratic mode energy."""
    r=response(k,p)
    ell=-3*r["relative_fraction"]*u
    q=u+ell/3
    Erel=.5*r["K"]*r["omega02"]*abs(ell)**2
    Ecurv=r["C"]*abs(q)**2
    total=Erel+Ecurv
    theory=r["C"]*r["transmission"]*abs(u)**2
    return dict(ell=ell,q=q,relative=Erel,curvature=Ecurv,
                total=total,expected=theory,
                on_shell_residual=total-theory)

def off_shell_energy(k,u,ell,p):
    r=response(k,p)
    return (.5*r["K"]*r["omega02"]*abs(ell)**2
             +r["C"]*abs(u+ell/3)**2)

def density_mode(positions,k,masses=None,length=2*np.pi,
                 core_mass=0.,core_position=0.,core_width=0.):
    """Nonzero periodic x-direction density contrast mode.

    A localized remaining core may coexist with redistributed particles;
    a Gaussian core window exp[-(k sigma)^2/2] is optional.
    Materials always sum to Mu; the core is not extra mass.
    """
    x=np.asarray(positions,dtype=float)
    if x.ndim!=1 or not np.isfinite(x).all():
        raise ValueError("finite vector of material positions required")
    if masses is None:
        masses=np.ones(len(x))
    mass=np.asarray(masses,dtype=float)
    if mass.shape!=x.shape or not np.isfinite(mass).all() or np.any(mass<0):
        raise ValueError("nonnegative masses of identical shape")
    if core_mass<0 or core_width<0 or length<=0 or k<=0:
        raise ValueError("physical nonnegative masses/scales and positive k,L")
    Mu=float(mass.sum()+core_mass)
    if Mu<=0: raise ValueError("positive conserved Mu needed")
    kphysical=2*np.pi*k/length
    distributed=np.sum(mass*np.exp(-1j*kphysical*x))
    core=core_mass*np.exp(-1j*kphysical*core_position)
    core*=np.exp(-.5*(kphysical*core_width)**2)
    return (distributed+core)/Mu

def effective_count(masses):
    masses=np.asarray(masses,dtype=float)
    if masses.ndim!=1 or np.any(masses<0) or masses.sum()<=0:
        raise ValueError("positive weighted distribution required")
    return float(masses.sum()**2/np.dot(masses,masses))

def iid_mode_moments(masses,mean_single):
    """Exact first and second Fourier-moment for iid positions."""
    m=np.asarray(masses,dtype=float)
    Neff=effective_count(m)
    mu=complex(mean_single)
    return dict(mean_real=float(mu.real),mean_imag=float(mu.imag),
           power=float(abs(mu)**2+(1-abs(mu)**2)/Neff),
           variance=float((1-abs(mu)**2)/Neff),
           effective_count=Neff)

def correlated_positions(count,epsilon,rng):
    """Exact rejection-free inverse CDF for p(x)=(1+eps cos x)/(2pi).

    CDF=(x+eps sin x)/(2pi), strictly monotonic for |eps|<1.
    Positions are independent even though the one-body
    spatial distribution is nonuniform; the phase-coherent
    expectation at k=1 is eps/2.
    """
    if count<1 or not (-1<epsilon<1):
        raise ValueError("count>=1 and abs(epsilon)<1")
    z=2*np.pi*rng.random(count)
    x=z.copy()
    for _ in range(12):
        x-=(x+epsilon*np.sin(x)-z)/(1+epsilon*np.cos(x))
    return x

def poisson_source(k,delta,A_m):
    """Subhorizon Jordan Poisson approximation, neglecting mapper m_k.

    The sign convention u=-Phi_J and -k² Phi_J=A_m delta
    makes u=A_m delta/k²; sign does not alter energy.
    """
    if k<=0 or A_m<0:raise ValueError("positive k and A_m>=0")
    return A_m*delta/k**2

def spectrum_response(deltas,k_values,p,A_m):
    """Energy sum and modes for positive k; no independent abundance."""
    ks=np.asarray(k_values,dtype=float)
    modes=np.asarray(deltas)
    if ks.ndim!=1 or modes.shape!=ks.shape or len(ks)==0:
        raise ValueError("matched nonempty array of Fourier modes")
    rows=[]
    for k,d in zip(ks,modes):
        u=poisson_source(float(k),d,A_m)
        s=source_energy(float(k),u,p)
        r=response(float(k),p)
        rows.append(dict(k=float(k),
                power=float(abs(d)**2),
                transmission=r["transmission"],
                induced_relative=float(abs(s["ell"])),
                energy=float(s["total"]),
                predicted_cancelled_k4_energy=float(
                   4*p.W*p.eta/p.aS**4*A_m**2*
                   r["transmission"]*abs(d)**2)))
    return dict(modes=rows,E_sum=float(sum(x["energy"] for x in rows)),
                observational_gravity_proven=False)

def activation_exponents(p_v):
    """Conditional Gamma_E=E_min/rho_v, W,T and P_delta histories.

    d ln Gamma_E /dln a=-2p_v-2 + dln W + dln T + dln P_delta
    and desired dust slope is 2p_v-3.
    """
    if p_v<=0:raise ValueError("positive structural tracking exponent")
    return dict(p_v=float(p_v),
          dust_required_gamma_slope=float(2*p_v-3),
          baseline_gamma_slope=float(-2*p_v-2),
          required_power_slope_if_W_T_constant=float(4*p_v-1),
          required_delta_amplitude_slope=float(2*p_v-.5))

def diagnostic(seed=20261102):
    rng=np.random.default_rng(seed)
    count=4096
    power={}
    for name,epsilon in [("uniform",0.),("density_correlated",.12)]:
        amp=[]
        for _ in range(128):
            positions=correlated_positions(count,epsilon,rng)
            amp.append(density_mode(positions,k=1))
        amp=np.asarray(amp)
        expected=iid_mode_moments(np.ones(count),epsilon/2)
        power[name]=dict(N=count,epsilon=epsilon,
            empirical_mean_real=float(amp.mean().real),
            empirical_mean_imag=float(amp.mean().imag),
            predicted_mean=epsilon/2,
            empirical_power=float(np.mean(abs(amp)**2)),
            predicted_power=expected["power"])
    pars=SelectorParameters()
    shape=spectrum_response([.05,.015,.004],[1,2,4],pars,A_m=.01)
    zero=spectrum_response([.05,.015,.004],[1,2,4],
          SelectorParameters(W=0.),A_m=.01)
    return dict(scope="conditional matter spectrum -> quasi-static structural energy; not actual FLRW rho_g",
         random_modes=power,
         spectrum=shape,
         selector_off=zero,
         required_activation={
            "radiation":activation_exponents(2.),
            "matter":activation_exponents(1.5),
            "coasting":activation_exponents(1.)},
         omitted=["full Jordan F/G2/G3 Poisson constraints",
                  "stochastic initial conditions from Planck-to-radiation dynamics",
                  "matter-structure transfer and photon-baryon pressure",
                  "metric potential feedback and observable effective gravity",
                  "physical global rho_g,P_g,J_g and amplitude normalization"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
