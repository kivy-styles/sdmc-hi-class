"""Quadratic SDMC supporting-space source: mode coupling and Jordan averages.

This is the *reduced* scalar-trace sector of Manuscript B's ADM candidate,
not a covariant, fully constrained perturbation theory. It derives a positive
quadratic proxy for structural energy and its inhomogeneous mode coupling.

Important:
 - rho_g(action) here is computed at fixed ADM/Jordan lapse coefficients.
 - spatial averaging does NOT establish physical pressure, momentum, or Weyl.
 - S, structural time, and a_J remain separate; a_s=S/S0.
 - output units depend on the action coefficients; not matched to omega_c.
"""
from dataclasses import dataclass
from typing import Optional
import json
import numpy as np

@dataclass(frozen=True)
class TraceCoefficients:
    """Constant coefficients at one Jordan epoch (symbolic physical units)."""
    bf2: float = 1.0
    c_l2: float = 1.0
    m_l2: float = 1.0
    w_eta: float = 1.0
    a_j: float = 1.0
    a_s: float = 1.0
    length: float = 2*np.pi

    def __post_init__(self):
        if self.bf2<=0 or self.a_j<=0 or self.a_s<=0 or self.length<=0:
            raise ValueError("bf2, Jordan/structural scales and length must be positive")
        if min(self.c_l2,self.m_l2,self.w_eta)<0:
            raise ValueError("positive quadratic energy branch only")

    @property
    def K(self): return self.bf2/3
    @property
    def C(self): return 4*self.w_eta/self.a_s**4

def _array(x,n=None):
    y=np.asarray(x,dtype=float)
    if y.ndim!=1 or y.size<16 or not np.all(np.isfinite(y)):
        raise ValueError("need finite 1d periodic field with at least 16 grid points")
    if n is not None and len(y)!=n: raise ValueError("field length mismatch")
    return y

def derivative(field,length,order):
    """Spectral derivative on a periodic comoving x domain."""
    x=_array(field)
    if order not in (1,2,4): raise ValueError("order must be 1,2,4")
    wave=2*np.pi*np.fft.fftfreq(x.size,d=length/x.size)
    val=np.fft.ifft((1j*wave)**order*np.fft.fft(x))
    return val.real

def gaussian_coarse_grain(field,length,width):
    """Gaussian low-pass in a periodic Jordan comoving coordinate."""
    x=_array(field)
    if width<0: raise ValueError("width >= 0")
    if width==0: return x.copy()
    wave=2*np.pi*np.fft.fftfreq(x.size,d=length/x.size)
    filtered=np.fft.ifft(
        np.fft.fft(x)*np.exp(-.5*(wave*width)**2)).real
    return filtered

def quadratic_energy_terms(ell,coeff,ell_dot=None,m=None,phi=None):
    """Jordan-volume local energy proxy in the fixed-coefficient trace sector.

    q=ell/3 + m - Phi_J.
    E_kin=K/2 dot(ell)^2
    E_mass=K*m_L^2/2 ell^2
    E_grad=K*c_L^2/(2a_J^2)(partial_x ell)^2
    E_curv=4 W eta_K/a_s^4 (partial_x^2 q)^2.

    Scalar q geometry is a 1d periodic profile embedded in 3 spatial
    dimensions. The Riemann-squared contraction gives 8(q_xx)^2/a_s^4.
    """
    ell=_array(ell)
    n=ell.size
    zero=np.zeros(n)
    vel=_array(ell_dot,n) if ell_dot is not None else zero
    mm=_array(m,n) if m is not None else zero
    pp=_array(phi,n) if phi is not None else zero
    q=ell/3+mm-pp
    return {
        "kinetic":.5*coeff.K*vel**2,
        "mass":.5*coeff.K*coeff.m_l2*ell**2,
        "gradient":.5*coeff.K*coeff.c_l2/coeff.a_j**2*
                    derivative(ell,coeff.length,1)**2,
        "curvature":coeff.C*derivative(q,coeff.length,2)**2,
    }

def quadratic_energy_density(ell,coeff,ell_dot=None,m=None,phi=None):
    terms=quadratic_energy_terms(ell,coeff,ell_dot,m,phi)
    return sum(terms.values())

def first_variation_energy(ell,dell,coeff,
                           ell_dot=None,dell_dot=None,
                           m=None,dm=None,phi=None,dphi=None):
    """First functional variation at FIXED coefficients and a_j,a_s.

    Vanishes about a perfectly locked-flat unexcited background.
    Nonzero for a pre-existing inhomogeneous structural state.
    Neither a physical gauge invariant delta-rho nor complete T_mu_nu.
    """
    ell=_array(ell); n=ell.size
    de=_array(dell,n); z=np.zeros(n)
    vel=_array(ell_dot,n) if ell_dot is not None else z
    dv=_array(dell_dot,n) if dell_dot is not None else z
    mm=_array(m,n) if m is not None else z
    dmm=_array(dm,n) if dm is not None else z
    pp=_array(phi,n) if phi is not None else z
    dpp=_array(dphi,n) if dphi is not None else z
    q=ell/3+mm-pp
    dq=de/3+dmm-dpp
    return (coeff.K*(vel*dv+coeff.m_l2*ell*de+
        coeff.c_l2/coeff.a_j**2*
        derivative(ell,coeff.length,1)*derivative(de,coeff.length,1))
        +2*coeff.C*derivative(q,coeff.length,2)*
                    derivative(dq,coeff.length,2))

def averaged_energy(ell,coeff,ell_dot=None,m=None,phi=None):
    return float(quadratic_energy_density(ell,coeff,ell_dot,m,phi).mean())

def coarse_grained_contrast(ell,coeff,width,ell_dot=None,m=None,phi=None):
    """Contrast relative to nonzero action energy mean, NOT the CDM delta."""
    rho=quadratic_energy_density(ell,coeff,ell_dot,m,phi)
    avg=float(np.mean(rho))
    if avg<=0: raise ValueError("zero physical energy: no contrast exists")
    return gaussian_coarse_grain(rho,coeff.length,width)/avg-1

def energy_mode(ell,ell_dot,u,k,coeff):
    """Single scalar-trace mode's quadratic energy coefficient.

    u = delta ln M - Phi_J; all values at one Jordan epoch.
    The E is not an independently conserved global geometric energy.
    """
    k=float(k)
    omega2=coeff.m_l2+coeff.c_l2*k*k/coeff.a_j**2
    Ck=coeff.C*k**4
    q=ell/3+u
    return .5*coeff.K*(ell_dot**2+omega2*ell**2)+Ck*q*q

def mode_energy_rate(ell,ell_dot,u,u_dot,k,coeff,H_j,H_s):
    """On-shell dE/dt_c for fixed bf2,c_l2,m_l2,w_eta.

    Uses the action-derived reduced trace EOM (including its 3H friction).
    Derivation:
    dE/dt=-3H_J K ldot^2
          -H_J K (c_l^2 k^2/a_J^2)ell^2
          -4H_s C_k q^2
          +2 C_k q u_dot
    with H_s=d ln a_s/dt_c. Last term is external metric/mapper work.
    """
    k=float(k); q=ell/3+u; Ck=coeff.C*k**4
    return (-3*H_j*coeff.K*ell_dot**2
            -H_j*coeff.K*coeff.c_l2*k*k/coeff.a_j**2*ell**2
            -4*H_s*Ck*q*q+2*Ck*q*u_dot)


def restricted_isotropic_pressure(ell,coeff,ell_dot=None,m=None,phi=None):
    """Restricted homogeneous Jordan-a variation at FIXED structural S and fields.

    The trace action has L=Ekin-Emass-Egrad-Ecurv, with
    Egrad proportional to a_J^-2 and Ecurv proportional to a_s^-4.
    When varying a_J at fixed S (hence fixed a_s), the spatial scale
    derivative of a_J^3 L gives:
        P_red=Ekin-Emass-Egrad/3-Ecurv.
    This is NOT the complete pressure from varying all coupled ADM fields
    or from an averaged inhomogeneous cosmological spacetime.
    """
    t=quadratic_energy_terms(ell,coeff,ell_dot,m,phi)
    return t["kinetic"]-t["mass"]-t["gradient"]/3-t["curvature"]

def restricted_isotropic_pressure_mean(ell,coeff,ell_dot=None,m=None,phi=None):
    return float(restricted_isotropic_pressure(ell,coeff,ell_dot,m,phi).mean())

def integrated_static_kernel(k,coeff):
    """Eliminate the quadratic structural trace at FIXED metric/mapper driver.

    E_static=K*omega0^2/2*ell^2+C_k*(u+ell/3)^2
    => E_eff(u)=Keff*u^2 where
       Keff=C_k*omega0^2/(omega0^2+2*C_k/(9*K)).
    No Poisson source or gravitational mu is inferred from Keff alone.
    """
    k=float(k)
    if k<0: raise ValueError("k>=0")
    omega2=coeff.m_l2+coeff.c_l2*k*k/coeff.a_j**2
    ck=coeff.C*k**4
    stiff=2*ck/(9*coeff.K)
    if omega2+stiff<=0: return 0.
    return ck*omega2/(omega2+stiff)

def integrated_static_energy(u,k,coeff):
    return integrated_static_kernel(k,coeff)*u*u

def harmonic_demo():
    """Deterministic 1d short-wave / long-envelope example only."""
    n=1024; length=2*np.pi
    x=np.arange(n)*length/n
    short,large,epsilon=32,2,.1
    base=np.cos(short*x)
    envelope=1+epsilon*np.cos(large*x)
    coeff=TraceCoefficients(length=length, bf2=3.,c_l2=1.,
                             m_l2=.05,w_eta=.01,a_j=1.,a_s=1.)
    zero=np.zeros(n)
    rho0=averaged_energy(base,coeff)
    rho=quadratic_energy_density(envelope*base,coeff)
    width=.16
    contrast=gaussian_coarse_grain(rho,length,width)/rho.mean()-1
    signal=float(2*np.mean(contrast*np.cos(large*x)))
    theoretical=2*epsilon*np.exp(-.5*(large*width)**2)
    return {"scope":"quadratic structural trace energy, 1D illustrative support",
      "mean_energy_unmodulated_arbitrary_units":rho0,
      "mean_energy_modulated_arbitrary_units":float(rho.mean()),
      "long_density_contrast_cosine_amplitude":signal,
      "leading_separation_approximation":float(theoretical),
      "long_wavenumber":large,"short_wavenumber":short,
      "not_tested":["full ADM lapse/shift constraints","Phi_J and Psi_J response",
                    "observable gravitational lensing","physical isotropic 3D ensemble",
                    "independent rho_g normalization"]}

if __name__=="__main__":
    print(json.dumps(harmonic_demo(),indent=2))
