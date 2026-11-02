"""Post-handoff structural intrinsic spatial curvature and physical averaging.

This module is a 3-D conformal-scalar exact kinematic calculation,
NOT a microphysical initial-condition generator or global flatness proof.

Manuscript B Part III §102:
   gamma_ij=a_s² exp(2q_s(x)) delta_ij
   q_s=m+ell/3-Phi_J.

Exact intrinsic scalar spatial curvature of conformal Euclidean 3-metric:
   R3=-a_s^-2 exp(-2q)[4 Laplacian(q)+2 |grad q|²].

Structural proper volume: dV_s=a_s³ exp(3q) d³x.
For periodic geometry:
   <R3>_V= (2/a_s²)* <exp(q)|grad q|²>_coord /
                         <exp(3q)>_coord >=0.
This follows from integration by parts and includes second order.
The coordinate mean of *linear* curvature is zero:
   <delta R3_lin> =-4/a_s² <Laplacian(q)> =0.
The positive full volume mean does NOT assert K_FLRW>0,
since backreaction, slicing, Q_ij, actual 4D dynamics are omitted.

For a one-dimensional sinusoidal mode q=eps cos(kx) on a
three-dimensional torus (independent of y,z), the exact
average is:
 <R3>_V=eps² k²/a_s² *
   [I0(eps)-I2(eps)]/I0(3eps),
where I_n are modified Bessel functions.
Its linear RMS is sqrt(8)*abs(eps)*k²/a_s².

For a weak statistically isotropic periodic scalar field,
 <R3>_V = (2/a_s²) integral k² P_q(k) d³k/(2pi)³ + O(q³)
 and linear variance = (16/a_s⁴) integral k⁴ P_q(k)d³k/(2pi)³.
These calculations cannot determine the actual primordial q_k.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import iv

@dataclass(frozen=True)
class CurvatureSummary:
    epsilon: float
    mode: int
    structural_scale:float
    coordinate_mean_linear:float
    physical_mean_exact:float
    physical_mean_integral_identity:float
    linear_rms:float
    center_curvature:float
    trough_curvature:float
    nonlinear_variance:float

def sinusoid_exact(epsilon=.1, mode=1, a_s=1., samples=16384):
    if not np.isfinite(epsilon) or not np.isfinite(a_s) or a_s<=0 or \
       int(mode)!=mode or mode<1 or samples<256:
        raise ValueError("finite epsilon, positive structural scale, integer mode, >=256 samples")
    grid=2*np.pi*np.arange(samples)/samples
    q=epsilon*np.cos(mode*grid)
    gradient=-epsilon*mode*np.sin(mode*grid)
    laplacian=-epsilon*mode*mode*np.cos(mode*grid)
    lin_curv=-4*laplacian/a_s**2
    full_curv=-np.exp(-2*q)*(4*laplacian+2*gradient**2)/a_s**2
    weights=np.exp(3*q)
    mean=float(np.sum(weights*full_curv)/np.sum(weights))
    gradient_identity=float(2*np.mean(np.exp(q)*gradient**2)/
                            (a_s**2*np.mean(weights)))
    variance=float(np.sum(weights*(full_curv-mean)**2)/np.sum(weights))
    center=float(4*epsilon*mode**2*np.exp(-2*epsilon)/a_s**2)
    trough=float(-4*epsilon*mode**2*np.exp(2*epsilon)/a_s**2)
    return CurvatureSummary(float(epsilon),int(mode),float(a_s),
      float(np.mean(lin_curv)),mean,gradient_identity,
      float(np.sqrt(np.mean(lin_curv**2))),
      center,trough,variance)

def sinusoid_bessel_mean(epsilon=.1,mode=1,a_s=1.):
    if a_s<=0 or mode<=0:return ValueError("positive scale and mode needed")
    return float(epsilon**2*mode**2/a_s**2*
            (iv(0,epsilon)-iv(2,epsilon))/iv(0,3*epsilon))

def periodic_fourier_metrics(amplitudes,modes,phases=None,a_s=1.,
                              samples=16384):
    """Many locally correlated Fourier modes, 3D conformal metric q(x).

    Mathematical model is 3D yet q varies along one periodic axis.
    q(x)=sum eps_i cos(k_i x+phase_i).
    Only spatial scalar strain; no traceless Q_ij.
    """
    amps=np.atleast_1d(np.asarray(amplitudes,dtype=float))
    ks=np.atleast_1d(np.asarray(modes,dtype=int))
    if phases is None:
        phases=np.zeros_like(amps)
    phases=np.atleast_1d(np.asarray(phases,dtype=float))
    if amps.ndim!=1 or amps.size==0 or amps.shape!=ks.shape or \
       amps.shape!=phases.shape or np.any(ks<1) or \
       a_s<=0 or samples<256 or \
       not np.isfinite(amps).all() or not np.isfinite(phases).all():
        raise ValueError("valid finite periodic modes")
    x=np.arange(samples)*2*np.pi/samples
    theta=ks[:,None]*x[None,:]+phases[:,None]
    q=np.sum(amps[:,None]*np.cos(theta),axis=0)
    grad=-np.sum((amps*ks)[:,None]*np.sin(theta),axis=0)
    lap=-np.sum((amps*ks**2)[:,None]*np.cos(theta),axis=0)
    curv=-np.exp(-2*q)*(4*lap+2*grad**2)/a_s**2
    vweights=np.exp(3*q)
    mean=float(np.sum(curv*vweights)/np.sum(vweights))
    mean_identity=float(2*np.mean(np.exp(q)*grad**2)/
                        (a_s*a_s*np.mean(vweights)))
    # Weak-amplitude spectral identities for *distinct* Fourier harmonics
    if len(np.unique(ks))!=len(ks):
        weak_predicted_mean=None
        weak_lin_variance=None
    else:
        weak_predicted_mean=float(np.sum((amps*ks)**2)/a_s**2)
        weak_lin_variance=float(8*np.sum((amps*ks**2)**2)/a_s**4)
    return dict(n_modes=int(amps.size),
       coordinate_mean_laplacian=float(np.mean(lap)),
       exact_volume_mean_curvature=mean,
       integration_by_parts_mean=mean_identity,
       weak_mode_power_mean=weak_predicted_mean,
       weak_linear_curvature_variance=weak_lin_variance,
       physical_volume_weighted_variance=float(
           np.sum((curv-mean)**2*vweights)/np.sum(vweights)),
       free_mode_initial_spectrum_not_derived=True,
       homogeneous_FLRW_K_not_inferred=True)

def example():
    xs=[sinusoid_exact(x) for x in (1e-3,.01,.05,.1,.2)]
    multi=periodic_fourier_metrics(
        amplitudes=[.03,.01,.007],
        modes=[1,3,7],phases=[.2,1.1,2.3])
    return dict(scope="kinematic post-handoff 3D structural metric, not dynamical Planck-to-radiation solution",
      examples=[vars(v) for v in xs],
      multimode=multi,
      missing=["actual primordial q_k from full prethermal field equations",
               "traceless strain Q_ij and structural spatial gauge",
               "nonflat background K and tensor ADM constraints",
               "physical mode occupation and Jordan frame pressure",
               "genuine thermal transfer and causal matching"])
if __name__=="__main__":
    import json
    print(json.dumps(example(),indent=2))
