"""Exact conformally flat structural Riemann² and the first cubic metric vertex.

For the scalar trace of SDMC's supporting metric,
 gamma_ij=a_s² e^{2q} delta_ij, q=ell/3+zeta,
the FULL 3D conformal Riemann squared invariant is

 Riem^2=4 a_s^-4 e^-4q [
   H_ij H_ij + (Delta q)^2
  +2((Delta q)|grad q|²-H_ij q_i q_j)
  +|grad q|^4 ].

The spatial curvature ADM action on a Jordan slice is
 L_K= -N a_J³ e^{3zeta} (Weta/2) Riem^2.

This code keeps exact nonlinear 3D conformal curvature but only
the SCALAR TRACE and constant-background Weta, F/G2/G3 fields.
It derives the cubic *Jordan spatial-curvature metric source*
that a purely quadratic trace calculation about ell=0 misses.
Lapse/shift and other supporting-space vertices of the FULL
covariant action are not computed by this module.
"""
from __future__ import annotations
import json
import numpy as np
from sdmc.pg_quadratic_source_3d import geometry_derivatives
from sdmc.pg_energy_exchange_3d import _cube
from sdmc.pg_quadratic_source import TraceCoefficients

def _derivatives(q,L):
    grad,H,lap=geometry_derivatives(q,L)
    w=sum(v*v for v in grad)
    T=sum(H[i][j]*H[i][j] for i in range(3)
                    for j in range(3))+lap*lap
    Hvv=sum(H[i][j]*grad[i]*grad[j]
                    for i in range(3) for j in range(3))
    U=2*(lap*w-Hvv)
    return dict(grad=grad,H=H,lap=lap,w=w,T=T,U=U,V=w*w)

def exact_conformal_riemann2(q,coefficient):
    """Nonlinear R_ijkl R^ijkl of gamma_ij=a_s² e^{2q}delta_ij."""
    q=_cube(q)
    g=_derivatives(q,coefficient.length)
    I=4*(g["T"]+g["U"]+g["V"])
    value=np.exp(-4*q)*I/coefficient.a_s**4
    if np.min(value)<-1e-8:
        raise ArithmeticError("nonnegative 3D Riemann-squared invariant violated")
    return value

def exact_curvature_action(ell,zeta,coefficient):
    """Exact conformal spatial Riemann² Jordan action / a_J³.

    N=1, no traceless L, metric h_ij=a_J² e^{2zeta}delta_ij.
    """
    e=_cube(ell)
    z=_cube(zeta,len(e))
    q=e/3+z
    R2=exact_conformal_riemann2(q,coefficient)
    return float(np.mean(-.5*coefficient.w_eta*np.exp(3*z)*R2))

def quadratic_curvature_action(ell,zeta,coefficient):
    """Second-order curvature action, only T(q)."""
    e=_cube(ell);z=_cube(zeta,len(e))
    q=e/3+z
    g=_derivatives(q,coefficient.length)
    D=2*coefficient.w_eta/coefficient.a_s**4
    return -D*float(g["T"].mean())

def through_cubic_curvature_action(ell,zeta,coefficient):
    """Full curvature spatial action to cubic order in ell,zeta.

    Since exp(3zeta-4q)=exp(-4ell/3-zeta), cubic action is
       -D [T(q)+U(q)-(4ell/3+zeta)T(q)].
    """
    e=_cube(ell);z=_cube(zeta,len(e))
    core=e/3; q=core+z
    g=_derivatives(q,coefficient.length)
    D=2*coefficient.w_eta/coefficient.a_s**4
    return -D*float(np.mean(g["T"]+g["U"]-
                                 (4*core+z)*g["T"]))

def cubic_spatial_metric_source(ell,coefficient):
    """Analytic O(ell²) metric source δL_K^(3)/δzeta at zeta=0.

    R_zeta^(2)/D =
      T+8 ∂ij(e H_ij)+8 Δ(e Δe)
       -2 Δ|grad e|²
       +2 ∂ij(e_i e_j)
       +4 ∂i[(Δe)e_i-H_ij e_j],

    e=ell/3. This follows by integration by parts on a
    periodic 3D volume. Its mean is D<T>, as all other
    terms are total spatial divergences.

    It is a spatial-metric tadpole of the curvature operator,
    NOT the complete physical Jordan pressure/shear.
    """
    e=_cube(ell)/3
    n=e.shape[0]
    L=coefficient.length
    g=_derivatives(e,L)
    grad=g["grad"];H=g["H"];lap=g["lap"];w=g["w"]
    waves=2*np.pi*np.fft.fftfreq(n,d=L/n)
    ks=np.meshgrid(waves,waves,waves,indexing="ij")
    k2=sum(k*k for k in ks)
    def lapf(x):
        return np.fft.ifftn(-k2*np.fft.fftn(x)).real
    def double_div(arr):
        fx=np.zeros(e.shape,dtype=complex)
        for i in range(3):
            for j in range(3):
                fx-=ks[i]*ks[j]*np.fft.fftn(arr[i,j])
        return np.fft.ifftn(fx).real
    def divf(vec):
        ft=sum(1j*ks[i]*np.fft.fftn(vec[i]) for i in range(3))
        return np.fft.ifftn(ft).real
    eH=np.asarray([[e*H[i][j] for j in range(3)]
                                      for i in range(3)])
    gv=np.asarray([[grad[i]*grad[j] for j in range(3)]
                                      for i in range(3)])
    v=np.asarray([lap*grad[i]-
                  sum(H[i][j]*grad[j] for j in range(3))
                                     for i in range(3)])
    D=2*coefficient.w_eta/coefficient.a_s**4
    return D*(g["T"]+8*double_div(eH)+8*lapf(e*lap)-
            2*lapf(w)+2*double_div(gv)+4*divf(v))

def quadratic_spatial_metric_source(ell,coefficient):
    """O(ell) spatial metric source -4D Delta²(ell/3)."""
    e=_cube(ell)/3
    _,_,lap=geometry_derivatives(e,coefficient.length)
    _,_,lap2=geometry_derivatives(lap,coefficient.length)
    D=2*coefficient.w_eta/coefficient.a_s**4
    return -4*D*lap2

def conformal_demo():
    n=24;L=2*np.pi
    a=np.arange(n)*L/n
    x,y,z=np.meshgrid(a,a,a,indexing="ij")
    ell=.11*(np.cos(3*x)+.7*np.cos(4*y)+.4*np.sin(5*z))
    zeros=np.zeros_like(ell)
    c=TraceCoefficients(bf2=3.,c_l2=1.,m_l2=.2,
        w_eta=.015,a_j=1.,a_s=1.,length=L)
    g=_derivatives(ell/3,L)
    D=2*c.w_eta/c.a_s**4
    mc=cubic_spatial_metric_source(ell,c)
    R2=exact_conformal_riemann2(ell/3,c)
    h=1e-4
    derivative=(exact_curvature_action(ell,h*zeros+ h*np.ones_like(ell),c)-
      exact_curvature_action(ell,-h*np.ones_like(ell),c))/(2*h)
    return dict(
      scope="exact conformal 3D supporting curvature + cubic Jordan spatial metric vertex",
      mean_quadratic_curvature_energy=float(D*np.mean(g["T"])),
      mean_cubic_metric_source=float(mc.mean()),
      exact_riemann2_mean=float(R2.mean()),
      exact_uniform_spatial_metric_derivative=float(derivative),
      cubic_vs_quadratic_mean_identity_residual=float(
             mc.mean()-D*np.mean(g["T"])),
      missing=["full nonlinear trace/ADM kinetic and gradient metric stress",
         "traceless structural mode Q_ij",
         "rolling F/G2/G3 and physical homogeneous backreaction",
         "complete cosmological scalar/tensor perturbation action",
         "no-CDM CMB and lensing likelihood"])
if __name__=="__main__":
    print(json.dumps(conformal_demo(),indent=2))
