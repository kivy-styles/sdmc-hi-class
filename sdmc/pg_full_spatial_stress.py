"""SDMC Part III: 3D spatial-metric stress of an excited trace AND traceless well.

This extends Manuscript B Part III Section 106 supporting operator

  S_g = int N sqrt(h) {B f_L²/2[(D_n L)²-c_L²(DL)²-m_L² L²]
                       -W eta_K/2 Riem[gamma]²},

with a concrete symmetric 3-metric completion:
  h_ij = a_J²[exp(2epsilon)]_ij,
  gamma_ij = a_s² exp(2ell/3) [
                   exp(epsilon) exp(2Q) exp(epsilon)]_ij.

epsilon_ij(x) is an externally varied Jordan spatial strain; Q_ij is
a traceless symmetric structural relative strain and ell=tr L.
This metric congruence is positive definite and reduces exactly to the
manuscript's scalar conformal gamma_ij=a_s² exp(2q)delta_ij at Q=0.
The tensor completion beyond the flat/background reference frame is a
SPECIFIED covariantization for this calculation, not yet a uniquely
derived parent 4D action.

The 3D curvature is computed from the EXACT spatial metric via the
Christoffel connection and the identity Riem²=4Ricci²-R². Arbitrary
inhomogeneous Jordan spatial shear strains may be applied and their
action variation computed by centered finite differences.

The scalar kinetic/gradient part uses a frozen mixed-frame trace and
traceless tensor L, at zero ADM shift; its full time-dependent spatial
connection and shear/shift Lie transport require the parent action.
No conventional CDM species or density normalization is included.

All derivatives and background scales use Jordan spatial coordinates,
Jordan time t_c and structural a_s=S/S0 separately.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np

def _cube(x,n=None):
    x=np.asarray(x,dtype=float)
    if x.ndim!=3 or len(set(x.shape))!=1 or x.shape[0]<16:
        raise ValueError("periodic cubic 3D field with >=16 points")
    if n is not None and x.shape!=(n,n,n):
        raise ValueError("inconsistent periodic field shape")
    if not np.isfinite(x).all():
        raise ValueError("finite structural field required")
    return x

def _tensor(x,n):
    x=np.asarray(x,dtype=float)
    if x.shape!=(3,3,n,n,n):
        raise ValueError("tensor must have (3,3,n,n,n) shape")
    if not np.isfinite(x).all() or not np.allclose(x,x.swapaxes(0,1)):
        raise ValueError("real finite symmetric tensor required")
    return x

def _diff(x,axis,L):
    n=x.shape[-1]
    wave=2*np.pi*np.fft.fftfreq(n,d=L/n)
    fft=np.fft.fftn(x,axes=(-3,-2,-1))
    weights=(1j*wave).reshape(tuple(n if j==axis else 1
                                             for j in range(3)))
    return np.fft.ifftn(fft*weights,axes=(-3,-2,-1)).real

def _mat_exp_symmetric(x):
    """Batched symmetric matrix exp; input (3,3,n,n,n)."""
    arr=np.moveaxis(x,[0,1],[-2,-1])
    eig,v=np.linalg.eigh(arr)
    out=np.einsum("...ik,...k,...jk->...ij",
                    v,np.exp(eig),v,optimize=True)
    return np.moveaxis(out,[-2,-1],[0,1])

def _strain(strain,n):
    if strain is None:
        return np.zeros((3,3,n,n,n))
    x=np.asarray(strain,dtype=float)
    if x.shape==(3,3):
        x=np.broadcast_to(x[:,:,None,None,None],
                          (3,3,n,n,n)).copy()
    return _tensor(x,n)

def exact_riemann_squared_3d(gamma,L):
    """Exact Ricci and Riemann² of an arbitrary positive spatial metric.

    3D identity R_ijkl R^ijkl=4R_ij R^ij-R². Fourier derivatives
    on a periodic domain; nonlinear aliasing decreases with resolution.
    """
    gamma=np.asarray(gamma,dtype=float)
    if gamma.ndim!=5 or gamma.shape[0:2]!=(3,3):
        raise ValueError("metric requires 3x3xNxNxN array")
    n=gamma.shape[-1]
    gamma=_tensor(gamma,n)
    gm=np.moveaxis(gamma,[0,1],[-2,-1])
    if np.min(np.linalg.eigvalsh(gm))<=0:
        raise ValueError("gamma must be symmetric positive definite")
    inv=np.moveaxis(np.linalg.inv(gm),[-2,-1],[0,1])
    dgam=np.array([_diff(gamma,i,L) for i in range(3)])
    T=np.empty((3,3,3,n,n,n))
    for i in range(3):
        for j in range(3):
            for l in range(3):
                T[i,j,l]=dgam[i,j,l]+dgam[j,i,l]-dgam[l,i,j]
    Christoffel=.5*np.einsum("kl...,ijl...->kij...",
                             inv,T,optimize=True)
    deriv=np.empty((3,3,n,n,n))
    other=np.empty_like(deriv)
    for i in range(3):
        for j in range(3):
            deriv[i,j]=sum(_diff(Christoffel[k,i,j],k,L)
                                                  for k in range(3))
            other[i,j]=sum(_diff(Christoffel[k,i,k],j,L)
                                                  for k in range(3))
    trace_gamma=np.asarray([
        sum(Christoffel[l,k,l] for l in range(3))
                                   for k in range(3)])
    bilinear=np.einsum("kij...,k...->ij...",
                      Christoffel,trace_gamma,optimize=True)
    nonlin=np.einsum("kil...,ljk...->ij...",
                    Christoffel,Christoffel,optimize=True)
    Ricci=deriv-other+bilinear-nonlin
    R=np.einsum("ij...,ij...->...",inv,Ricci,optimize=True)
    ricci2=np.einsum("ik...,jl...,ij...,kl...->...",
              inv,inv,Ricci,Ricci,optimize=True)
    riem2=4*ricci2-R*R
    return dict(riemann_squared=riem2,Ricci=Ricci,R=R)

@dataclass(frozen=True)
class SupportCoefficients:
    bf2:float=3.
    c2:float=1.
    m2:float=.2
    weta:float=.01
    a_j:float=1.
    a_s:float=1.
    length:float=2*np.pi

    def __post_init__(self):
        vals=np.array(list(vars(self).values()),dtype=float)
        if not np.isfinite(vals).all() or min(vals)<=0:
            raise ValueError("all coefficients positive in this stability test")

def supporting_gamma(ell,Q,coeff,strain=None):
    ell=_cube(ell)
    n=ell.shape[0]
    q=_tensor(Q,n)
    e=_strain(strain,n)
    hhalf=_mat_exp_symmetric(e)
    eq=_mat_exp_symmetric(2*q)
    core=np.einsum("ia...,ab...,bj...->ij...",
                   hhalf,eq,hhalf,optimize=True)
    return coeff.a_s**2*core*np.exp(2*ell/3)[None,None]

def action_spatial_parts(ell,ell_dot,Q,Qdot,coeff,strain=None):
    """Spatially averaged ADM action per a_J³ at N=1, shift=0.

    The trace/traceless kinetic, restoring and gradient operators are
    the quadratic mixed reference-frame L² contractions.
    The structural curvature is computed exactly from gamma.
    epsilon enters the Jordan volume, h^{-1} gradients, and gamma.
    """
    ell=_cube(ell);n=ell.shape[0]
    vel=_cube(ell_dot,n)
    Q=_tensor(Q,n);Qdot=_tensor(Qdot,n)
    eps=_strain(strain,n)
    eI=np.eye(3)[:,:,None,None,None]
    field=Q+(ell/3)[None,None]*eI
    speed=Qdot+(vel/3)[None,None]*eI
    grad=np.array([_diff(field,i,coeff.length)
                                       for i in range(3)])
    invh=_mat_exp_symmetric(-2*eps)/coeff.a_j**2
    vol=np.exp(np.einsum("ii...->...",eps,optimize=True))
    kin=.5*coeff.bf2*np.einsum("ij...,ji...->...",
                            speed,speed,optimize=True)
    mass=.5*coeff.bf2*coeff.m2*np.einsum("ij...,ji...->...",
                            field,field,optimize=True)
    grad_contr=np.array([
        [np.einsum("ab...,ba...->...",grad[i],grad[j],
                   optimize=True) for j in range(3)] for i in range(3)])
    gradval=.5*coeff.bf2*coeff.c2*np.einsum(
                  "ij...,ij...->...",invh,grad_contr,optimize=True)
    gamma=supporting_gamma(ell,Q,coeff,eps)
    R2=exact_riemann_squared_3d(gamma,coeff.length)["riemann_squared"]
    curvature=.5*coeff.weta*R2
    parts=dict(
      kinetic=float(np.mean(vol*kin)),
      mass=float(np.mean(vol*mass)),
      gradient=float(np.mean(vol*gradval)),
      curvature=float(np.mean(vol*curvature)))
    parts["action"]=parts["kinetic"]-parts["mass"]-parts["gradient"]-parts["curvature"]
    parts["density"]=sum(parts[k] for k in
                             ("kinetic","mass","gradient","curvature"))
    return parts

def spatial_metric_projection(ell,ell_dot,Q,Qdot,coeff,basis,h=2e-5):
    """Off-shell action derivative with respect to a specified ε_ij(x).

    h_ij=a_J²exp(2ε)_ij; at ε=0 the output is
        (1/a_J³V) δS/δε[basis] = <T^i_j basis^j_i>
    within the chosen supporting-space metric completion.
    A full gauge-invariant total stress requires G3, mapper,
    lapse/shift transport and higher-derivative ADM completion.
    """
    if h<=0:raise ValueError("positive difference step")
    p=action_spatial_parts(ell,ell_dot,Q,Qdot,coeff,h*basis)
    m=action_spatial_parts(ell,ell_dot,Q,Qdot,coeff,-h*basis)
    return (p["action"]-m["action"])/(2*h)

def pressure_fixed_mapper(parts):
    """Exact uniform-isotropic stress projection in this metric completion.

    Under h_ij -> e^{2ε}h_ij at *fixed mapper*, all modes of gamma
    scale e^{2ε}; hence curvature energy has P_K=+E_K/3.
    """
    return (parts["kinetic"]-parts["mass"]-parts["gradient"]/3+
            parts["curvature"]/3)

def stress_projection_diagnostic(ell,ell_dot,Q,Qdot,coeff):
    """Compute 6 independent average Jordan spatial metric projections."""
    shape=ell.shape
    vals={}
    for i in range(3):
        basis=np.zeros((3,3))
        basis[i,i]=1
        vals[f"diag_{i}"]=spatial_metric_projection(
                    ell,ell_dot,Q,Qdot,coeff,basis)
    for i,j in [(0,1),(0,2),(1,2)]:
        basis=np.zeros((3,3))
        basis[i,j]=basis[j,i]=.5
        vals[f"offdiag_{i}{j}"]=spatial_metric_projection(
                    ell,ell_dot,Q,Qdot,coeff,basis)
    p=action_spatial_parts(ell,ell_dot,Q,Qdot,coeff)
    meanpressure=sum(vals[f"diag_{i}"] for i in range(3))/3
    return dict(energy=p["density"],parts=p,
        stress_projections=vals,pressure_from_metric=meanpressure,
        pressure_analytic_fixed_mapper=pressure_fixed_mapper(p),
        pressure_check_error=meanpressure-pressure_fixed_mapper(p),
        anisotropic_diagonal_spread=max(vals[f"diag_{i}"] for i in range(3))-
                                 min(vals[f"diag_{i}"] for i in range(3)),
        missing=["complete 4D parent action and mapper transformation",
                 "time-dependent covariant transport of Q",
                 "full Jordan lapse/shift and F/G2/G3 mixing",
                 "physical no-CDM FLRW background and matter perturbations"])

def sample():
    n=20;L=2*np.pi
    a=np.arange(n)*L/n
    x,y,z=np.meshgrid(a,a,a,indexing="ij")
    ell=.055*(np.cos(2*x)+.6*np.cos(3*y)+.2*np.sin(4*z))
    ell_dot=.018*np.sin(2*x)
    Q=np.zeros((3,3,n,n,n))
    Q[0,0]=.025*np.cos(3*x)
    Q[1,1]=-.025*np.cos(3*x)
    qdot=np.zeros_like(Q)
    c=SupportCoefficients(bf2=3.,c2=.7,m2=.25,weta=.025,
                          a_j=1.2,a_s=1.4,length=L)
    return dict(scope="nonlinear spatial Riemann² + trace/traceless reference-frame strain; mapper fixed",
                **stress_projection_diagnostic(ell,ell_dot,Q,qdot,c))

if __name__=="__main__":
    print(json.dumps(sample(),indent=2))
