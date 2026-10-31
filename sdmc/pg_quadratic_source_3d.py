"""3D scalar conformal curvature and averaged source for SDMC trace geometry.

Directly restricts the ADM gamma_ij=M^2 h_ik exp(2L)^k_j to the
scalar conformal perturbation q=m+ell/3-Phi_J in 3D flat background.
This is a second-order geometrical *energy proxy*, not a full
constraint-closed stress tensor, Weyl potential, or Boltzmann model.
"""
from dataclasses import dataclass
import json
import numpy as np
from sdmc.pg_quadratic_source import TraceCoefficients

def _cube(arr,n=None):
    v=np.asarray(arr,dtype=float)
    if v.ndim!=3 or min(v.shape)<16 or len(set(v.shape))!=1:
        raise ValueError("need cubic, periodic, 3-dimensional grid of size >=16")
    if n is not None and v.shape!=(n,n,n): raise ValueError("inconsistent 3D shape")
    if not np.isfinite(v).all(): raise ValueError("nonfinite source field")
    return v

def _waves(n,length):
    w=2*np.pi*np.fft.fftfreq(n,d=length/n)
    return tuple(np.asarray(k) for k in np.meshgrid(w,w,w,indexing='ij'))

def geometry_derivatives(field,length):
    """All spatial first/second derivatives on a 3D comoving periodic grid."""
    x=_cube(field)
    ks=_waves(x.shape[0],length)
    fx=np.fft.fftn(x)
    grad=[np.fft.ifftn(1j*ks[i]*fx).real for i in range(3)]
    hess=[[np.fft.ifftn(-ks[i]*ks[j]*fx).real for j in range(3)]
          for i in range(3)]
    lap=hess[0][0]+hess[1][1]+hess[2][2]
    return grad,hess,lap

def curvature_energy_3d(q,coeff):
    """V_K = 2 W eta/a_s^4 [sum_ij(q_,ij)^2 + (Lap q)^2].

    For linear conformally flat 3-geometry:
    R_ij=-q_,ij-delta_ij Lap(q), R=-4 a_s^-2 Lap(q),
    Riem^2=4 a_s^-4 [sum(q_,ij)^2 + (Lap q)^2].
    """
    _,H,lap=geometry_derivatives(q,coeff.length)
    hess_sq=sum(H[i][j]**2 for i in range(3) for j in range(3))
    return (2*coeff.w_eta/coeff.a_s**4)*(hess_sq+lap**2)

def quadratic_energy_3d(ell,coeff,ell_dot=None,m=None,phi=None):
    """Reduced physical Jordan-volume density for x,y,z trace profiles."""
    ell=_cube(ell)
    n=ell.shape[0]
    z=np.zeros_like(ell)
    vel=_cube(ell_dot,n) if ell_dot is not None else z
    mm=_cube(m,n) if m is not None else z
    pp=_cube(phi,n) if phi is not None else z
    q=ell/3+mm-pp
    grad,_,_=geometry_derivatives(ell,coeff.length)
    grad_sq=sum(g*g for g in grad)
    kinetic=.5*coeff.K*vel**2
    mass=.5*coeff.K*coeff.m_l2*ell**2
    spatial=.5*coeff.K*coeff.c_l2/coeff.a_j**2*grad_sq
    curv=curvature_energy_3d(q,coeff)
    return dict(kinetic=kinetic,mass=mass,gradient=spatial,curvature=curv)

def mean_quadratic_source_3d(ell,coeff,ell_dot=None,m=None,phi=None):
    v=quadratic_energy_3d(ell,coeff,ell_dot,m,phi)
    return sum(float(x.mean()) for x in v.values())

def gaussian_coarse_grain_3d(field,length,width):
    x=_cube(field)
    if width<0: raise ValueError("width>=0")
    ks=_waves(x.shape[0],length)
    k2=sum(kk**2 for kk in ks)
    return np.fft.ifftn(np.fft.fftn(x)*np.exp(-.5*width**2*k2)).real

def long_mode_contrast_3d(ell,coeff,k_long=1,width=.25):
    """Cosine contrast along x of the coarse grained quadratic energy."""
    terms=quadratic_energy_3d(ell,coeff)
    rho=sum(terms.values())
    mean=rho.mean()
    if mean<=0: raise ValueError("no geometric energy to normalize")
    contrast=gaussian_coarse_grain_3d(rho,coeff.length,width)/mean-1
    xx=np.arange(rho.shape[0])*coeff.length/rho.shape[0]
    basis=np.cos(2*np.pi*k_long*xx/coeff.length)[:,None,None]
    amp=2*float(np.mean(contrast*basis))
    return amp,contrast

def gradient_anisotropy_proxy_3d(ell,length):
    """Traceless <ell_,i ell_,j> correlation.

    This is not physical gravitational anisotropic stress until
    the full off-shell Jordan metric variation is carried out.
    """
    gr,_,_=geometry_derivatives(ell,length)
    tensor=np.array([[np.mean(gr[i]*gr[j]) for j in range(3)]
                     for i in range(3)],dtype=float)
    return tensor-np.trace(tensor)*np.eye(3)/3

def three_d_demo():
    n=48; L=2*np.pi
    axis=np.arange(n)*L/n
    xx,yy,zz=np.meshgrid(axis,axis,axis,indexing='ij')
    s=9; kL=1; eps=.1; sigma=.32
    short=(np.cos(s*xx)+np.cos(s*yy)+np.cos(s*zz))/np.sqrt(3)
    mod=(1+eps*np.cos(kL*xx))*short
    coeff=TraceCoefficients(bf2=3,c_l2=1,m_l2=.05,
                    w_eta=.01,a_j=1,a_s=1,length=L)
    amp,contrast=long_mode_contrast_3d(mod,coeff,k_long=kL,width=sigma)
    leading=2*eps*np.exp(-.5*(kL*sigma)**2)
    return {
       "scope":"periodic 3D scalar-trace energy proxy; not a physical no-CDM model",
       "grid":[n,n,n],"short_wave":s,"long_wave":kL,
       "short_field_mean":float(short.mean()),
       "mean_action_energy_arbitrary_units":mean_quadratic_source_3d(mod,coeff),
       "long_density_contrast_cosine_amplitude":amp,
       "separated_scale_approximation":float(leading),
       "unmodulated_gradient_anisotropy_norm":float(
           np.linalg.norm(gradient_anisotropy_proxy_3d(short,L))),
       "not_derived":["complete Jordan stress tensor","physical primordial 3D power spectrum",
                      "recombination transfer functions","geometric momentum/shear",
                      "observable potentials or lensing","rho_g cosmological normalization"]
       }

if __name__=="__main__": print(json.dumps(three_d_demo(),indent=2))
