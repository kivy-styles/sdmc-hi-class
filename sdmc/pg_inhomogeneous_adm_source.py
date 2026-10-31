"""SDMC structural wells: first off-shell ADM lapse and shift source in 3D.

The actual reduced scalar trace ADM action is used:
S_red=int dt d^3x N sqrt(h) [
  K/2 ((dot ell-N^i partial_i ell)/N)^2
  -K m_L^2 ell^2/2
  -K c_L^2 h^(ij)partial_i ell partial_j ell/2
  -(W eta_K/2)*Riem[gamma]^2 ].

On a flat Jordan spatial slice and the quadratic structural-curvature
approximation, gamma_ij=a_s^2 exp[2(ell/3+zeta)]delta_ij, and
V_K = D [sum_ij(q_,ij)^2 + (Delta q)^2], D=2W eta_K/a_s^4.

The code varies *this restricted quadratic ADM functional* in the lapse,
shift and relative spatial-curvature variable. It DOES NOT represent
the full off-shell metric stress of the nonlinear structural action:
the terms cubic in trace/metric perturbations, tensor Q_ij and the
rolling G3 ADM constraints must still be included around an inhomogeneous
background. In particular a nonzero averaged source requires backreaction
on the zeroth-order Friedmann constraint and cannot be appended to a
background already solved with rho_g=0.

All time and physical measurements are Jordan t_c and a_J; a_s=S/S0.
No empirical CDM density or cluster-mass abundance is inserted.
"""
from __future__ import annotations
import json
import numpy as np
from sdmc.pg_quadratic_source import TraceCoefficients
from sdmc.pg_quadratic_source_3d import geometry_derivatives
from sdmc.pg_energy_exchange_3d import _cube

def fields(ell,ell_dot,coefficient,zeta=None):
    ell=_cube(ell)
    n=ell.shape[0]
    velocity=_cube(ell_dot,n)
    z=np.zeros_like(ell) if zeta is None else _cube(zeta,n)
    grad,_,_=geometry_derivatives(ell,coefficient.length)
    q=ell/3+z
    _,hq,lap=geometry_derivatives(q,coefficient.length)
    gradient2=sum(g*g for g in grad)
    curvature_shape=sum(hq[i][j]**2 for i in range(3)
                        for j in range(3))+lap**2
    D=2*coefficient.w_eta/coefficient.a_s**4
    return dict(ell=ell,velocity=velocity,grad=grad,q=q,
                gradient2=gradient2,curvature_shape=curvature_shape,D=D)

def local_reduced_sources(ell,ell_dot,coefficient,zeta=None):
    """Sources from lapse/shift variation of restricted quadratic ADM.

    δS/δN=-sqrt(h)*rho, δS/δN^i=-sqrt(h)*j_i on flat N=1,shift=0.
    The full physical metric stress also has cubic and higher vertices.
    """
    v=fields(ell,ell_dot,coefficient,zeta)
    K=coefficient.K
    Ekin=.5*K*v["velocity"]**2
    Emass=.5*K*coefficient.m_l2*v["ell"]**2
    Egrad=.5*K*coefficient.c_l2/coefficient.a_j**2*v["gradient2"]
    Ecurv=v["D"]*v["curvature_shape"]
    rho=Ekin+Emass+Egrad+Ecurv
    momentum=np.stack([K*v["velocity"]*g for g in v["grad"]],axis=0)
    return dict(rho=rho,momentum=momentum,
       kinetic=Ekin,mass=Emass,gradient=Egrad,curvature=Ecurv)

def lapse_shift_action(ell,ell_dot,coefficient,N=None,shift=None,zeta=None):
    """Action density averaged over 3D comoving period, per a_J^3.

    N perturbation and contravariant ADM shift N^i are exact within
    the *quadratic scalar trace* truncation, at fixed spatial h=a_J²δ.
    The zeta perturbation enters only q in the curvature term;
    the h-volume and h-gradient cubic vertices are omitted here.
    """
    v=fields(ell,ell_dot,coefficient,zeta)
    shape=v["ell"].shape
    lapse=np.ones(shape) if N is None else _cube(N,len(v["ell"]))
    if np.min(lapse)<=0: raise ValueError("ADM lapse must be positive")
    if shift is None:
        shiftvec=np.zeros((3,)+shape)
    else:
        shiftvec=np.asarray(shift,dtype=float)
        if shiftvec.shape!=(3,)+shape or not np.isfinite(shiftvec).all():
            raise ValueError("shift must be 3 contravariant components")
    normal=(v["velocity"]-sum(shiftvec[i]*v["grad"][i]
                for i in range(3)))/lapse
    K=coefficient.K
    L=.5*K*normal**2-.5*K*coefficient.m_l2*v["ell"]**2
    L-=.5*K*coefficient.c_l2/coefficient.a_j**2*v["gradient2"]
    L-=v["D"]*v["curvature_shape"]
    return float(np.mean(lapse*L))

def source_curvature_variation(ell,coefficient,zeta=None):
    """Local δL_curv/δzeta=-4D Delta²q from periodic integration by parts.

    This is only the Riemann²-related linear equation source,
    not the full nonlinear metric equation's curvature stress.
    """
    z=np.zeros_like(ell) if zeta is None else zeta
    q=_cube(ell)/3+_cube(z,len(ell))
    _,_,lap=geometry_derivatives(q,coefficient.length)
    _,_,lap2=geometry_derivatives(lap,coefficient.length)
    D=2*coefficient.w_eta/coefficient.a_s**4
    return -4*D*lap2

def density_linearization(ell0,vel0,dell,dvel,coefficient):
    """Exact Gateaux linearization of the REDUCED quadratic energy density.

    δrho=K(vel0 dvel +m² ell0 dell+
           c²/a² grad(ell0).grad(dell))
        +2D[Hq0:Hdq+(Delta q0)(Delta dq)]
    with q0=ell0/3, dq=dell/3, at zeta=0.
    """
    e=_cube(ell0)
    n=e.shape[0]
    v=_cube(vel0,n)
    de=_cube(dell,n)
    dv=_cube(dvel,n)
    K=coefficient.K
    grad0,H0,lap0=geometry_derivatives(e,coefficient.length)
    grad1,H1,lap1=geometry_derivatives(de,coefficient.length)
    D=2*coefficient.w_eta/coefficient.a_s**4
    out=K*(v*dv+coefficient.m_l2*e*de+
       coefficient.c_l2/coefficient.a_j**2*
       sum(grad0[i]*grad1[i] for i in range(3)))
    out+=2*D/9*(sum(H0[i][j]*H1[i][j]
             for i in range(3) for j in range(3))+lap0*lap1)
    return out

def restricted_spatial_pressure(ell,ell_dot,coefficient):
    """Fixed-S homogeneous Jordan scale response from the reduced action.

    Not the physical global FLRW pg pressure; the full tensor spatial
    metric variation around nonzero inhomogeneous ell is not present.
    """
    s=local_reduced_sources(ell,ell_dot,coefficient)
    return s["kinetic"]-s["mass"]-s["gradient"]/3-s["curvature"]

def spectral_cosine_amplitude(field,length,k_long=1,width=0.):
    """Cosine amplitude of a periodic 3D field along x, optionally coarse grained."""
    f=_cube(field)
    n=f.shape[0]
    if width<0: raise ValueError("width must be >=0")
    xx=np.arange(n)*length/n
    basis=np.cos(2*np.pi*k_long*xx/length)[:,None,None]
    raw=2*float(np.mean(f*basis))
    damp=np.exp(-.5*(2*np.pi*k_long*width/length)**2)
    return raw*damp

def structural_source_demo():
    n=32;L=2*np.pi
    x=np.arange(n)*L/n
    xx,yy,zz=np.meshgrid(x,x,x,indexing="ij")
    k_short=6;k_long=1;eps=.08
    short=(np.cos(k_short*xx)+np.cos(k_short*yy)
                +np.cos(k_short*zz))/np.sqrt(3)
    ell=.18*(1+eps*np.cos(k_long*xx))*short
    vel=.03*np.sin(k_short*xx)
    coef=TraceCoefficients(bf2=3.,c_l2=1.,m_l2=.2,w_eta=.01,
                           a_j=1.,a_s=1.,length=L)
    s=local_reduced_sources(ell,vel,coef)
    rho=s["rho"]
    mean=float(rho.mean())
    amp=spectral_cosine_amplitude(rho,L,k_long)
    momentum=s["momentum"]
    z_source=source_curvature_variation(ell,coef)
    if mean<=0: raise RuntimeError("nonpositive test density")
    return dict(
        scope="restricted quadratic inhomogeneous ADM structural source, not full rho_g",
        grid=[n,n,n],
        mean_trace=float(ell.mean()),
        mean_reduced_density=mean,
        reduced_pressure_mean=float(restricted_spatial_pressure(
                         ell,vel,coef).mean()),
        long_density_amplitude=amp,
        long_density_contrast_amplitude=amp/mean,
        mean_momentum=[float(x.mean()) for x in momentum],
        momentum_rms=float(np.sqrt(np.mean(momentum**2))),
        structural_curvature_metric_source_mean=float(z_source.mean()),
        structural_curvature_metric_source_rms=float(
                    np.sqrt(np.mean(z_source**2))),
        mean_sources_nonzero_but_scalar_only_background_unmodified=False,
        missing=[
            "complete cubic Jordan metric stress of the excited state",
            "nonlinear Riemann squared and traceless structural metric",
            "background backreaction in F/G2/G3 ADM constraints",
            "gauge-invariant geometric density/momentum/shear",
            "initial well power spectrum and cosmological abundance",
            "no-CDM CMB and lensing predictions"
        ]
    )

if __name__=="__main__":
    print(json.dumps(structural_source_demo(),indent=2))
