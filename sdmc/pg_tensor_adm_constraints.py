"""SDMC nonzero trace+traceless supporting strain: ADM energy and shift sources.

This is the existing quadratic tensor-strain kinetic/restoring/gradient
operator of Manuscript B Part III Section 106, with the *exact spatial
curvature invariant* of the positive-definite spatial tensor completion
defined by sdmc.pg_full_spatial_stress.

L^i_j=ell delta^i_j/3+Q^i_j, Tr Q=0.
At a flat Jordan spatial reference metric, use ADM Lie transport
  D_n L = [dot L - Lie_{Nshift} L]/N,
  (Lie_N L)^i_j =
     N^k partial_k L^i_j
      -L^k_j partial_k N^i+L^i_k partial_j N^k.

Vary the reduced action with respect to ADM lapse and shift:
  rho = (bf2/2) Tr(dotL²)+(bf2 m²/2)Tr(L²)
      +(bf2 c²/(2a_J²))Tr(grad L grad L)+Weta Riem²/2,
  j_i = bf2 [Tr(dot L partial_i L)
             -partial_j [dot L,L]^j_i].

The commutator source disappears in pure scalar trace models, but can
be nonzero for rotated traceless relative strain. These are exact
lapse/shift variations within the RESTRICTED supporting action at
N=1, shift=0 and flat Jordan h. The full parent action would also
include h-dependent frame connections, structural clock/mapper and
rolling F/G2/G3 constraints. No CDM is introduced.
"""
from __future__ import annotations
import json
import numpy as np
from sdmc.pg_full_spatial_stress import (
    SupportCoefficients,_cube,_tensor,_diff,
    supporting_gamma,exact_riemann_squared_3d)

def fields(ell,ell_dot,Q,Qdot,coef):
    ell=_cube(ell);n=ell.shape[0]
    vel=_cube(ell_dot,n);Q=_tensor(Q,n);Qdot=_tensor(Qdot,n)
    eye=np.eye(3)[:,:,None,None,None]
    L=Q+(ell/3)[None,None]*eye
    dotL=Qdot+(vel/3)[None,None]*eye
    grad=np.array([_diff(L,i,coef.length) for i in range(3)])
    return L,dotL,grad

def local_3d_energy(ell,ell_dot,Q,Qdot,coef):
    """Exact quadratic trace/tensor ADM lapse source at background N=1."""
    L,dotL,grad=fields(ell,ell_dot,Q,Qdot,coef)
    kin=.5*coef.bf2*np.einsum("ij...,ji...->...",
                            dotL,dotL,optimize=True)
    mass=.5*coef.bf2*coef.m2*np.einsum("ij...,ji...->...",
                            L,L,optimize=True)
    G=sum(np.einsum("ab...,ba...->...",grad[i],grad[i],
                   optimize=True) for i in range(3))
    grad_energy=.5*coef.bf2*coef.c2/coef.a_j**2*G
    gamma=supporting_gamma(ell,Q,coef)
    R2=exact_riemann_squared_3d(gamma,coef.length)["riemann_squared"]
    curv=.5*coef.weta*R2
    return dict(
       kinetic=kin,mass=mass,gradient=grad_energy,curvature=curv,
       rho=kin+mass+grad_energy+curv)

def full_trace_tensor_shift_momentum(ell,ell_dot,Q,Qdot,coef):
    """ADM canonical momentum including Q commutator Lie-transport term."""
    L,dotL,grad=fields(ell,ell_dot,Q,Qdot,coef)
    n=len(ell)
    canonical=np.array([
       np.einsum("ab...,ba...->...",dotL,grad[i],optimize=True)
                                        for i in range(3)])
    dotL_L=np.einsum("ik...,kj...->ij...",dotL,L,optimize=True)
    L_dotL=np.einsum("ik...,kj...->ij...",L,dotL,optimize=True)
    comm=dotL_L-L_dotL
    derivative=np.array([
       sum(_diff(comm[j,i],j,coef.length) for j in range(3))
                                    for i in range(3)])
    return dict(momentum=coef.bf2*(canonical-derivative),
             canonical=coef.bf2*canonical,
             commutator_divergence=coef.bf2*derivative,
             tensor_commutator=comm)

def action_3d_lapse_shift(ell,ell_dot,Q,Qdot,coef,lapse=None,shift=None):
    """Full N/shift dependence within fixed-reference supporting ADM action.

    All V_restoring terms use fixed flat h and exact Riemann2[gamma];
    N and N^i enter through N sqrt(h) and D_n L as prescribed.
    """
    ell=_cube(ell);n=ell.shape[0]
    L,dotL,grad=fields(ell,ell_dot,Q,Qdot,coef)
    N=np.ones_like(ell) if lapse is None else _cube(lapse,n)
    if np.min(N)<=0:raise ValueError("positive ADM lapse required")
    if shift is None:
        sh=np.zeros((3,n,n,n))
    else:
        sh=np.asarray(shift,dtype=float)
        if sh.shape!=(3,n,n,n) or not np.isfinite(sh).all():
            raise ValueError("3 components of finite ADM shift required")
    adv=sum(sh[k][None,None]*grad[k] for k in range(3))
    dshift=np.array([[_diff(sh[i],j,coef.length)
                                      for j in range(3)] for i in range(3)])
    move_down=np.einsum("kj...,ik...->ij...",L,dshift,optimize=True)
    move_up=np.einsum("ik...,kj...->ij...",L,dshift,optimize=True)
    lie=adv-move_down+move_up
    normal=(dotL-lie)/N[None,None]
    energy=local_3d_energy(ell,ell_dot,Q,Qdot,coef)
    kinetic=.5*coef.bf2*np.einsum("ij...,ji...->...",
                   normal,normal,optimize=True)
    potential=energy["mass"]+energy["gradient"]+energy["curvature"]
    return float(np.mean(N*(kinetic-potential)))

def demo():
    n=20;lenbox=2*np.pi
    a=np.arange(n)*lenbox/n
    x,y,z=np.meshgrid(a,a,a,indexing="ij")
    ell=.05*(np.cos(2*x)+.5*np.cos(3*y))
    vel=.025*np.sin(2*x)
    Q=np.zeros((3,3,n,n,n))
    Qdot=Q.copy()
    Q[0,0]=.03*np.cos(2*x)
    Q[1,1]=-Q[0,0]
    Qdot[0,1]=.02*np.sin(3*y)
    Qdot[1,0]=Qdot[0,1]
    c=SupportCoefficients(bf2=3.,c2=.7,m2=.2,
                   weta=.015,a_j=1.2,a_s=1.4,length=lenbox)
    energy=local_3d_energy(ell,vel,Q,Qdot,c)
    mom=full_trace_tensor_shift_momentum(ell,vel,Q,Qdot,c)
    return dict(scope="restricted ADM supporting trace+traceless lapse and shift source",
       mean_energy=float(energy["rho"].mean()),
       mean_momentum=[float(s.mean()) for s in mom["momentum"]],
       momentum_rms=float(np.sqrt(np.mean(mom["momentum"]**2))),
       commutator_divergence_rms=float(np.sqrt(np.mean(
                      mom["commutator_divergence"]**2))),
       source_requires_full_mapper_constraint=True,
       missing=["rolling G3/ADM scalar plus supporting frame",
           "full physical gauge invariant Tmu_nu",
           "nonzero on-shell geometric cosmological background",
           "primordial geometric power and no-CDM likelihoods"])
if __name__=="__main__":
    print(json.dumps(demo(),indent=2))
