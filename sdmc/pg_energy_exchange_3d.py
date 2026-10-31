"""3D Jordan-averaged energy and exchange of the reduced SDMC structural trace.

This is a spatially averaged action stress diagnostic, not full ADM+Horndeski
physical energy transfer. Full metric/lapse/shift constraints remain unknown.
"""
from __future__ import annotations
import json
import numpy as np
from sdmc.pg_energy_exchange import Epoch
from sdmc.pg_quadratic_source_3d import geometry_derivatives

def _cube(x,n=None):
    a=np.asarray(x,dtype=float)
    if a.ndim!=3 or len(set(a.shape))!=1 or min(a.shape)<16:
        raise ValueError("cubic 3D periodic fields required")
    if n is not None and a.shape!=(n,n,n): raise ValueError("shape mismatch")
    if not np.isfinite(a).all(): raise ValueError("nonfinite field")
    return a

def _geometry(ell,vel,u,udot,length):
    ell=_cube(ell); n=len(ell)
    vel=_cube(vel,n); u=_cube(u,n); udot=_cube(udot,n)
    grad_e,_,lap_e=geometry_derivatives(ell,length)
    q=ell/3+u
    Hq=geometry_derivatives(q,length)[1]
    lap_q=Hq[0][0]+Hq[1][1]+Hq[2][2]
    Hu=geometry_derivatives(udot,length)[1]
    lap_u=Hu[0][0]+Hu[1][1]+Hu[2][2]
    A=float(np.mean(ell**2)); V=float(np.mean(vel**2))
    G=float(np.mean(sum(g*g for g in grad_e)))
    Q=float(np.mean(sum(Hq[i][j]**2 for i in range(3)
                        for j in range(3))+lap_q**2))
    QU=float(np.mean(sum(Hq[i][j]*Hu[i][j] for i in range(3)
                        for j in range(3))+lap_q*lap_u))
    return dict(A=A,V=V,G=G,Q=Q,QU=QU,
                ell=ell,vel=vel,q=q,lap_e=lap_e,lap_q=lap_q)

def averaged_energy_pressure_3d(ell,vel,u,e,length):
    """3D reduced mean energy and restricted Jordan isotropic pressure."""
    g=_geometry(ell,vel,u,np.zeros_like(ell),length)
    D=2*e.weta/e.a_s**4
    kin=.5*e.K*g["V"]
    mass=.5*e.K*e.mass2*g["A"]
    grad=.5*e.K*e.c2/e.a_j**2*g["G"]
    curv=D*g["Q"]
    return dict(E=kin+mass+grad+curv,
                P_restricted=kin-mass-grad/3-curv,
                kinetic=kin,mass=mass,gradient=grad,curvature=curv,
                curvature_shape=g["Q"])

def trace_eom_3d(ell,vel,u,e,length):
    """Euler-Lagrange equation on periodic 3D grid; u external."""
    ell=_cube(ell); vel=_cube(vel,len(ell));u=_cube(u,len(ell))
    _,_,lap_e=geometry_derivatives(ell,length)
    q=ell/3+u
    _,_,lap_q=geometry_derivatives(q,length)
    _,_,lap2_q=geometry_derivatives(lap_q,length)
    acc=(-(3*e.H_j+e.Kdot/e.K)*vel
         -e.mass2*ell+e.c2/e.a_j**2*lap_e
         -(8*e.weta/(3*e.K*e.a_s**4))*lap2_q)
    return acc

def mean_jordan_exchange_3d(ell,vel,u,udot,e,length):
    """Spatially averaged Jred=E_dot+3H_J(E+P_red), on shell.

    The explicit geometry/clock and coupling source is

      .5 Kdot*(mass2<A>+c2/a²<G>-<v²>)
     +.5 K*(mass2dot<A>+c2dot/a²<G>)
     + Ddot <Hq²+(Lap q)²>
     +2D <Hq:Hu_dot+Lapq Lap(u_dot)>.

    D=2Weta/a_s⁴, Ddot=2/a_s⁴(Weta_dot-4H_s Weta).
    """
    g=_geometry(ell,vel,u,udot,length)
    D=2*e.weta/e.a_s**4
    Ddot=2*(e.wetadot-4*e.H_s*e.weta)/e.a_s**4
    return (
      .5*e.Kdot*(e.mass2*g["A"]+
                 e.c2/e.a_j**2*g["G"]-g["V"])
      +.5*e.K*(e.mass2dot*g["A"]+
               e.c2dot/e.a_j**2*g["G"])
      +Ddot*g["Q"]+2*D*g["QU"])

def demo():
    n=24;length=2*np.pi
    x=np.arange(n)*length/n
    xx,yy,zz=np.meshgrid(x,x,x,indexing="ij")
    ell=(np.cos(2*xx)+np.sin(3*yy)+np.cos(4*zz))*.2
    z=np.zeros_like(ell)
    u=.1*np.cos(xx)
    e=Epoch(a_j=1.,a_s=1.,H_j=.05,H_s=.075,
         K=1.,c2=.7,mass2=.2,weta=.03)
    energy=averaged_energy_pressure_3d(ell,z,u,e,length)
    J=mean_jordan_exchange_3d(ell,z,u,z,e,length)
    return {"scope":"3D periodic reduced-action average, not physical rho_g",
        "mean_energy":energy["E"],"mean_pressure_restricted":energy["P_restricted"],
        "exchange_residual":J,"relative_exchange_per_H_E":J/(e.H_j*energy["E"]),
        "p_v":e.p_v,
        "missing":["ADM constraints","scalar-Jordan interaction transfer",
                   "physical full momentum and shear","pressureless attractor",
                   "CMB and weak lensing"]}
if __name__=="__main__": print(json.dumps(demo(),indent=2))
