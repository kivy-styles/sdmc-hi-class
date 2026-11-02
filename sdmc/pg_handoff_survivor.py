"""SDMC reduced structural trace: radiation-tracking handoff excitation survival.

A temporary structural curvature selector W=4B(1-B) can excite the
existing trace ell through an externally specified weak Jordan
perturbation u(t), even if ell and dot ell start at zero. As W
switches off, the *forced static response* tends to zero, but the
homogeneous free solution need not be zero. In a stable fast
mass-restoring mode, its adiabatic energy approaches a_J^-3.

This is a well-defined reduced-system result, not a proof of a
primordial dark-matter source, causal formation of the matter
spectrum, or of an observational Jordan gravitational enhancement.

The toy background is *decelerating radiation tracking*:
 a_J=(1+t/t0)^1/2, H_J=1/[2(t+t0)], a_s=a_J^p.
For p=2, a_s=1+t/t0 and thus d²a_s/dt²=0 (structural coasting).

B=1/(1+exp[-(t-t_center)/width]), W=4B(1-B);
B>0 on every finite t; the singular exact B=0 boundary is NOT used.

K=B f_L²/3; omega0²=m_L²+c_L² k²/a_J²;
Lambda=8 W eta_K k⁴/(3 B f_L² a_s⁴).

ddot ell+[3H_J+dot B/B] dot ell
         +(omega0²+Lambda)ell = -3 Lambda u(t).

The illustrative prescribed u(t)=u0/a_J is compatible with the
simple subhorizon Newtonian potential of a frozen contrast in
conserved nonrelativistic material; it is not a fully self-consistent
radiation-era Jordan potential. All microscopic SDMC initial data
and the response of the metric to ell remain undetermined.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class HandoffScenario:
    t0:float=5.
    center:float=5.
    width:float=.8
    p:float=2.
    k:float=2.
    eta:float=.025
    fL:float=1.
    m2:float=1.
    c2:float=.1
    u0:float=.01
    source_on:bool=True
    selector_on:bool=True

    def __post_init__(self):
        v=np.asarray([self.t0,self.center,self.width,self.p,
          self.k,self.eta,self.fL,self.m2,self.c2,self.u0])
        if not np.isfinite(v).all():
            raise ValueError("finite toy parameters required")
        if min(self.t0,self.width,self.p,self.k,self.eta,self.fL)<=0 or \
              self.m2<0 or self.c2<0 or self.m2+self.c2<=0:
            raise ValueError("positive regular stable branch and time scales")

def background(t,s:HandoffScenario):
    if t<0: raise ValueError("post-transition t>=0 toy domain only")
    B=float(1/(1+np.exp(np.clip(-(t-s.center)/s.width,-700,700))))
    W=4*B*(1-B) if s.selector_on else 0.
    aJ=float(np.sqrt(1+t/s.t0))
    aS=float(aJ**s.p)
    H=float(1/(2*(t+s.t0)))
    bdot_over_b=float((1-B)/s.width)
    Kl=B*s.fL*s.fL/3
    C=4*W*s.eta*s.k**4/aS**4
    lam=2*C/(9*Kl)
    omega0=s.m2+s.c2*s.k*s.k/aJ**2
    u=s.u0/aJ if s.source_on else 0.
    return dict(B=B,W=W,aJ=aJ,aS=aS,H=H,bdot_over_b=bdot_over_b,
               Kl=Kl,C=C,Lambda=lam,omega0=omega0,u=u,
               jordan_deceleration=1.)

def rhs(t,state,s):
    ell,velocity=state
    b=background(t,s)
    return np.array([velocity,
           -(3*b["H"]+b["bdot_over_b"])*velocity
           -(b["omega0"]+b["Lambda"])*ell
           -3*b["Lambda"]*b["u"]])

def energies(t,ell,velocity,s):
    b=background(t,s)
    free=.5*b["Kl"]*(velocity**2+b["omega0"]*ell**2)
    forced=b["C"]*abs(b["u"]+ell/3)**2
    return dict(free=float(free),selector=float(forced),
                total=float(free+forced),
                comoving_free_action=float(free*b["aJ"]**3))

def evolve(s=HandoffScenario(),t_end=60.,sample=None):
    if t_end<=0:raise ValueError("positive time interval required")
    if sample is None:
        sample=(0.,1.,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,40.,50.,60.)
    x=np.array(sample,dtype=float)
    if x.ndim!=1 or min(x)<0 or max(x)>t_end or len(x)<2 or \
        np.any(np.diff(x)<=0):
        raise ValueError("sorted sample times inside evolution interval")
    sol=solve_ivp(lambda t,y:rhs(t,y,s),(0.,t_end),[0.,0.],
          method="DOP853",dense_output=True,rtol=5e-11,atol=2e-12,
          max_step=.045)
    if not sol.success:raise RuntimeError(sol.message)
    records=[]
    for t in x:
        ell,vel=sol.sol(t)
        b=background(t,s)
        E=energies(t,ell,vel,s)
        records.append(dict(t=float(t),ell=float(ell),
              velocity=float(vel),B=b["B"],W=b["W"],
              Lambda=b["Lambda"],aJ=b["aJ"],aS=b["aS"],
              **E))
    late=np.array([a["comoving_free_action"] for a in records
                             if a["t"]>=20])
    ratio=float(np.ptp(late)/np.mean(late))
    return dict(scope="decaying-W forced trace creates free massive-mode tail on a toy decelerating Jordan radiation tracker",
      records=records,
      late_comoving_action_spread_fraction=ratio,
      free_tail_nonzero=bool(records[-1]["free"]>0),
      observational_gravitational_extra_source_verified=False)

def diagnostic():
    main=evolve()
    noforce=evolve(HandoffScenario(source_on=False))
    noselector=evolve(HandoffScenario(selector_on=False))
    double=evolve(HandoffScenario(u0=.02))
    return dict(
       main=main,
       no_perturbation_max_abs_ell=max(abs(x["ell"]) for x in noforce["records"]),
       no_selector_max_abs_ell=max(abs(x["ell"]) for x in noselector["records"]),
       double_source_free_energy_ratio=(double["records"][-1]["free"]/
                                       main["records"][-1]["free"]),
       missing=["actual Planck-to-radiation kinetic coefficient history",
               "derived post-handoff perturbation spectrum u_k",
               "physical full mapper and Jordan gravitational potential backreaction",
               "initial stochastic trace statistics and amplitude",
               "3D gauge-invariant averaged density, pressure and lensing source"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
