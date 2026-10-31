"""SDMC supporting-space trace: retarded response and generic Jordan coupling.

This continues the *reduced* ADM supporting-space trace action with
K = B f_L^2/3, omega0^2=m_L^2+c_L^2 k^2/a_J^2,
C=4 W eta_K k^4/a_s^4, a_s=S/S0.

The trace equation uses Jordan cosmic time t_c:
  ell_ddot + Gamma ell_dot + (omega0^2+Lambda)ell = -3 Lambda u
where u=delta ln(M)-Phi_J, Gamma=3H_J+Kdot/K,
Lambda=2C/(9K).

In a LOCAL FROZEN-COEFFICIENT approximation with harmonic u~e^(-i omega t):
  ell/u = -3 Lambda/(omega0^2+Lambda-omega^2-i Gamma omega)
  q/u = (omega0^2-omega^2-i Gamma omega)/D
  Kresp(omega,k) = C*q/u.

This is the causal forcing-response coefficient; an effective local
action at nonzero frequency requires treating time-nonlocal kernels and
the full dissipative/constraint system, not merely inserting a complex
spring constant in an instantaneous Lagrangian.

The exact rank-one algebra Hfull=H0+lambda vv^T is implemented for any
SUPPLIED nonsingular complex symmetric linear constraint operator H0.
The accepted rolling G3 Horndeski+ADM operator remains uncomputed.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np
from scipy.integrate import solve_ivp
from sdmc.pg_constraint_schur import TraceAction, trace_stiffness, effective_kernel

@dataclass(frozen=True)
class ResponseEpoch:
    """Freeze local coefficients at a Jordan epoch. Units must be consistent."""
    act: TraceAction
    H_j: float
    Kdot_over_K: float = 0.

    def __post_init__(self):
        if not np.isfinite([self.H_j,self.Kdot_over_K]).all():
            raise ValueError("finite H_j and Kdot/K required")

    @property
    def Gamma(self): return 3*self.H_j+self.Kdot_over_K

    def coefficients(self,k):
        v=trace_stiffness(k,self.act)
        return dict(
            omega0_squared=v["omega2"],Lambda=v["Lambda"],
            C=v["C"],omegaeff_squared=v["trace_frequency_squared"],
            Gamma=self.Gamma)

def trace_response(k,frequency,epoch):
    """Forced Fourier response, under frozen coefficient local approximation.

    harmonic convention: exp(-i omega t_c), signed Gamma >=0 is damping.
    At poles with Gamma=0 and real frequency, the response is undefined.
    """
    c=epoch.coefficients(k)
    if not np.isfinite([frequency]).all():
        raise ValueError("finite frequency")
    d=(c["omegaeff_squared"]-complex(frequency)**2
       -1j*c["Gamma"]*complex(frequency))
    if abs(d)<1e-13:
        raise ValueError("trace mode pole: static elimination invalid")
    ell_over_u=-3*c["Lambda"]/d
    q_over_u=1+ell_over_u/3
    return dict(
        ell_over_u=ell_over_u,
        q_over_u=q_over_u,
        K_response=c["C"]*q_over_u,
        denominator=d,
        response_frequency=frequency,
        static_kernel=effective_kernel(k,epoch.act),
        coefficients=c)

def static_adequacy(k,frequency,epoch):
    """Compare the exact local retarded response with static trace.

    Relative error normalized by max(C,1e-12); also reports dynamical
    ratios |omega|/omega_eff and |Gamma omega|/omega_eff^2.
    These are diagnostic, not complete QS validity conditions.
    """
    c=epoch.coefficients(k)
    r=trace_response(k,frequency,epoch)
    w2=c["omegaeff_squared"]
    if w2<=0: raise ValueError("zero restoring frequency")
    return dict(
        frequency_ratio=abs(frequency)/np.sqrt(w2),
        damping_ratio=abs(epoch.Gamma*frequency)/w2,
        absolute_kernel_difference=abs(r["K_response"]-r["static_kernel"]),
        scaled_kernel_difference=abs(r["K_response"]-r["static_kernel"])/
                max(c["C"],1e-12),
        dynamic_kernel_complex=r["K_response"],
        static_kernel=r["static_kernel"])


def local_trace_stability(k,epoch):
    """Frozen-coefficient homogeneous trace stability diagnostic.

    ell_ddot+Gamma ell_dot+omegaeff² ell=0 has roots
      lambda±=(-Gamma±sqrt(Gamma²-4omegaeff²))/2.
    This is a necessary local trace-sector gate. It cannot establish
    stability of the full G3+ADM metric/structural coupled system.
    """
    v=epoch.coefficients(k)
    Gamma=v["Gamma"];om2=v["omegaeff_squared"]
    roots=np.array([
       (-Gamma+np.lib.scimath.sqrt(Gamma*Gamma-4*om2))/2,
       (-Gamma-np.lib.scimath.sqrt(Gamma*Gamma-4*om2))/2],dtype=complex)
    max_growth=float(np.max(roots.real))
    p_K=epoch.Kdot_over_K/epoch.H_j if epoch.H_j!=0 else None
    return dict(
       Gamma=Gamma,
       omegaeff_squared=om2,
       roots=roots,
       max_local_growth_rate=max_growth,
       no_exponential_trace_growth=max_growth<=1e-12,
       Jordan_kinetic_running=p_K,
       expanding_epoch_bound_on_running=(
           "d ln K/d ln a_J >= -3" if epoch.H_j>0 else "not expanding"),
       necessary_only=True)

def step_trace_solution(times,k,epoch,u0):
    """Evolve causal trace with u(t>=0)=u0 and initial ell=ell_dot=0.

    This is a finite-time homogeneous-background test: coefficients are
    kept locally frozen, including the damping rate Gamma.
    """
    ts=np.asarray(times,dtype=float)
    if ts.ndim!=1 or ts.size<2 or np.any(np.diff(ts)<=0) or ts[0]<0:
        raise ValueError("ascending times, starting >=0")
    c=epoch.coefficients(k)
    if c["omegaeff_squared"]<=0: raise ValueError("restoring frequency required")
    if not np.isfinite(u0): raise ValueError("finite driver")
    def rhs(t,y):
        return [y[1],-epoch.Gamma*y[1]-
                c["omegaeff_squared"]*y[0]-3*c["Lambda"]*u0]
    sol=solve_ivp(rhs,(0,float(ts[-1])),[0.,0.],t_eval=ts,
                  method="DOP853",atol=2e-12,rtol=1e-11,max_step=.04)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(sol.message)
    ell=sol.y[0]
    q=u0+ell/3
    return dict(t=ts,ell=ell,ell_dot=sol.y[1],q=q,
        ell_static=-3*c["Lambda"]*u0/c["omegaeff_squared"],
        q_static=u0*c["omega0_squared"]/c["omegaeff_squared"])

def analytic_step_trace(times,k,epoch,u0):
    """Underdamped analytic retarded response. Gamma>=0 and Omega^2>0."""
    ts=np.asarray(times,dtype=float)
    c=epoch.coefficients(k)
    G=c["Gamma"];om2=c["omegaeff_squared"]-G*G/4
    if G<0 or om2<=0: raise ValueError("underdamped stable regime required")
    omega=np.sqrt(om2)
    lstat=-3*c["Lambda"]*u0/c["omegaeff_squared"]
    return lstat*(1-np.exp(-.5*G*ts)*
      (np.cos(omega*ts)+(G/(2*omega))*np.sin(omega*ts)))

def solve_rank_one(H0,v,b,lam):
    """Exact rank-one Dyson/Schur update on a supplied full linear operator.

    Hfull=H0+lam*v*v.T. Do NOT conj(v) in an analytic Fourier
    complex-symmetric linear-response action. lambda may be complex.

    Physical Phi/Psi response requires a known complete H0 from the
    rolling G3 + lapse/shift constrained action.
    """
    H=np.asarray(H0,dtype=complex)
    v=np.asarray(v,dtype=complex)
    b=np.asarray(b,dtype=complex)
    n=len(v)
    if H.shape!=(n,n) or b.shape!=(n,) or not np.allclose(H,H.T,
                           rtol=1e-11,atol=1e-11):
        raise ValueError("H0 must be square complex-symmetric; v,b vectors")
    if not (np.isfinite(H).all() and np.isfinite(v).all()
            and np.isfinite(b).all() and np.isfinite(lam)):
        raise ValueError("finite operator input required")
    x0=np.linalg.solve(H,b)
    y=np.linalg.solve(H,v)
    den=1+lam*np.dot(v,y)
    if abs(den)<1e-12:
        raise ValueError("modified constraint near a dynamical pole")
    x=x0-lam*np.dot(v,x0)/den*y
    return dict(response=x,baseline=x0,pole_denominator=den,
                trace_feedback=lam*np.dot(v,y),
                updated_operator=H+lam*np.outer(v,v))

def structural_jordan_rank_one(k,frequency,epoch,background_hessian,
                               source,relative_direction):
    """Couple action-derived retarded kernel to user-supplied H0.

    L_support=-Kresp(u)^2. Hence lambda=-2*Kresp for the
    signed Jordan Lagrangian Hessian H_J0, u=v.T x.
    This is a response diagnostic, not closed cosmological equations.
    """
    r=trace_response(k,frequency,epoch)
    coupled=solve_rank_one(background_hessian,relative_direction,
                           source,-2*r["K_response"])
    return dict(**coupled,structural_response=r)

def example():
    epoch=ResponseEpoch(TraceAction(K=1.,c2=.7,mass2=.4,
                                    weta=.035,a_j=1.,a_s=1.),H_j=.06)
    k=1.2
    coeff=epoch.coefficients(k)
    omega=np.sqrt(coeff["omegaeff_squared"])
    sample=[0.,.05*omega,.2*omega,.5*omega,omega,2*omega,10*omega]
    entries=[]
    for w in sample:
        r=trace_response(k,w,epoch)
        q=static_adequacy(k,w,epoch)
        entries.append(dict(frequency=w,
           frequency_over_omegaeff=w/omega,
           response_real=r["K_response"].real,
           response_imag=r["K_response"].imag,
           static_kernel=r["static_kernel"],
           scaled_static_difference=q["scaled_kernel_difference"]))
    t=np.linspace(0,40,400)
    step=step_trace_solution(t,k,epoch,u0=.35)
    return {
       "scope":"retarded reduced structural trace / supplied Jordan H0; NOT full rolling SDMC",
       "epoch":dict(K=epoch.act.K,H_J=epoch.H_j,Gamma=epoch.Gamma,
                    omega_eff=omega),
       "frequency_scan":entries,
       "step_driver":.35,
       "step_trace_at_end":float(step["ell"][-1]),
       "static_trace_target":float(step["ell_static"]),
       "missing":["full rolling G3 and G2 lapse-shift kernel",
                  "physical pg stress and Jordan metric source",
                  "no-CDM sigma/mode normalization",
                  "physical k,H units and slow-background validity",
                  "CMB and weak-lensing likelihood"]
    }
if __name__=="__main__":
    print(json.dumps(example(),indent=2))
