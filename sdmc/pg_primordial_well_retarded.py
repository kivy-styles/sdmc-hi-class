"""SDMC source-supported primordial structural well: retarded response and shutdown.

This directly implements Manuscript B Part III Section 106 Eqs (106.3-106.5):
    ell_ddot+(3 H_J+Kdot/K)ell_dot
       +(omega0²+Lambda_k)ell = -3 Lambda_k (m_k-Phi_J),
    K=B f_L²/3, Lambda_k=8 W eta_K k^4/(9 K a_s^4).

There is NO direct material-trace source in this accepted REDUCED action
once W=0. The prescribed metric mismatch u=m-Phi_J is an external input:
it has NOT been solved from the full F/G2/G3+matter constraints.

A material Fourier diagnostic conserves M_u during redistribution.
For a normalized wrapped Gaussian in a periodic box L with comoving
width sigma, its cosine density mode is
  rho_k=2 M_u/(a_J³ L³) exp[-(k sigma)²/2].
With M_remaining=(1-f)M_u, M_inferred=fM_u, a two-width model gives
the exact corresponding sum of those two normalized profiles.

The dimensionless imposed u(x) is taken proportional to the
NORMALIZED material Fourier SHAPE solely as an illustrative geometric
mismatch history. The proportionality is NOT a derived coupling.

A compact-support W(x) turns off exactly: an excited trace thereafter
obeys the unforced damped equation, unless another material coupling
is supplied by a COMPLETE physical parent action.

Static field+curvature energy for the COSINE mode:
  E(ell)=K omega0² ell²/4+C (u+ell/3)²/2,
  C=4 W eta_K k^4/a_s^4,
  ell_eq=-3 Lambda u/(omega0²+Lambda),
  E_min=C omega0²/(2(omega0²+Lambda)) u².

At constant source and coefficients, the minimized bare field energy
must not be counted as extra gravitational mass without interaction
stress. If W->0 and omega0²>0, ell_eq and E_min both ->0.

The homogeneous tail can remain dynamically excited for finite time
but a damped stable free mode cannot maintain a permanent static
primordial depression without persistent forcing or nonlinear binding.
"""
from dataclasses import dataclass
import json
import numpy as np
from scipy.integrate import solve_ivp

def _sigmoid(x):
    return .5*(1+np.tanh(x/2))

@dataclass(frozen=True)
class PrimordialWell:
    H0:float=.035
    h:float=1.5
    p_v:float=1.5
    B_floor:float=.10
    B_width:float=.18
    f_L2:float=3.
    m2:float=1.
    c2:float=.4
    k:float=1.2
    W_peak:float=.20
    eta_K:float=1.
    pulse_center:float=0.
    pulse_halfwidth:float=.27
    transfer_center:float=0.
    transfer_width:float=.15
    sigma_core:float=.20
    sigma_distributed:float=1.80
    mismatch_initial:float=.075
    M_u:float=1.
    box_length:float=2*np.pi

    def __post_init__(self):
        arr=np.asarray(list(vars(self).values()),dtype=float)
        if not np.isfinite(arr).all():raise ValueError("finite parameters required")
        positives=[self.H0,self.h,self.p_v,self.B_width,self.f_L2,self.m2,
          self.pulse_halfwidth,self.transfer_width,self.sigma_core,
          self.sigma_distributed,self.M_u,self.box_length]
        if min(positives)<=0 or not 0<self.B_floor<1:
            raise ValueError("positive well parameters and 0<B_floor<1")
        if min(self.k,self.c2,self.W_peak,self.eta_K,self.mismatch_initial)<0:
            raise ValueError("nonnegative source and curvature parameters")

def material_fraction(e,x):
    return _sigmoid((x-e.transfer_center)/e.transfer_width)

def material_ledger(e,x):
    f=material_fraction(e,x)
    remaining=(1-f)*e.M_u
    inferred=f*e.M_u
    # Each Gaussian has unit integrated mass on the periodic torus.
    c=np.exp(-.5*(e.k*e.sigma_core)**2)
    d=np.exp(-.5*(e.k*e.sigma_distributed)**2)
    a=np.exp(x)
    fourier_density=2*(remaining*c+inferred*d)/(a**3*e.box_length**3)
    return dict(f=f,M_remaining=remaining,M_inferred=inferred,
      M_total=remaining+inferred,
      material_density_cosine=fourier_density,
      normalized_mode_shape=(1-f)*c+f*d,
      core_shape=c,distributed_shape=d)

def _B(e,x):
    v=_sigmoid(x/e.B_width)
    B=e.B_floor+(1-e.B_floor)*v
    Bprime=(1-e.B_floor)*v*(1-v)/e.B_width
    return B,Bprime

def _W(e,x):
    z=(x-e.pulse_center)/e.pulse_halfwidth
    if abs(z)>=1:return 0.
    return e.W_peak*(1-z*z)**2

def coefficients(e,x):
    B,Bprime=_B(e,x)
    K=B*e.f_L2/3
    a=np.exp(x)
    a_s=np.exp(e.p_v*x)
    H=e.H0*np.exp(-e.h*x)
    W=_W(e,x)
    C=4*W*e.eta_K*e.k**4/a_s**4
    lam=2*C/(9*K)
    omega0_sq=e.m2+e.c2*e.k**2/a**2
    mismatch=e.mismatch_initial*material_ledger(e,x)["normalized_mode_shape"]
    total=omega0_sq+lam
    return dict(K=K,Krun=Bprime/B,H=H,W=W,a=a,a_s=a_s,
        C=C,Lambda=lam,omega0_sq=omega0_sq,
        total_omega_sq=total,u=mismatch,
        source_rhs=-3*lam*mismatch,
        static_ell=-3*lam*mismatch/total,
        static_q=omega0_sq*mismatch/total,
        C_geom_effective=C*omega0_sq/total)

def static_minimum(e,x):
    c=coefficients(e,x)
    ell=c["static_ell"];u=c["u"]
    Ebare=c["K"]*c["omega0_sq"]*ell*ell/4
    Ecurv=c["C"]*(u+ell/3)**2/2
    predicted=.5*c["C_geom_effective"]*u*u
    return dict(E_bare=Ebare,E_curv=Ecurv,E_min=Ebare+Ecurv,
        predicted=predicted,minimum_residual=Ebare+Ecurv-predicted)

def driven_acceleration(e,x,A,A_x):
    c=coefficients(e,x)
    return -(3-e.h+c["Krun"])*A_x -c["total_omega_sq"]*A/c["H"]**2 +c["source_rhs"]/c["H"]**2

def evolve(e,x_i=-.50,x_f=.85,n=451,excited=False):
    """Solve a driven trace from the period before W switches on.

    Both cases start at the equilibrium of the zero-selector theory,
    ell=0. The 'excited' case inserts a controlled initial ell velocity,
    displaying a distinct free homogeneous contribution.
    """
    if not x_f>x_i or n<25:raise ValueError("ordered simulation grid")
    x=np.linspace(x_i,x_f,n)
    c0=coefficients(e,x_i)
    y0=[float(c0["static_ell"]), float(.08/c0["H"] if excited else 0.)]
    def rhs(x,y):
        return [y[1],driven_acceleration(e,x,y[0],y[1])]
    sol=solve_ivp(rhs,(x_i,x_f),y0,t_eval=x,
           rtol=2e-10,atol=3e-11,method="DOP853",max_step=.0015)
    if not sol.success:raise RuntimeError(sol.message)
    A=sol.y[0];Ax=sol.y[1]
    eq=np.array([coefficients(e,t)["static_ell"] for t in x])
    C=np.array([coefficients(e,t)["C"] for t in x])
    source=np.array([coefficients(e,t)["source_rhs"] for t in x])
    mat=np.array([material_ledger(e,t)["M_total"] for t in x])
    K=np.array([coefficients(e,t)["K"] for t in x])
    H=np.array([coefficients(e,t)["H"] for t in x])
    omega0=np.array([coefficients(e,t)["omega0_sq"] for t in x])
    u=np.array([coefficients(e,t)["u"] for t in x])
    Ekin=K*H**2*Ax**2/4
    Emass=K*omega0*A*A/4
    Ecurv=C*(u+A/3)**2/2
    return dict(x=x,ell=A,ell_x=Ax,equilibrium=eq,
        mismatch=A-eq,curvature_coefficient=C,source_rhs=source,
        material_total=mat,energy_kin=Ekin,energy_strain=Emass,
        energy_curvature=Ecurv,
        B=np.array([_B(e,t)[0] for t in x]),
        H=H,K=K)

def followup_free_tail(e,sol):
    """Check the homogeneous ODE only after compact selector is EXACTLY off."""
    x=np.asarray(sol["x"])
    start=e.pulse_center+e.pulse_halfwidth
    mask=x>start+1e-3
    if not np.any(mask):raise ValueError("need post-selector interval")
    return dict(
      all_W_exactly_zero=bool(np.all(sol["curvature_coefficient"][mask]==0)),
      all_source_exactly_zero=bool(np.all(sol["source_rhs"][mask]==0)),
      static_well_exactly_zero=bool(np.all(sol["equilibrium"][mask]==0)),
      initial_tail_ell=float(sol["ell"][np.where(mask)[0][0]]),
      end_tail_ell=float(sol["ell"][-1]),
      initial_tail_energy=float((sol["energy_kin"]+
                   sol["energy_strain"])[np.where(mask)[0][0]]),
      final_tail_energy=float((sol["energy_kin"]+
                   sol["energy_strain"])[-1]),
      geometric_energy_still_positive=bool((sol["energy_kin"]+
                   sol["energy_strain"])[-1]>0))


def locked_jordan_well_curvature(e,x,ell_override=None):
    """Physical local 3-curvature of gamma on a SCALAR Fourier slice.

    In Newtonian gauge, q_k=m_k+ell_k/3-Phi_J(k)=u_k+ell_k/3,
    delta R_s(k)=4k² q_k/a_s², Manuscript B Eq(102.5).
    Thus at W=0 and ell=0, a nonzero Jordan potential can still
    create local structural curvature through gamma=M²h.

    The curvature-selector ENERGY is different:
      E_K=(C/2)*q_k², and vanishes with W.
    A nonzero R_s is not by itself a new dark gravitating mass.
    """
    c=coefficients(e,x)
    ell=c["static_ell"] if ell_override is None else float(ell_override)
    q=c["u"]+ell/3
    return dict(q=q,structural_ricci_mode=4*e.k**2*q/c["a_s"]**2,
          selector_energy=.5*c["C"]*q*q,
          relative_strain=ell,
          locked_strain_zero=ell==0,
          W=c["W"],
          has_local_structural_curvature=abs(4*e.k**2*q/c["a_s"]**2)>0)

def amplitude_curvature_source(e,x):
    """The effective curvature mismatch kernel exactly vanishes with W."""
    c=coefficients(e,x)
    return dict(W=c["W"],C=c["C"],Lambda=c["Lambda"],
       driven_amplitude=c["static_ell"],
       curvature_mismatch=c["static_q"],geometric_response=c["C_geom_effective"])

def response_memory_demo():
    e=PrimordialWell()
    r=evolve(e)
    out=followup_free_tail(e,r)
    xpeak=e.pulse_center
    z1=coefficients(e,xpeak)
    zoff=coefficients(e,e.pulse_center+1.1*e.pulse_halfwidth)
    check=static_minimum(e,xpeak)
    energy=r["energy_kin"]+r["energy_strain"]+r["energy_curvature"]
    ledgers=[material_ledger(e,x) for x in [-.5,0,.85]]
    return dict(
      scope="Manuscript B Eq106.4 driven curvature-mismatch source with prescribed metric u, not closed matter/Jordan constraints",
      origin="accepted W-dependent source -3 Lambda (m-Phi), no direct material-to-L source after W off",
      initial_x=-.5,end_x=.85,
      peak=dict(W=z1["W"],K=z1["K"],Lambda=z1["Lambda"],
         u=z1["u"],static_ell=z1["static_ell"],
         static_energy=check["E_min"],
         static_energy_residual=check["minimum_residual"]),
      final=dict(W=zoff["W"],Lambda=zoff["Lambda"],
          forcing=zoff["source_rhs"],static_ell=zoff["static_ell"],
          ell=float(r["ell"][-1]),energy=float(energy[-1]),
          B=float(r["B"][-1])),
      material_ledger=ledgers,
      locked_jordan_well_after_W_off=locked_jordan_well_curvature(e,.75),
      post_handoff=out,
      max_ell=float(np.max(abs(r["ell"]))),
      late_ell_residual=float(abs(r["ell"][-1])),
      missing=[
        "physical primordial coupling u=m-Phi from total metric constraints",
        "Planck-boundary initial canonical source and strain correlators",
        "late material-sourced L term independent of W, if any",
        "fully coupled structural and material stress and 4D mapper",
        "global rho_g, cold clustering and gravitational lensing predictions"
      ])

if __name__=="__main__":
    print(json.dumps(response_memory_demo(),indent=2))
