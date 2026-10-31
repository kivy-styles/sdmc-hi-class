"""SDMC locked-flat structural trace + accepted rolling G3 scalar: 2-field action.

On the exactly locked-flat *zero structural-background energy* branch, the
supporting-space trace ell is an existing structural gravitational degree
of freedom. Its quadratic ADM kinetic and curvature terms couple to the
unitary-gauge Jordan scalar curvature zeta.

The scalar-only rolling Horndeski lapse/shift constraints yield Qs, Fs.
The structural supporting action adds K/2*(D_n ell)^2 and
-W eta_K/2 Riem[gamma]^2; around Lbar=0 and flat gamma background
the lapse/shift first-order constraints are unchanged.

For each real cosine Fourier mode (spatial average cos^2=1/2):
  L2 = (Qs/2) dot(zeta)^2 +(K/4)dot(ell)^2
     -(Fs/2)(k/a_J)^2 zeta^2
     -(K/4)omega0^2 ell^2
     -(C/2)(zeta+ell/3)^2,
  C=4W eta k^4/a_s^4, omega0^2=m_L^2+c_L^2 k^2/a_J^2.

Unitary gauge has delta sigma=0 so the *relative spatial curvature*
is q=zeta+ell/3, since gamma_ij=M^2 h_ik (exp2L)^k_j and
h_ij=a_J^2 exp(2zeta)delta_ij.
This identity uses Jordan time, not a structural Hubble observed clock.

This is a quadratic locked-flat branch with NO independent homogeneous
rho_g. It is not the nonlinear inhomogeneously averaged source and not
the physical cosmological no-CDM Boltzmann system.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np
from scipy.integrate import solve_ivp
from sdmc.pg_constraint_schur import TraceAction, effective_kernel,trace_stiffness
from sdmc.pg_rolling_g3_adm import (
  RollingAction,kinetic_from_constraints,isolated_scalar_gradient,
  no_slip_gradient_identity)

@dataclass(frozen=True)
class CoupledBranch:
    rolling:RollingAction
    trace:TraceAction
    Hdot:float
    alpha_B_dot:float

    def __post_init__(self):
        if not np.isfinite([self.Hdot,self.alpha_B_dot]).all():
            raise ValueError("finite Hdot and alpha_Bdot")
        if abs(self.rolling.a-self.trace.a_j)>1e-12:
            raise ValueError("Jordan scale a must match in both sectors")
        if self.trace.K<=0:
            raise ValueError("positive structural trace kinetic coefficient")

    def scalar(self):
        return isolated_scalar_gradient(
          self.rolling,self.Hdot,self.alpha_B_dot)

def kinetic_and_stiffness_matrices(branch,k):
    """Matrices for L2=1/2 qdot^T A qdot-1/2 q^T M q, q=(zeta,ell).

    This is the scalar-only Horndeski + existing structural trace
    principal quadratic action about the unexcited locked-flat branch.
    Physical matter/pg-background terms would change both blocks.
    """
    if k<0 or not np.isfinite(k):raise ValueError("finite k>=0")
    scalar=branch.scalar()
    Qs=scalar["Qs"]; Fs=scalar["Fs"]
    trace=branch.trace
    v=trace_stiffness(k,trace)
    omega0=v["omega2"]
    C=v["C"]
    A=np.diag([Qs,trace.K/2])
    M=np.array([
        [Fs*(k/trace.a_j)**2 +C,C/3],
        [C/3,.5*trace.K*omega0+C/9]
      ],dtype=float)
    return dict(A=A,M=M,Qs=Qs,Fs=Fs,
        cs2_scalar=scalar["cs2"],C=C,omega0_squared=omega0,
        relative_curvature="zeta + ell/3",
        scope="isolated rolling scalar + locked-flat structural trace; no rho_g background")

def quadratic_mode_lagrangian(branch,k,zeta,ell,zeta_dot,ell_dot):
    v=kinetic_and_stiffness_matrices(branch,k)
    q=np.array([zeta,ell])
    dot=np.array([zeta_dot,ell_dot])
    return float(.5*dot@v["A"]@dot-.5*q@v["M"]@q)

def dispersion(branch,k):
    """Local frozen-coefficient normal modes and principal stability gate.

    Stability of the *full* nonlinear action, cosmological pg
    background and matter-coupled system is NOT established.
    """
    d=kinetic_and_stiffness_matrices(branch,k)
    A=d["A"];M=d["M"]
    if np.any(np.diag(A)<=0):
        return dict(k=k,stable_kinetic=False,
                    eigenvalues_omega_squared=None,
                    stable_frozen_modes=False,
                    Qs=d["Qs"],Fs=d["Fs"])
    inv=np.diag(1/np.sqrt(np.diag(A)))
    B=inv@M@inv
    omegas2=np.linalg.eigvalsh(B)
    return dict(
      k=k,Qs=d["Qs"],Fs=d["Fs"],C=d["C"],
      omega0_squared=d["omega0_squared"],
      kinetic_eigenvalues=list(np.diag(A)),
      eigenvalues_omega_squared=omegas2,
      stable_kinetic=True,
      stable_frozen_modes=bool(np.all(omegas2>=-1e-11)),
      determinant_stiffness=float(np.linalg.det(M)),
      trace_stiffness=float(np.trace(M)),
      caution="local frozen coefficients and isolated scalar background only")

def static_trace_elimination(branch,k):
    """Schur eliminate ell from *joint* structural/Jordan stiffness M.

    At fixed zeta, ell*=-M12/M22*zeta and the leftover stiffness
    is Fs k²/a² + Keff, where Keff is exactly the previous
    static supporting-space Schur kernel.
    """
    d=kinetic_and_stiffness_matrices(branch,k)
    M=d["M"]
    if M[1,1]<=0:
        raise ValueError("structural trace static stiffness vanishes")
    response=-M[0,1]/M[1,1]
    eff=M[0,0]-M[0,1]**2/M[1,1]
    predicted=d["Fs"]*(k/branch.trace.a_j)**2+effective_kernel(k,branch.trace)
    return dict(ell_per_zeta=response,
      scalar_static_stiffness=eff,
      predicted_stiffness=predicted,
      schur_identity_residual=eff-predicted,
      dynamic_ell_not_eliminated=True)

def twofield_energy(branch,k,q,qdot):
    d=kinetic_and_stiffness_matrices(branch,k)
    q=np.asarray(q,dtype=float);dot=np.asarray(qdot,dtype=float)
    return float(.5*dot@d["A"]@dot+.5*q@d["M"]@q)

def integrate_locked_flat_modes(branch,k,times,q0,qdot0):
    """Exact frozen-matrix two-mode Hamiltonian control (no damping).

    Tests eigenfrequencies and conservation of an autonomous quadratic
    Hamiltonian. Cosmological background time evolution is omitted.
    """
    d=kinetic_and_stiffness_matrices(branch,k)
    A=d["A"];M=d["M"]
    if np.any(np.diag(A)<=0):raise ValueError("ghost kinetic block")
    t=np.asarray(times,dtype=float)
    if t.ndim!=1 or len(t)<2 or not np.all(np.diff(t)>0) or t[0]<0:
        raise ValueError("increasing nonnegative times")
    q0=np.asarray(q0,dtype=float);dot=np.asarray(qdot0,dtype=float)
    if q0.shape!=(2,) or dot.shape!=(2,):
        raise ValueError("need two-mode initial amplitudes")
    KminusM=np.linalg.solve(A,M)
    def rhs(_,y):
        return np.concatenate((y[2:],-KminusM@y[:2]))
    sol=solve_ivp(rhs,(t[0],t[-1]),np.concatenate((q0,dot)),
                  t_eval=t,method="DOP853",rtol=2e-11,atol=2e-12,
                  max_step=.025)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(sol.message)
    energies=np.asarray([
        twofield_energy(branch,k,sol.y[:2,i],sol.y[2:,i])
        for i in range(len(t))])
    return dict(t=t,zeta=sol.y[0],ell=sol.y[1],
                zeta_dot=sol.y[2],ell_dot=sol.y[3],
                energy=energies,
                fractional_energy_range=float(
                   np.ptp(energies)/max(abs(energies[0]),1e-15)))

def example():
    rolling=RollingAction(F=1.05,F_sigma=.12,H=.8,sigma_dot=.4,
              k1=1.7,k2=.3,g=-.75,g_sigma=.15)
    trace=TraceAction(K=.6,c2=.7,mass2=.2,weta=.012,
                      a_j=1.,a_s=1.)
    b=CoupledBranch(rolling=rolling,trace=trace,
            Hdot=-.05,alpha_B_dot=-.08)  # No-Slip: alpha_Mdot=.04
    rows=[]
    for k in [.01,.1,.5,1.,2.,5.]:
        d=dispersion(b,k)
        s=static_trace_elimination(b,k)
        rows.append(dict(k=k,omega_squared=[
                   float(x) for x in d["eigenvalues_omega_squared"]],
                   stable=d["stable_frozen_modes"],
                   Schur_residual=s["schur_identity_residual"]))
    sol=integrate_locked_flat_modes(b,k=1.,
                   times=np.linspace(0,30,300),
                   q0=[.2,.03],qdot0=[0.,0.])
    return dict(scope="quadratic locked-flat geometric trace plus complete rolling G3 ADM kinetic",
        scalar_Qs=b.scalar()["Qs"],scalar_Fs=b.scalar()["Fs"],
        scalar_isolated_cs2=b.scalar()["cs2"],
        geometric_trace_K=trace.K,
        dispersion=rows,
        energy_fractional_range=sol["fractional_energy_range"],
        missing=["nontrivial inhomogeneous structural background",
                 "physical rho_g and its momentum, pressure and shear",
                 "matter and radiation perturbation blocks",
                 "full non-quasistatic lapse-shift time-dependent constraints",
                 "physical structural gravitational response and likelihoods"])
if __name__=="__main__":
    print(json.dumps(example(),indent=2))
