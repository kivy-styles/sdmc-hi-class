"""SDMC Part III: adiabatic geometric well dilution from the ACTUAL reduced trace.

No CDM fluid is inserted. This implements the existing supporting-space
trace mode in an expanding Jordan background, *not* a global FLRW pg source.

For a real comoving cosine mode ell=A(t_c) cos(kx), with fixed background
metric driver u=0 and no inhomogeneous Jordan potential, the reduced
spatially averaged quadratic Lagrangian is

  L_bar = K/4 [ A_dot^2 - Omega^2 A^2],
  Omega²=m_L² + (c_L² k²/a_J²) + 2 C/(9K),
  C = 4 W eta_K k^4/a_s^4.

The structural and Jordan rulers are distinct:
 a_s=S/S0, d ln a_s/d ln a_J=p_v,
 while all time derivatives use observable Jordan t_c.

The restricted homogeneous-Jordan-scale stress is
  Ekin=K A_dot²/4,
  Emass=K m_L² A²/4,
  Egrad=K (c_L² k²/a_J²) A²/4,
  Ecurv=C A²/18 = K Lambda A²/4,
  E=Ekin+Emass+Egrad+Ecurv,
  P_red=Ekin-Emass-Egrad/3-Ecurv.

The exact ON-SHELL reduced energy exchange is
  Jred/H_J = (dlnK/dlna)(Emass+Egrad-Ekin)
           +(dln(Weta)/dlna-4p_v) Ecurv
for constant mass² and gradient coefficient c_L².
This is NOT the physical matter energy transfer of the complete
Jordan FLRW geometry.

In the high-frequency, slowly varying WKB regime, a uniform two-phase
ensemble has E~a_J^-3 Omega and P/E~G/(3 Omega²),
G=c_L²k²/a_J². Then for curvature domination and constant couplings,
Omega~a_s^-2, E~a_J^(-3-2p_v), NOT a_J^-3.

The conditional curvature-sector WKB dust scaling instead requires
d ln(Weta/K)/d ln a_J = 4p_v, so Omega stays constant. This is a
necessary tuning relation, NOT a result of the existing SDMC action.
"""
from dataclasses import dataclass,asdict
from math import exp
import json
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class WellEra:
    K0:float=1.
    m2:float=0.
    G0:float=0.
    Lambda0:float=1.
    H0:float=.003
    h:float=1.5
    p_v:float=1.5
    k_run:float=0.
    w_run:float=0.

    def __post_init__(self):
        vals=np.asarray(list(vars(self).values()),dtype=float)
        if not np.isfinite(vals).all():raise ValueError("finite era parameters")
        if self.K0<=0 or self.H0<=0:
            raise ValueError("K0,H0 must be positive")
        if min(self.m2,self.G0,self.Lambda0)<0:
            raise ValueError("positive stable structural frequency branch")
        if self.m2+self.G0+self.Lambda0<=0:
            raise ValueError("nonzero physical restoring frequency required")
        if self.h<=0 or self.p_v<=0:
            raise ValueError("expanding Jordan and structural eras expected")

def coefficients(e,x):
    """x=ln a_J, H_J=H0 exp(-h x), a_s/a_s0=exp(p_v x)."""
    K=e.K0*np.exp(e.k_run*x)
    H=e.H0*np.exp(-e.h*x)
    G=e.G0*np.exp(-2*x)
    Lam=e.Lambda0*np.exp((e.w_run-4*e.p_v-e.k_run)*x)
    Om2=e.m2+G+Lam
    return dict(K=K,H=H,G=G,Lambda=Lam,Omega2=Om2,
                a_j=np.exp(x),a_s_relative=np.exp(e.p_v*x),
                omega=np.sqrt(Om2),
                omega_over_H=np.sqrt(Om2)/H)

def well_x_acceleration(e,x,A,A_x):
    """Reduced EL equation in Jordan e-fold x, not structural time.

      A_xx+(3-h+k_run) A_x+Omega²/H² A=0
    """
    c=coefficients(e,x)
    return -(3-e.h+e.k_run)*A_x-c["Omega2"]/c["H"]**2*A

def mode_averaged_terms(e,x,A,A_x):
    c=coefficients(e,x)
    vel=c["H"]*A_x
    coeff=c["K"]/4
    Ekin=coeff*vel**2
    Emass=coeff*e.m2*A**2
    Egrad=coeff*c["G"]*A**2
    Ecurv=coeff*c["Lambda"]*A**2
    energy=Ekin+Emass+Egrad+Ecurv
    pressure=Ekin-Emass-Egrad/3-Ecurv
    J_over_H=(e.k_run*(Emass+Egrad-Ekin)
                +(e.w_run-4*e.p_v)*Ecurv)
    return dict(Ekin=Ekin,Emass=Emass,Egrad=Egrad,
        Ecurv=Ecurv,E=energy,P_restricted=pressure,
        Jred_over_H=J_over_H,H=c["H"],
        omega=c["omega"],omega_over_H=c["omega_over_H"])

def on_shell_analytic_energy_x_derivative(e,x,A,A_x):
    """Independent chain-rule derivative, using the EL acceleration."""
    c=coefficients(e,x)
    K,H,G,Lam=(c["K"],c["H"],c["G"],c["Lambda"])
    Omega2=c["Omega2"]
    ypp=well_x_acceleration(e,x,A,A_x)
    omega2prime=-2*G+(e.w_run-4*e.p_v-e.k_run)*Lam
    return (K/4)*(
       e.k_run*(H**2*A_x*A_x+Omega2*A*A)
       -2*e.h*H**2*A_x*A_x+2*H**2*A_x*ypp
       +omega2prime*A*A+2*Omega2*A*A_x)

def exact_reduced_ward_identity(e,x,A,A_x):
    terms=mode_averaged_terms(e,x,A,A_x)
    rate=on_shell_analytic_energy_x_derivative(e,x,A,A_x)
    residual=(rate+3*(terms["E"]+terms["P_restricted"])
              -terms["Jred_over_H"])
    return dict(residual=residual,energy_x_rate=rate,terms=terms)

def adiabatic_formula(e,x):
    """WKB averaged action invariant E a³/Omega = constant (mode ensemble).

    Works only if Omega/H >>1 and |dlnOmega/dx|*(H/Omega)<<1.
    """
    c=coefficients(e,x)
    c0=coefficients(e,0)
    return np.exp(-3*x)*c["omega"]/c0["omega"]

def wkb_equation_of_state(e,x):
    c=coefficients(e,x)
    return c["G"]/(3*c["Omega2"])

def wkb_exchange_over_H_E(e,x):
    c=coefficients(e,x)
    return ((e.w_run-4*e.p_v-e.k_run)*c["Lambda"]
            /(2*c["Omega2"]))

def wkb_log_dilution_slope(e,x):
    c=coefficients(e,x)
    # E ∝ a^-3 Omega, Omega²=m²+G+Lambda
    return (-3+(-2*c["G"]+
         (e.w_run-4*e.p_v-e.k_run)*c["Lambda"])
            /(2*c["Omega2"]))

def curvature_dust_running_condition(e):
    """Necessary in curvature-dominated WKB limit; NOT derived coupling flow."""
    return 4*e.p_v+e.k_run

def evolve_two_phase_ensemble(e,x_end=1.1,num=201):
    """Integrate two quadrature phases with equal initial mode energy.

    Taking a two-phase average removes rapid phase oscillations from
    rho and P to leading WKB order; these are two realizations of the
    SAME trace mode, not different physical species.
    """
    if x_end<=0 or num<20:raise ValueError("positive interval and >=20 samples")
    xpoints=np.linspace(0,x_end,num)
    c0=coefficients(e,0)
    om0=c0["omega"];h0=c0["H"]
    # Phase 1: A=1, A_dot=0; phase 2: A=0, A_dot=Omega0.
    y0=np.array([1.,0.,0.,om0/h0])
    def rhs(x,y):
        A1,V1,A2,V2=y
        return [V1,well_x_acceleration(e,x,A1,V1),
                V2,well_x_acceleration(e,x,A2,V2)]
    sol=solve_ivp(rhs,(0.,x_end),y0,t_eval=xpoints,
        method="DOP853",rtol=2e-10,atol=2e-11,
        max_step=.003)
    if not sol.success:raise RuntimeError(sol.message)
    out={key:[] for key in ["x","E_ratio","P_over_E",
         "J_over_H_E","WKB_E_ratio","WKB_P_over_E",
         "WKB_J_over_H_E","omega_over_H","ward_max_error"]}
    E0=None
    for j,x in enumerate(sol.t):
        t1=mode_averaged_terms(e,x,sol.y[0,j],sol.y[1,j])
        t2=mode_averaged_terms(e,x,sol.y[2,j],sol.y[3,j])
        E=.5*(t1["E"]+t2["E"])
        P=.5*(t1["P_restricted"]+t2["P_restricted"])
        J=.5*(t1["Jred_over_H"]+t2["Jred_over_H"])
        if E0 is None:E0=E
        res1=exact_reduced_ward_identity(e,x,sol.y[0,j],sol.y[1,j])
        res2=exact_reduced_ward_identity(e,x,sol.y[2,j],sol.y[3,j])
        out["x"].append(float(x))
        out["E_ratio"].append(float(E/E0))
        out["P_over_E"].append(float(P/E))
        out["J_over_H_E"].append(float(J/E))
        out["WKB_E_ratio"].append(float(adiabatic_formula(e,x)))
        out["WKB_P_over_E"].append(float(wkb_equation_of_state(e,x)))
        out["WKB_J_over_H_E"].append(float(wkb_exchange_over_H_E(e,x)))
        out["omega_over_H"].append(float(t1["omega_over_H"]))
        out["ward_max_error"].append(
           float(max(abs(res1["residual"]),abs(res2["residual"]))/
                 max(E,.000000000000000001)))
    return out

def demo():
    eras={
        "mass_restoring_matter":WellEra(m2=1.,G0=0.,Lambda0=0.),
        "spatial_gradient_matter":WellEra(m2=0.,G0=1.,Lambda0=0.),
        "structural_curvature_matter":WellEra(m2=0.,G0=0.,Lambda0=1.),
        "structural_curvature_radiation":WellEra(
                  m2=0.,G0=0.,Lambda0=1.,h=2.,p_v=2.,
                  H0=.001),
        "structural_curvature_compensated":WellEra(
                  m2=0.,G0=0.,Lambda0=1.,w_run=6.),
    }
    rows=[]
    for name,e in eras.items():
        r=evolve_two_phase_ensemble(e,x_end=1.1,num=201)
        j=-1
        rows.append(dict(case=name,
             p_v=e.p_v,coupling_w_run=e.w_run,
             end_a_ratio=float(np.exp(r["x"][j])),
             final_omega_over_H=r["omega_over_H"][j],
             E_ratio=r["E_ratio"][j],
             WKB_E_ratio=r["WKB_E_ratio"][j],
             P_over_E=r["P_over_E"][j],
             WKB_P_over_E=r["WKB_P_over_E"][j],
             J_over_H_E=r["J_over_H_E"][j],
             WKB_J_over_H_E=r["WKB_J_over_H_E"][j],
             exact_ward_max_error=max(r["ward_max_error"]),
             dustlike_coupling_running_condition=
                 curvature_dust_running_condition(e)))
    return {
       "scope":"local reduced supporting-space trace, not physical global rho_g",
       "results":rows,
       "principle":"rho_curvature,WKB ~ a_J^-3 a_s^-2 at fixed couplings",
       "not_derived":["global FLRW pressureless geometric fluid",
            "physical nonzero inhomogeneous structural background",
            "full Jordan ADM metric and scalar constraints around that state",
            "initial geometric mode spectrum and abundance",
            "CMB and weak-lensing no-CDM likelihood"]
    }

if __name__=="__main__":
    print(json.dumps(demo(),indent=2))
