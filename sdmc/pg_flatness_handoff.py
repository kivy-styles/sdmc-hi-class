"""SDMC Part III: coexistence test of early flatness selection and late geometric wells.

This module directly uses Manuscript B, Part I Sections 53-56 candidate UV action.
The structural/Jordan spatial lock and transient curvature selector are treated
as an EARLY BRANCH-CHOICE functional, not an unphysical continuous evolution
of the discrete mature FLRW curvature index K_J.

A positive curvature selector has cost
  V_K(κ,x) = 6 eta_K W(x) κ² / a_s(x)^4,
with K_J=0 its *energetic* homogeneous minimum when W>0.

After W→0, the existing relative strain L_ij still has a mass,
kinetic and gradient action proportional to B f_L². A nontrivial
occupied coherent mode need not be zero even if the ground-state
locking configuration is L=0.

For the scalar trace with a single real spatial cosine ell=A(t)cos(kx),
Ktrace(x)=B(x) fL²/3, Omega²=mL²+cL² k²/a_J²+
                          8 eta_K W(x) k⁴/(9 Ktrace a_s⁴).
The exact frozen-background Jordan mode equation in x=ln a_J is

  A_xx+(3-H_power + dlnKtrace/dx)A_x+Omega²/H_J² A=0.

This is an illustrative POST-Planck handoff evolution: B/W activation
are not predicted by SDMC's present 4D parent action.

A coherent rapidly oscillating mass-restoring mode can satisfy
 rho_ell≈Ktrace mL² A_envelope²/4 ≈ a_J^-3,
 even though A_envelope decreases when the locking coefficient B rises.

Neither the amount of initial excitation nor the physical observable
Jordan perturbation source/late300 likelihood is predicted.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class Handoff:
    B_floor:float=.05
    x_handoff:float=0.
    width:float=.20
    W_peak:float=.20
    fL2:float=3.
    etaK:float=1.
    mL:float=1.4
    cL2:float=.5
    k:float=.30
    H0:float=.012
    H_power:float=1.5
    p_v:float=1.5

    def __post_init__(self):
        v=np.asarray(list(vars(self).values()),dtype=float)
        if not np.all(np.isfinite(v)):raise ValueError("finite parameters")
        if not (0<self.B_floor<1) or self.width<=0:
            raise ValueError("finite positive early B floor")
        if self.W_peak<0 or self.etaK<0 or self.k<0 or self.cL2<0:
            raise ValueError("nonnegative structural selector and gradient")
        if min(self.fL2,self.mL,self.H0,self.H_power,self.p_v)<=0:
            raise ValueError("positive scalar kinetic, mass, expansion")

    def values(self,x):
        x=np.asarray(x,dtype=float)
        u=1/(1+np.exp(-np.clip(
            (x-self.x_handoff)/self.width,-300,300)))
        B=self.B_floor+(1-self.B_floor)*u
        Bprime=(1-self.B_floor)*u*(1-u)/self.width
        W=4*self.W_peak*u*(1-u)
        K=B*self.fL2/3
        Krun=Bprime/B
        ascale=np.exp(self.p_v*x)
        a=np.exp(x)
        H=self.H0*np.exp(-self.H_power*x)
        C=4*self.etaK*W*self.k**4/ascale**4
        Lambda=2*C/(9*K)
        m2=self.mL**2
        G=self.cL2*self.k**2/a**2
        Omega2=m2+G+Lambda
        return dict(u=u,B=B,Bprime=Bprime,K=K,Krun=Krun,W=W,
          a=a,a_s=ascale,H=H,C=C,Lambda=Lambda,
          gradient_frequency_squared=G,
          Omega=np.sqrt(Omega2),Omega2=Omega2,
          frequency_to_H=np.sqrt(Omega2)/H,
          curvature_selector_intensity=6*self.etaK*W/ascale**4)

    def hom_curvature_cost(self,x,kappa):
        return self.values(x)["curvature_selector_intensity"]*np.asarray(kappa)**2

def equation(e,x,y):
    v=e.values(x)
    friction=3-e.H_power+v["Krun"]
    return [y[1],-friction*y[1]-v["Omega2"]/v["H"]**2*y[0],
            y[3],-friction*y[3]-v["Omega2"]/v["H"]**2*y[2]]

def two_phase(e,x_i=-.6,x_f=.8,n=241):
    """Direct Jordan-time oscillator integration; two orthogonal phases.

    Initial amplitudes (1,0) and (0,omega/H) have equal energy and
    represent a phase average of the SAME mode, not extra species.
    """
    if x_f<=x_i or n<30:raise ValueError("ordered grid")
    initial=e.values(x_i)
    samples=np.linspace(x_i,x_f,n)
    y0=[1.,0.,0.,float(initial["Omega"]/initial["H"])]
    sol=solve_ivp(lambda x,y:equation(e,x,y),
          (x_i,x_f),y0,t_eval=samples,rtol=2e-9,atol=2e-10,
          method="DOP853",max_step=.002)
    if not sol.success:raise RuntimeError(sol.message)
    energy=[];potential=[];kinetic=[]
    envelope=[];invariant=[];expected=[];mass=[];curv=[]
    for j,x in enumerate(samples):
        v=e.values(x);a1,ax1,a2,ax2=sol.y[:,j]
        E_kin=v["K"]*v["H"]**2*(ax1**2+ax2**2)/8
        E_mass=v["K"]*e.mL**2*(a1**2+a2**2)/8
        E_gradient=v["K"]*v["gradient_frequency_squared"]*(a1**2+a2**2)/8
        E_curv=v["K"]*v["Lambda"]*(a1**2+a2**2)/8
        E=E_kin+E_mass+E_gradient+E_curv
        energy.append(E)
        mass.append(E_mass)
        curv.append(E_curv)
        kinetic.append(E_kin)
        potential.append(E_mass+E_gradient+E_curv)
        envelope.append(np.sqrt(a1*a1+a2*a2))
        invariant.append(E*v["a"]**3/v["Omega"])
        expected.append((initial["a"]/v["a"])**3*v["Omega"]/initial["Omega"])
    inv0=invariant[0];E0=energy[0]
    envelope0=envelope[0]
    predicted_envelope=(
        (initial["a"]**3*initial["K"]*initial["Omega"])/
        (e.values(x_f)["a"]**3*e.values(x_f)["K"]*
         e.values(x_f)["Omega"]))**.5
    return dict(
       x=samples,
       B=[float(e.values(x)["B"]) for x in samples],
       W=[float(e.values(x)["W"]) for x in samples],
       E_ratio=np.asarray(energy)/E0,
       E_mass_fraction=np.asarray(mass)/np.asarray(energy),
       E_curv_fraction=np.asarray(curv)/np.asarray(energy),
       pressure_ratio_reduced=(np.asarray(kinetic)-np.asarray(potential))/
                                np.asarray(energy),
       amplitude_envelope=np.asarray(envelope)/envelope0,
       WKB_E_ratio=np.asarray(expected),
       WKB_envelope_final=float(predicted_envelope),
       action_invariant=np.asarray(invariant)/inv0,
       H_over_frequency=np.array(
          [1/e.values(x)["frequency_to_H"] for x in samples]))

def formal_flatness(e,x,kappa_values=(-1,0,1)):
    """Only the curvature-sector energy gap, not a branch transition.

    A mature fixed FLRW curvature index κ is an invariant parameter,
    not a continuous dynamical oscillator in the late FLRW ansatz.
    Energetic flat minimum at W>0 is not a proof of dynamical
    branch selection or preferred topology.
    """
    v=e.values(x)
    gaps={str(k):float(e.hom_curvature_cost(x,k)) for k in kappa_values}
    if 0 not in kappa_values: raise ValueError("include exact flat branch")
    return dict(x=x,B=float(v["B"]),W=float(v["W"]),
      curvature_costs=gaps,
      flat_energy_minimum=all(v>=gaps["0"] for v in gaps.values()),
      homogeneous_κ_time_evolution_not_derived=True)

def amplitude_for_fraction(f,F,H,K,m):
    """Required *occupied mode* peak amplitude, not abundance prediction.

    ell=A_peak cos(mx)cos(kx), spatial and temporal average:
    <rho>=K m² A_peak²/4, rho_crit=3F H².
    Hence A_peak²=12 f F H²/(K m²).
    """
    if not all(np.isfinite([f,F,H,K,m])) or min(F,H,K,m)<=0 or f<0:
        raise ValueError("positive physical inputs and nonnegative fraction")
    return np.sqrt(12*f*F/K)*(H/m)

def post_handoff_consistency(e,x_end,frac=.25,ratio_mH=1e3,
                             F=1,K=1):
    """Formal occupancy consistency, not a computed late300 source."""
    v=e.values(x_end)
    tiny=amplitude_for_fraction(frac,F,1.,K,ratio_mH)
    return dict(W_fraction_of_peak=float(v["W"]/e.W_peak)
               if e.W_peak else 0.,
       B_final=float(v["B"]),
       flat_selector_off=float(v["W"])<.005*max(e.W_peak,.001),
       required_peak_ell_for_assumed_fraction=tiny,
       assumed_energy_fraction=frac,
       assumed_m_over_H=ratio_mH,
       empirical_abundance_fitted_or_predicted=False,
       note="illustrative amplitude inversion; full background and EFT cutoff remain untested")

def results():
    # These are illustrative post-Planck handoff profiles, NOT a prediction.
    m=Handoff(k=.05,W_peak=0.)
    with_selector=Handoff(k=.30,W_peak=.20)
    runs={}
    for name,e in [("mass_only",m),("transient_selector",with_selector)]:
        d=two_phase(e)
        j=-1
        runs[name]=dict(x_start=float(d["x"][0]),x_end=float(d["x"][-1]),
           B_start=d["B"][0],B_end=d["B"][j],
           W_start=d["W"][0],W_end=d["W"][j],
           energy_ratio=float(d["E_ratio"][j]),
           WKB_energy_ratio=float(d["WKB_E_ratio"][j]),
           invariant_end=float(d["action_invariant"][j]),
           mode_amplitude_ratio=float(d["amplitude_envelope"][j]),
           WKB_amplitude_ratio=d["WKB_envelope_final"],
           mass_energy_fraction_final=float(d["E_mass_fraction"][j]),
           curvature_energy_fraction_final=float(d["E_curv_fraction"][j]),
           final_pressure_over_E_restricted=float(d["pressure_ratio_reduced"][j]),
           initial_H_over_frequency=float(d["H_over_frequency"][0]))
    e=with_selector
    selection=[formal_flatness(e,x) for x in [-.8,0,.8,1.6]]
    return dict(scope="conditional prethermal-to-Jordan handoff of candidate SDMC spatial EFT",
      runs=runs,homogeneous_curvature_selector=selection,
      occupancy_inversion=post_handoff_consistency(e,x_end=1.6),
      core_finding="early positive curvature selector and late surviving mass-restoring well energies can coexist mathematically",
      missing=[
          "independent physical derivation of B(x), W(x) and the cutoff",
          "actual dynamical selection of K_J from nonflat initial space",
          "production spectrum of L_ij excitations and adiabatic survival",
          "full field equations with geometric background in Jordan frame",
          "observable no-CDM clustering, lensing and CMB likelihood"
      ])

if __name__=="__main__":
    print(json.dumps(results(),indent=2))
