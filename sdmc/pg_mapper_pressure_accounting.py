"""SDMC Part III: mapper-fixed versus structural-ruler-fixed stress closure.

Two distinct *restricted off-shell* Jordan spatial variations have been
used in the SDMC supporting-space investigation:

(A) structural-ruler-fixed / fixed-geometry scale diagnostic:
    P_A = Ekin-Emass-Egrad/3-Ecurv.

(B) mapper-fixed spatial ADM metric variation, with
    gamma_ij=M² h_ik exp(2L)^k_j and M held fixed:
    P_B = Ekin-Emass-Egrad/3+Ecurv/3.

(B) is independently verified by varying the EXACT nonlinear 3D
supporting curvature metric, including its traceless strain.

The local reduced energy E evolves with the *physical* structural
history a_s=S/S0 and the Jordan time t_c. Consequently
   J_B=dot E+3H_J(E+P_B)
      =J_A+4 H_J Ecurv.

None of these partial-stress choices is the FULL covariant physical
P_g, because varying M, sigma, lapse, shift, h and supporting degrees
together requires a complete parent action and consistent background.

In WKB oscillatory curvature dominance:
  Ekin≈Ecurv≈E/2,
  w_B≈2/3,
  J_B/(H_J E)≈ [w_run-k_run-4(p_v-1)]/2.
For constant couplings, -2(p_v-1); E~a_J^(-3-2p_v).

Even an imposed frequency-compensation w_run=4p_v+k_run yields
E~a_J^-3 with w_B≈2/3 and J_B/(H_J E)≈+2. This is NOT
independently conserved pressureless matter. The mass-restoring
adiabatic regime alone has both w≈0 and J≈0 under these assumptions.
"""
from __future__ import annotations
import json
import numpy as np
from sdmc.pg_adiabatic_dilution import (
    WellEra, coefficients,mode_averaged_terms,
    on_shell_analytic_energy_x_derivative,
    wkb_log_dilution_slope,evolve_two_phase_ensemble)

def two_restricted_stresses(Ekin,Emass,Egrad,Ecurv):
    E=Ekin+Emass+Egrad+Ecurv
    A=Ekin-Emass-Egrad/3-Ecurv
    B=Ekin-Emass-Egrad/3+Ecurv/3
    return dict(E=E,P_ruler_fixed=A,P_mapper_fixed=B,
                pressure_difference=B-A,
                prediction=4*Ecurv/3)

def on_shell_mapper_exchange(era,x,A,Ax):
    """Exact local Ward identity in Jordan e-fold x for stress B."""
    terms=mode_averaged_terms(era,x,A,Ax)
    p=two_restricted_stresses(
       terms["Ekin"],terms["Emass"],terms["Egrad"],terms["Ecurv"])
    j_ruler=terms["Jred_over_H"]
    j_mapper=j_ruler+4*terms["Ecurv"]
    derivative=on_shell_analytic_energy_x_derivative(era,x,A,Ax)
    residual=derivative+3*(terms["E"]+p["P_mapper_fixed"])-j_mapper
    return dict(E=terms["E"],
          P_mapper=p["P_mapper_fixed"],
          P_ruler=p["P_ruler_fixed"],
          J_mapper_over_H=j_mapper,
          J_ruler_over_H=j_ruler,
          J_difference_over_H=j_mapper-j_ruler,
          J_difference_expected=4*terms["Ecurv"],
          action_ward_residual=residual)

def wkb_mapper_stress(era,x):
    """Phase-averaged fixed-mapper stress and exchange.

    Assumes Omega/H>>1 and adiabatic trace oscillations. No physical
    no-CDM abundance or late300 coupling law inserted.
    """
    c=coefficients(era,x)
    frac_g=c["G"]/c["Omega2"]
    frac_curv=c["Lambda"]/c["Omega2"]
    w=(frac_g+2*frac_curv)/3
    JoverH=((era.w_run-4*era.p_v-era.k_run+4)*
                   frac_curv/2)
    slope=-3*(1+w)+JoverH
    return dict(
       E_gradient_fraction=frac_g/2,
       E_curvature_fraction=frac_curv/2,
       P_over_E=w,
       J_over_H_E=JoverH,
       dilution_log_slope=slope,
       WKB_slope=wkb_log_dilution_slope(era,x),
       identity_residual=slope-wkb_log_dilution_slope(era,x),
       physical_global_dust=False)

def combine_modes_adiabatically(eras,weights,x):
    """Illustrative positive weighted superposition, not a fitted spectrum."""
    if len(eras)!=len(weights) or not eras:
        raise ValueError("matching nonempty spectrum/weights required")
    if np.any(np.array(weights)<0) or sum(weights)<=0:
        raise ValueError("positive nonzero energy weights")
    from sdmc.pg_adiabatic_dilution import adiabatic_formula
    E=np.array([w*adiabatic_formula(e,x)
                 for e,w in zip(eras,weights)])
    p=np.array([wkb_mapper_stress(e,x)["P_over_E"]
                 for e in eras])
    j=np.array([wkb_mapper_stress(e,x)["J_over_H_E"]
                 for e in eras])
    rho=float(E.sum())
    return dict(E_total=rho,
       w_effective=float(np.dot(E,p)/rho),
       J_over_H_E=float(np.dot(E,j)/rho),
       curvature_fraction=float(np.dot(E,
          [wkb_mapper_stress(e,x)["E_curvature_fraction"]
                       for e in eras])/rho),
       indeterminate_abundance=True,
       notes="all weights are free initial mode energies, not derived cosmological values")

def demo():
    samples={
       "mass_restoring":WellEra(m2=1.,G0=0.,Lambda0=0.),
       "gradient":WellEra(m2=0.,G0=1.,Lambda0=0.),
       "curvature_matter":WellEra(m2=0.,G0=0.,Lambda0=1.,
                                  p_v=1.5),
       "curvature_radiation":WellEra(m2=0.,G0=0.,Lambda0=1.,
                                    p_v=2.,h=2.),
       "curvature_frequency_compensated":WellEra(
             m2=0.,G0=0.,Lambda0=1.,p_v=1.5,w_run=6.)
    }
    rows={}
    for name,e in samples.items():
        rows[name]=wkb_mapper_stress(e,0)
    mixture=combine_modes_adiabatically(
        [samples["mass_restoring"],samples["curvature_matter"]],
        [.6,.4],x=.8)
    return dict(
      scope="fixed mapper vs fixed ruler stress: two restricted action variations, not full physical pg",
      wkb_cases=rows,
      arbitrary_two_mode_mixture=mixture,
      missing=["mapper field equation and transformation under ADM variations",
       "nonlinear trace and traceless ADM transport",
       "rolling G3 and structural background constraint closure",
       "observational normalization of mode distribution",
       "physical P_g and J_g in global Jordan FLRW"])

if __name__=="__main__":
    print(json.dumps(demo(),indent=2))
