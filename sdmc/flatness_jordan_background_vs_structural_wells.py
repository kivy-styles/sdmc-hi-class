"""SDMC: homogeneous Jordan FLRW background vs inhomogeneous structural spatial metric.

The *background* Jordan FLRW spatial metric may be flat even if the
independent supporting-space spatial metric has localized intrinsic
curvature. The FULL Jordan metric still supports local gravitational
perturbations. One cannot make the full Jordan metric exactly uniform
with nontrivial local material density without satisfying its
Hamiltonian, momentum and spatial stress constraints.

Manuscript B Part III §102 defines:
 gamma_ij=M(x,t)^2 h_ik exp(2L)^k_j,
 L^i_j=ell(x,t) delta^i_j/3+Q^i_j;
 for scalar-only Q=0:
 h_ij=aJ² exp(-2 Phi_J(x)) delta_ij,
 gamma_ij=aS² exp[2(m+ell/3-Phi_J)]delta_ij,
 where M(x,t)=Mbar(t)*exp(m(x,t)), aS=Mbar*aJ.

A singular mode m(x)=mA*cos(kx), ell(x)=ellA*cos(kx),
 Phi_J(x)=phiA*cos(kx) gives
 structural q amplitude eps_s=mA+ellA/3-phiA;
 Jordan q amplitude eps_J=-phiA.

Three exact illustrative possibilities:
 (i) Jordan EXACT homogeneous phiA=0, mA=0, ellA!=0
     -> nonuniform gamma, but cannot be exactly L=0 locked.
 (ii) exact spatial locking L=0, spatially uniform mapper m=0,
      phiA!=0 -> both Jordan and gamma have local perturbations,
      though each FLRW BACKGROUND is flat.
 (iii) exact lock L=0, phiA=0, mA!=0
       -> nonuniform gamma due to local mapper, which must be a
          dynamically allowed physical scalar perturbation.

The coordinate average of gamma is aS² I0(2eps_s)delta_ij,
whose scalar curvature is precisely zero (constant-coefficient
metric). But the *physical volume average of structural R3* is

 <R3[gamma]> =
 eps_s²*k²/aS² *
  [I0(eps_s)-I2(eps_s)]/I0(3eps_s) >=0.

The noncommutation R3[<gamma>]=0 != <R3[gamma]>_V
demonstrates why assigning global Jordan flatness by averaging
the structural metric is a CHOICE OF COARSE GRAINING, not a
derived equivalence in the current SDMC action.

Furthermore coordinate-metric and volume-matched effective
rulers differ:
 aJ_coord/aJ=sqrt(I0(2eps_s));
 aJ_volume/aJ=I0(3eps_s)**(1/3).
The correct physically observed background metric requires
full Jordan and supporting-space constraints and an averaging
prescription, not arbitrary smoothing.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import math
import numpy as np
from scipy.special import iv
from .flatness_posthandoff_curvature import sinusoid_exact

@dataclass(frozen=True)
class ScalarFrameMode:
    mapper_amplitude:float=0.
    strain_trace_amplitude:float=.3
    jordan_potential_amplitude:float=0.
    k:int=1
    aJ:float=1.
    Mbar:float=1.

    def __post_init__(self):
        floats=[self.mapper_amplitude,self.strain_trace_amplitude,
                self.jordan_potential_amplitude,self.aJ,self.Mbar]
        if not all(math.isfinite(x) for x in floats) or \
           self.aJ<=0 or self.Mbar<=0 or self.k<1 or \
           int(self.k)!=self.k:
            raise ValueError("regular Jordan/structural scalar Fourier data")

    @property
    def structural_amplitude(self):
        return (self.mapper_amplitude+self.strain_trace_amplitude/3
                -self.jordan_potential_amplitude)

    @property
    def jordan_amplitude(self):
        return -self.jordan_potential_amplitude

    @property
    def aS(self):
        return self.Mbar*self.aJ

def frame_curvature(case:ScalarFrameMode):
    eps=case.structural_amplitude
    s=sinusoid_exact(eps,mode=case.k,a_s=case.aS)
    j=sinusoid_exact(case.jordan_amplitude,mode=case.k,a_s=case.aJ)
    return dict(
        scenario=vars(case),
        background_jordan_K=0,
        jordan_full_exactly_homogeneous=case.jordan_potential_amplitude==0,
        local_structural_3R_mean=s.physical_mean_exact,
        local_structural_linear_R_rms=s.linear_rms,
        local_jordan_3R_mean=j.physical_mean_exact,
        local_jordan_linear_R_rms=j.linear_rms,
        structural_local_geometry_nontrivial=abs(eps)>1e-14,
        lock_L_is_zero=abs(case.strain_trace_amplitude)<1e-14,
        mapper_spatially_uniform=abs(case.mapper_amplitude)<1e-14,
        coordinate_metric_mean_factor=float(iv(0,2*eps)),
        averaged_metric_scalar_curvature=0.,
        coordinate_coarse_jordan_scale_ratio=float(np.sqrt(iv(0,2*eps))),
        volume_coarse_jordan_scale_ratio=float(iv(0,3*eps)**(1/3)),
        physical_average_curvature_not_equal_to_metric_average=bool(
            s.physical_mean_exact>1e-18),
        physical_Jordan_constraint_solved=False)

def diagnostic():
    independent=frame_curvature(ScalarFrameMode())
    lock=frame_curvature(ScalarFrameMode(
           mapper_amplitude=0,strain_trace_amplitude=0,
           jordan_potential_amplitude=.1))
    mapper=frame_curvature(ScalarFrameMode(
           mapper_amplitude=.1,strain_trace_amplitude=0,
           jordan_potential_amplitude=0))
    trivial=frame_curvature(ScalarFrameMode(
           mapper_amplitude=0,strain_trace_amplitude=0,
           jordan_potential_amplitude=0))
    return dict(
        scope="exact scalar spatial geometry and averaging, not 4D cosmological background or observed curvature prediction",
        exactly_homogeneous_jordan_independent_structural_well=independent,
        locked_with_jordan_local_perturbation=lock,
        locked_with_local_mapper=mapper,
        exact_lock_uniform_mapper_exact_jordan=trivial,
        missing=["full Jordan and structural spatial Einstein/ADM constraints",
                 "primordial microscopic q and ell amplitudes from thermal handoff",
                 "physical 3D coarse-graining/gauge and extrinsic curvature",
                 "actual observational Jordan Omega_k from nonflat dynamics"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
