"""Joint SDMC flatness and necessary local-well checks with one sample measure.

IMPORTANT SCOPE:
  This does NOT compute the probability of a physically allowed
  Planck-to-radiation initial configuration yielding an observed
  near-flat Jordan cosmology *and* a demonstrated CDM-like response.
  Neither a UV initial-state measure nor the full Jordan growth/
  lensing response is available from the previous flatness ensemble.

  It DOES compute the fraction passing BOTH *geometric prerequisites*
  within the SAME plane-symmetric conserved-material toy ensemble:
    A: weak average structural curvature R_after/R_before < .1
       (at least 90% reduction, not observed Omega_k,J).
    B: appreciable inhomogeneous locally concentrated density field,
       rms(delta_m) >=1, with at least 5% of density Fourier power
       in wavelengths k>=8 (short-scale spatial structure).

  B is NOT sufficient to establish gravitational dust:
    rho_g>0; w_g≈0; rho_g∝a_J^-3;
    small physical sound speed, linear growth and Weyl lensing
    response; correct CMB/BAO spectra and independent abundance.
  All these need a self-consistent physical stress/metric solution.

  The same explicit priors and pseudorandom stream as the previous
  60k Monte Carlo are used: count=8 equal positive mass portions,
  sigma~Uniform(.1,.6), 20,000 draws per prior, seed=4102,
  n_harmonics=100, initial one-lump sigma=.6.
  Evenly spaced prior contains translational jitter sigma=.1.
  Clustered prior has common center and within-cluster offset=.2.
  No primordial data or separate CDM species are introduced.

  Conformal scalar-only flat-reference geometry has
   <R3_s>_V=(2/a_s²)
      integral exp(q)|grad q|² / integral exp(3q) >=0.
  Under EXACT locked conformal map to the same Jordan dust
  average this implies R_D>=0. The special positive-energy
  *Buchert dust branch* in Manuscript B Section 107 instead
  requires R_D=3Q_D<0, since rho_g,av=-R_D/(12pi G)>0.
  Thus zero nontrivial configurations in THIS restricted
  locked conformal scalar geometry can simultaneously
  realize that particular positive geometric-dust branch.
  This is NOT a no-go for SDMC's separate direct action-level
  supporting-space density; Buchert GR dust is not interchangeable
  with the distinct structural gamma_ij stress tensor.

  The exact physical requested probability remains UNDEFINED
  without a normalized measure over genuine UV constraint
  solutions and a map to the observed present Jordan Omega_k
  and CDM-like growth/lensing criteria.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import json
import numpy as np
from sdmc.flatness_probability_ensembles import (
    Ensemble,spatial_positions,wilson_interval)

@dataclass(frozen=True)
class JointThresholds:
    mean_curvature_ratio:float=.1
    density_rms:float=1.
    short_scale_power_fraction:float=.05
    short_scale_min_k:int=8

    def __post_init__(self):
        if self.mean_curvature_ratio<=0 or self.density_rms<0 or \
           self.short_scale_power_fraction<0 or \
           self.short_scale_power_fraction>1 or self.short_scale_min_k<1:
            raise ValueError("valid positive curvature/locality thresholds required")

def joint_modes(positions,widths,e:Ensemble,threshold=JointThresholds()):
    """Shared-realization geometric proxy statistics, not CDM observables.

    delta_k=exp(-k²sigma²/2) mean_i exp(-ikx_i) for k>=1.
    For a real field with these Fourier coefficients:
      sigma_delta²=2sum_{k>=1}|delta_k|².
    The weak proper-volume mean of spatial scalar R3 is
    proportional to sum_k |delta_k|²/k², for the same mapper.
    """
    x=np.asarray(positions,dtype=float)
    widths=np.asarray(widths,dtype=float)
    if x.ndim!=2 or x.shape[1]!=e.count or widths.shape!=(len(x),) or \
       not np.isfinite(x).all() or not np.isfinite(widths).all() or \
       np.any(widths<=0):
        raise ValueError("finite positions and positive widths for one common sample set")
    if threshold.short_scale_min_k>e.modes:
        raise ValueError("short-mode cutoff exceeds Fourier resolution")
    k=np.arange(1,e.modes+1,dtype=float)
    phase=np.exp(-1j*k[None,:,None]*x[:,None,:]).mean(axis=2)
    power=np.exp(-(widths[:,None]*k[None,:])**2)*abs(phase)**2
    initial_reference=float(np.sum(np.exp(-(k*e.initial_width)**2)/k**2))
    curv_ratio=np.sum(power/k[None,:]**2,axis=1)/initial_reference
    fluct_rms=np.sqrt(2*np.sum(power,axis=1))
    short_power=np.sum(power[:,threshold.short_scale_min_k-1:],axis=1)
    all_power=np.sum(power,axis=1)
    local_power_share=np.divide(short_power,all_power,
                                out=np.zeros_like(short_power),
                                where=all_power>0)
    flat=curv_ratio<threshold.mean_curvature_ratio
    wells=(fluct_rms>=threshold.density_rms)&(
           local_power_share>=threshold.short_scale_power_fraction)
    return dict(curvature_ratio=curv_ratio,density_contrast_rms=fluct_rms,
        short_scale_power_share=local_power_share,
        flat_proxy=flat,well_proxy=wells,joint_proxy=flat&wells)

def one_prior(kind,e=Ensemble(),th=JointThresholds(),chunk=400):
    rng=np.random.default_rng(e.seed)
    counts=dict(flat=0,well=0,both=0)
    rsum=0.; rms_sum=0.; power_sum=0.
    for start in range(0,e.draws,chunk):
        n=min(chunk,e.draws-start)
        x=spatial_positions(kind,e.count,rng,e,n)
        widths=rng.uniform(e.min_width,e.max_width,n)
        d=joint_modes(x,widths,e,th)
        counts["flat"]+=int(np.sum(d["flat_proxy"]))
        counts["well"]+=int(np.sum(d["well_proxy"]))
        counts["both"]+=int(np.sum(d["joint_proxy"]))
        rsum+=float(np.sum(d["curvature_ratio"]))
        rms_sum+=float(np.sum(d["density_contrast_rms"]))
        power_sum+=float(np.sum(d["short_scale_power_share"]))
    probs={key:dict(count=v,probability=v/e.draws,
                      wilson_95=wilson_interval(v,e.draws))
           for key,v in counts.items()}
    return dict(prior=kind,total=e.draws,
          thresholds=asdict(th),fractions=probs,
          conditional_well_given_flat=(
             counts["both"]/counts["flat"] if counts["flat"]>0 else None),
          averages=dict(curvature_ratio=rsum/e.draws,
                        density_rms=rms_sum/e.draws,
                        short_power_share=power_sum/e.draws),
          physical_jordan_OmegaK_computed=False,
          physical_CDM_like_growth_and_lensing_computed=False)

def buchert_dust_sign_test(averaged_R3,gravity=1.):
    """Averaged GR *dust comparison branch*, not SDMC UV action rho_g.

    On rho_g,av>0 and P_g,av=0, require R_D=3Q_D<0.
    For q-only compact periodic conformal geometry R>=0:
    no positive-energy Buchert dust branch in this same metric.
    """
    if gravity<=0:raise ValueError("positive G required")
    q=averaged_R3/3
    rho=-averaged_R3/(12*np.pi*gravity)
    return dict(R_D=float(averaged_R3),Q_D_if_dust=float(q),
      rho_g_av_if_dust=float(rho),
      positive_geometric_dust_compatible=bool(rho>0),
      is_Jordan_Buchert_only=bool(True),
      distinct_action_route_not_excluded=bool(True))

def diagnostic():
    e=Ensemble()
    results={name:one_prior(name,e) for name in
             ("independent_uniform","clustered","evenly_spaced")}
    sensitivity={}
    for l in (.5,1.,2.):
        th=JointThresholds(density_rms=l)
        sensitivity[str(l)]={name:one_prior(name,e,th)["fractions"]["both"]["probability"]
                for name in results}
    return dict(scope="conditional geometric prerequisite overlap, NOT physical nearflatness-and-CDM probability",
      priors_and_thresholds=asdict(e),
      one_common_sample_measure_per_prior=True,
      joint_geometric_proxies=results,
      local_density_threshold_sensitivity=sensitivity,
      special_locked_conformal_Buchert_dust_sign=buchert_dust_sign_test(.01),
      actual_physically_allowed_primordial_probability_known=False,
      missing=[
         "normalized gauge-fixed measure over admissible Planck-to-radiation initial states",
         "initial 3D isotropic geometry and material correlations from UV dynamics",
         "causal radiation/matter redistribution and structural metric response",
         "independent Jordan Omega_k,0 for each nonflat trajectory",
         "positive physical action-derived geometric density P≈0 and a^-3 conservation",
         "correct gravity slip, potential growth, lensing and CMB/BAO likelihoods",
         "strongest scale-dependent local structural wave response through handoff"
      ])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
