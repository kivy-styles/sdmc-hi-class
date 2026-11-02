"""Conditional SDMC near-flatness fraction: explicit priors, not a physical probability.

GOAL: Test whether many possible post-handoff material distributions
yield small *averaged structural intrinsic curvature*, and keep
that quantity separate from the actual observational Jordan Omega_k.

INPUT A: Manuscript B conditional flat-background curvature mapping
  |Omega_k,J,0| = 0.0917953 |Ck0 Omega_k,R,P|.
For Ck0=1 and uniform X=Omega_k,R,P in [-B,B], fraction meeting
|Omega_k,J,0|<epsilon is min(1,epsilon/(0.0917953 B)).
This is a conditional hypothesis test, not a distribution from SDMC.

INPUT B: A 3D conformal spatial metric q(x) depending only on x,
periodic in x,y,z. The nonrelativistic Poisson source is a set of
n=8 equal material concentrations, positive Gaussian spatial
profiles uniform in y and z, conserving M_u exactly:
  rho_k ∝ exp(-k² sigma²/2) * [sum_i exp(-ik x_i)]/n.

This is a PLANE-SYMMETRIC positive-density toy, not a real 3D
population of isolated primordial wells. If q=-Phi_J is determined
by Newtonian Poisson on a flat reference slice, the *leading*
second-order physical mean structural Ricci curvature obeys
  <R3_s>_V ≈2/a_s² <|grad q|²>
          ∝ sum_{k>=1} |rho_k|²/k².
Its ratio to an initial one-lump fixed-width reference cancels
overall G,M_u,a_J,a_s factors:
  Rratio = [sum exp(-k²sigma²)|mean_i exp(-ik x_i)|²/k²] /
           [sum exp(-k²sigma_initial²)/k²].

The Monte Carlo assigns three UNJUSTIFIED example priors:
  independent_uniform: each x_i uniform on periodic circle;
  clustered: all x_i close to one random common center
             with Gaussian positional offsets std=.2;
  evenly_spaced: x_i near 2pi*i/8+random translation,
             with positional jitter std=.1.
Each run also samples fragment width sigma~Uniform(.1,.6).
All scenarios preserve exactly the same material budget.
No SDMC Planck-to-radiation dynamics were used to select
any of these measures or the degree of correlation.

The fraction Rratio<0.25 means at least 75% SUPPRESSION of
AVERAGED STRUCTURAL R3 in this restricted weak-field toy.
The fraction Rratio<0.1 means at least 90% suppression.
Neither is the probability of observed |Omega_k,J,0|<0.004.

This experiment rejects a universal claim that narrow wells always
imply macroscopic flatness, and explicitly demonstrates
measure-dependent fractions. Its outputs cannot establish
"1 in 10^60" or any other physical primordial probability.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import numpy as np

ARCHIVED_PRESENT_CURVATURE_TRANSFER=.0917953

@dataclass(frozen=True)
class Ensemble:
    draws:int=20000
    seed:int=4102
    count:int=8
    modes:int=100
    initial_width:float=.6
    min_width:float=.1
    max_width:float=.6
    clustered_offset:float=.2
    regular_jitter:float=.1

    def __post_init__(self):
        if self.draws<100 or self.count<2 or self.modes<20 or \
           self.initial_width<=0 or self.min_width<=0 or \
           self.max_width<self.min_width or self.clustered_offset<0 or \
           self.regular_jitter<0: raise ValueError("nonnegative physical toy geometry required")

def conditional_observable_probability(initial_half_range, epsilon,
           transfer=ARCHIVED_PRESENT_CURVATURE_TRANSFER,Ck=1.):
    """Uniform prior X~U[-B,B], not physical primordial SDMC measure."""
    if min(initial_half_range,epsilon,transfer)<=0 or Ck<0:
        raise ValueError("positive prior range, threshold, transfer; nonnegative mapper ratio")
    if Ck==0:return 1.
    return min(1.,epsilon/(initial_half_range*transfer*Ck))

def weak_curvature_ratio(positions,widths,initial_width=.6,modes=100):
    """Many positive equal-material Gaussian profiles; 3D sheet symmetry.

    positions shape (draws,n), widths shape (draws,).
    """
    x=np.asarray(positions,float)
    s=np.asarray(widths,float)
    if x.ndim!=2 or x.shape[1]<2 or s.shape!=(len(x),) or \
       min(initial_width,s.min())<=0 or modes<20 or \
       not np.isfinite(x).all() or not np.isfinite(s).all():
        raise ValueError("finite positions, >1 material unit and positive widths")
    k=np.arange(1,modes+1,dtype=float)
    amp=np.exp(-1j*k[None,:,None]*x[:,None,:]).mean(axis=2)
    weight=np.exp(-(k[None,:]*s[:,None])**2)/(k[None,:]**2)
    numerator=np.sum(abs(amp)**2 * weight,axis=1)
    denominator=float(np.sum(np.exp(-(k*initial_width)**2)/k**2))
    return numerator/denominator

def spatial_positions(kind,count,rng,e:Ensemble,draws):
    L=2*np.pi
    if kind=="independent_uniform":
        return rng.random((draws,count))*L
    if kind=="clustered":
        centers=rng.random((draws,1))*L
        return (centers+rng.normal(0,e.clustered_offset,(draws,count)))%L
    if kind=="evenly_spaced":
        translations=rng.random((draws,1))*L
        return (translations+np.arange(count)[None,:]*L/count+
                rng.normal(0,e.regular_jitter,(draws,count)))%L
    raise ValueError("unknown ensemble prior")

def wilson_interval(successes,n,z=1.95996398454):
    """Approximate two-sided 95% binomial interval."""
    phat=successes/n
    den=1+z*z/n
    mid=(phat+z*z/(2*n))/den
    half=(z*np.sqrt(phat*(1-phat)/n+z*z/(4*n*n)))/den
    return [float(max(0,mid-half)),float(min(1,mid+half))]

def evaluate(kind,e=Ensemble(),chunk=400):
    rng=np.random.default_rng(e.seed)
    arr=np.empty(e.draws,dtype=float)
    for start in range(0,e.draws,chunk):
        n=min(chunk,e.draws-start)
        x=spatial_positions(kind,e.count,rng,e,n)
        sigma=rng.uniform(e.min_width,e.max_width,n)
        arr[start:start+n]=weak_curvature_ratio(
            x,sigma,e.initial_width,e.modes)
    estimates={}
    for reduction,threshold in [(75,.25),(90,.1),(99,.01)]:
        successes=int(np.sum(arr<threshold))
        estimates[str(reduction)+"pct_reduction"]=dict(
            count=successes,probability=successes/e.draws,
            wilson_95=wilson_interval(successes,e.draws))
    return dict(n_draws=e.draws,assumed_distribution=kind,
        mean_ratio=float(arr.mean()),median_ratio=float(np.median(arr)),
        ratio_p05=float(np.quantile(arr,.05)),
        ratio_p95=float(np.quantile(arr,.95)),
        estimates=estimates)

def fixed_width_mean_ratio(width, count=8,initial_width=.6,modes=100):
    """Exact expected Rratio for independent uniform positions."""
    k=np.arange(1,modes+1,dtype=float)
    return float(np.sum(np.exp(-(k*width)**2)/k**2)/
              (count*np.sum(np.exp(-(k*initial_width)**2)/k**2)))

def diagnostic():
    e=Ensemble()
    priors=[1.,.1,.01]
    thresholds=[.004,.001]
    return dict(scope="conditional priors and 3D plane-symmetric weak-field structural averaging; not physical SDMC primordial probability",
         observed_planck_BAO_reference="OmegaK0=0.0007 +/- 0.0019, 68% CL",
         user_requested_probability_physically_determined=False,
         conditional_Jordan_curvature=[
          dict(initial_uniform_half_range=B,threshold=eps,
               probability=conditional_observable_probability(B,eps))
             for B in priors for eps in thresholds],
         geometric_mass_ensembles={
           name:evaluate(name,e) for name in
             ("independent_uniform","clustered","evenly_spaced")},
         fixed_width_analytic_expected_ratio=[
             dict(width=w,expected_ratio=fixed_width_mean_ratio(w))
             for w in (.1,.3,.6)],
         missing=[
           "measure over physically allowed prethermal SDMC states",
           "initial density correlations and profile-width dynamics",
           "full 4D ADM Hamiltonian and momentum constraints",
           "source-derived structural and Jordan spatial metric perturbations",
           "radiation pressure and thermal transfer generation",
           "map from averaged structural Ricci scalar to observed Jordan OmegaK",
           "present-day nonflat late300 likelihood propagation"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
