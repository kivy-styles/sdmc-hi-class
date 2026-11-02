"""Conserved-material conformal initial-data ADM constraint test (GR comparison).

This goes beyond density Fourier proxies while remaining short of the
user's requested fully physical SDMC primordial joint probability.

Given the SAME ensembles of putative after-transition mass positions
and profile widths, create q(x)=-Phi_J from the linear weak-field
Jordan Poisson source. Treat this q as an initial spatial *template*,
not an exactly solved nonlinear matter source.

Take a plane-symmetric periodic conformal 3-metric
    h_ij=e^{2q(x)}delta_ij  on T^3, x,y,z in [0,2pi).
Choose *constant mean curvature* initial extrinsic geometry
    K_ij=-H h_ij  (H>0 constant, sigma_ij=0)
with the dust 4-velocity normal to the slice (J_i=0).

For c=1 and cosmological constant zero, the 4D GR initial-data
Hamiltonian and momentum constraints are
    R3 + K² - K_ij K^ij =16piG rho,
    D_j(K^{ij}-h^{ij}K)=8piG J^i.
The momentum constraint vanishes IDENTICALLY for K_ij=-Hh_ij
and constant H, J_i=0.

The Hamiltonian is exactly solved by reconstructing
    rho_phys(x)=[R3(x)+6 H²]/(16piG).
We additionally enforce the conserved total material inventory
    M_u = integral_T3 rho_phys(x) sqrt(h) d³x.
Thus
    H²=(16piG M_u/V_phys-<R3>_phys)/6,
where V_phys=(2pi)³<e^{3q}>_coord.
At any positive H², this is an exactly constrained GR dust
initial slice IF reconstructed rho_phys>=0 everywhere.
The resulting rho_phys is the true matter density of this
GR comparison slice, not necessarily the Gaussian seed
profile (which fixed only the trial q).

No classical FLRW curvature k is assigned from these slices.
No cosmological evolution or Planck-to-radiation kinetic
matter transfer is solved. These slices are not necessarily
allowed initial states of the FULL SDMC F/G2/G3+UV action,
and this is not a CDM-like Jordan gravitational source.

The exact physical volume mean of the intrinsic R3 is
  <R3>_V=2 <e^q(q_x)^2>/(<e^{3q}>) >=0 (for a_s=1).
The CMC expansion theta=3H uniform and shear=0 => Q_D=0.
The Buchert geometric effective source for this GR dust
initial slice is therefore
    rho_g,av=-<R3>/(16piG)<=0,
    P_g,av=<R3>/(48piG).
Even allowing general kinematic Q_D, the *special* positive
pressureless Buchert branch requires R_D=3Q_D<0, which is
impossible in a strictly locked periodic scalar conformal
spatial geometry with R_D>=0.

For any R_D>0 and rho_g,av>0, Q_D<-R_D. Then
    w_g,av=P_g,av/rho_g,av
          = 1+4R_D/(3[-R_D-Q_D]) >1,
rather than 0. This DOES NOT rule out a separate action
derived physical supporting-space geometric stress, other
slicings, nonconformal Q_ij, or genuinely curved initial
spatial reference.

The CMC comparison uses illustrative normalized G=.05,
M_u=1 and box length 2pi, no late300 parameters.
It imposes only *necessary within chosen GR CMC slice*
admissibility, not full UV initial-state physicality.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import json
import numpy as np
from sdmc.flatness_probability_ensembles import Ensemble,spatial_positions,wilson_interval
from sdmc.flatness_joint_cold_necessary import JointThresholds,joint_modes

@dataclass(frozen=True)
class CMCSlice:
    grid:int=256
    length:float=2*np.pi
    gravity:float=.05
    material_mass:float=1.
    density_floor_tolerance:float=1e-12

    def __post_init__(self):
        if self.grid<128 or self.grid%2 or self.length<=0 or \
           self.gravity<=0 or self.material_mass<=0 or \
           self.density_floor_tolerance<0:
            raise ValueError("regular positive GR initial-data geometry required")

def build_gr_cmc_slices(positions,widths,e:Ensemble,c:CMCSlice=CMCSlice()):
    x=np.asarray(positions,float)
    widths=np.asarray(widths,float)
    if x.ndim!=2 or x.shape[1]!=e.count or widths.shape!=(len(x),) or \
         np.any(widths<=0) or not np.isfinite(x).all() or \
         not np.isfinite(widths).all() or e.modes>=c.grid//2:
        raise ValueError("valid samples and resolved k modes")
    n=x.shape[0]
    k=np.arange(1,e.modes+1,dtype=float)
    # The initial gaussian sheet templates each carry total M_u/count.
    coeff=np.exp(-1j*k[None,:,None]*x[:,None,:]).mean(axis=2)
    coeff*=np.exp(-.5*(widths[:,None]*k[None,:])**2)
    average_density=c.material_mass/c.length**3
    # q=-Phi and Delta q=-4pi G (rho_template-rhobar)
    # in a_s=a_J=1 normalized test geometry.
    spectrum=np.zeros((n,c.grid//2+1),dtype=np.complex128)
    spectrum[:,1:e.modes+1]=4*np.pi*c.gravity*average_density*(
           coeff/k[None,:]**2)
    freq=np.arange(c.grid//2+1,dtype=float)
    q=np.fft.irfft(c.grid*spectrum,n=c.grid,axis=1)
    gradient=np.fft.irfft(c.grid*spectrum*(1j*freq[None,:]),
                         n=c.grid,axis=1)
    laplacian=np.fft.irfft(-c.grid*spectrum*(freq[None,:]**2),
                         n=c.grid,axis=1)
    R=-np.exp(-2*q)*(4*laplacian+2*gradient**2)
    vweight=np.exp(3*q)
    V=c.length**3*vweight.mean(axis=1)
    Ravg=(R*vweight).mean(axis=1)/vweight.mean(axis=1)
    gradient_avg=2*np.mean(np.exp(q)*gradient**2,axis=1)/vweight.mean(axis=1)
    Hsq=(16*np.pi*c.gravity*c.material_mass/V-Ravg)/6
    rho=(R+6*Hsq[:,None])/(16*np.pi*c.gravity)
    physical_mass=V*(rho*vweight).mean(axis=1)/vweight.mean(axis=1)
    good=(Hsq>0)&(rho.min(axis=1)>=-c.density_floor_tolerance)
    # exact CMC flow theta=3H uniform, shear=0, thus Q_D=0.
    Buchert_Q=np.zeros(n)
    rho_geom=-(Ravg+Buchert_Q)/(16*np.pi*c.gravity)
    pressure_geom=(Ravg-3*Buchert_Q)/(48*np.pi*c.gravity)
    return dict(
       geometry_R3_mean=Ravg,
       gradient_identity=gradient_avg,
       curvature_scalar_min=R.min(axis=1),
       curvature_scalar_max=R.max(axis=1),
       proper_volume=V,
       H_squared=Hsq,
       reconstructed_min_rho=rho.min(axis=1),
       conserved_physical_Mu=physical_mass,
       Hamiltonian_constraint_max_abs_residual=np.max(
            np.abs(R+6*Hsq[:,None]-16*np.pi*c.gravity*rho),axis=1),
       Momentum_constraint_identically_zero=True,
       GR_CMC_initial_dust_admissible=good,
       Buchert_Q=Buchert_Q,
       Buchert_geometric_rho=rho_geom,
       Buchert_geometric_pressure=pressure_geom,
       positive_pressureless_Buchert_dust=(rho_geom>0)&(
                        np.abs(pressure_geom)<1e-11),
       physical_SDMc_initial_state_verified=False,
       physical_CDM_like_Jordan_response_verified=False)

def one_prior(kind,e=Ensemble(),th=JointThresholds(),
              c=CMCSlice(),chunk=400):
    rng=np.random.default_rng(e.seed)
    count={k:0 for k in (
        "flat","wells","flat_and_wells","gr_cmc_admissible",
        "flat_and_gr_cmc","flat_wells_gr_cmc","positive_buchert_dust")}
    smallest_mass_error=0.
    max_constraint_error=0.
    max_mean_identity_error=0.
    for start in range(0,e.draws,chunk):
        n=min(chunk,e.draws-start)
        x=spatial_positions(kind,e.count,rng,e,n)
        w=rng.uniform(e.min_width,e.max_width,n)
        a=joint_modes(x,w,e,th)
        b=build_gr_cmc_slices(x,w,e,c)
        F=a["flat_proxy"]; W=a["well_proxy"]; C=b["GR_CMC_initial_dust_admissible"]
        counts_here={"flat":F,"wells":W,"flat_and_wells":F&W,
                     "gr_cmc_admissible":C,"flat_and_gr_cmc":F&C,
                     "flat_wells_gr_cmc":F&W&C,
                     "positive_buchert_dust":b["positive_pressureless_Buchert_dust"]}
        for key,v in counts_here.items(): count[key]+=int(np.sum(v))
        smallest_mass_error=max(smallest_mass_error,float(
            np.max(abs(b["conserved_physical_Mu"]-c.material_mass))))
        max_constraint_error=max(max_constraint_error,float(np.max(
                  b["Hamiltonian_constraint_max_abs_residual"])))
        max_mean_identity_error=max(max_mean_identity_error,float(np.max(
            abs(b["geometry_R3_mean"]-b["gradient_identity"]))))
    results={name:dict(count=num,proportion=num/e.draws,
            wilson_95=wilson_interval(num,e.draws))
           for name,num in count.items()}
    return dict(prior=kind,draws=e.draws,geometry=asdict(c),
         counts=results,
         joint_given_gr_cmc_admissible=(
           count["flat_wells_gr_cmc"]/count["gr_cmc_admissible"]
           if count["gr_cmc_admissible"] else None),
         max_mass_conservation_residual=smallest_mass_error,
         max_hamiltonian_constraint_residual=max_constraint_error,
         max_curvature_average_identity_residual=max_mean_identity_error,
         GR_CMC_initial_data_not_FULL_SDMC=True,
         observational_Jordan_flatness_not_computed=True,
         physical_CDM_like_gravity_not_computed=True)

def positive_R_buchert_bound(R,Q,G=1.):
    """Exact effective rho,w from Section 107's Buchert comparison."""
    if G<=0:raise ValueError("G>0 required")
    rho=-(R+Q)/(16*np.pi*G)
    P=(R-3*Q)/(48*np.pi*G)
    return dict(R=float(R),Q=float(Q),rho_geom=float(rho),
       pressure_geom=float(P),
       w=float(P/rho) if rho!=0 else None,
       positive_energy=bool(rho>0),
       near_pressureless_positive=bool(rho>0 and abs(P/rho)<.01))

def diagnostic():
    e=Ensemble()
    res={name:one_prior(name,e) for name in
        ("independent_uniform","clustered","evenly_spaced")}
    return dict(scope="GR constraint-satisfying CMC post-handoff toy; not physical Planck-to-radiation SDMC probability",
       physical_two_requirement_probability_known=False,
       same_60000_prior_samples=True,
       illustration_of_positive_R_and_negative_Q=positive_R_buchert_bound(
           .01,-.03),
       GR_CMC_prior_filtered_counts=res,
       restrictions=[
         "initial q inferred from weak Poisson Gaussian templates; exact GR material rho reconstructed from Hamiltonian and has same total Mu",
         "only irrotational dust GR CMC initial data; actual Planck-to-radiation includes pressure, SDMC scalar and separate structural gamma fields",
         "CMC Kij=-Hh has zero initial kinematic Buchert Q and no positive geometric dust",
         "no time evolution of initial configurations, structural energy amplitude, Jeans/sound-speed, growth or Weyl potential",
         "no physically justified gauge-fixed prior, no actual late300 nonflat Jordan OmegaK from CMB/BAO"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
