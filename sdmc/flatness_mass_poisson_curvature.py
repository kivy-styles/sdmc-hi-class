"""SDMC flatness: exactly conserved mass -> 3-D post-handoff structural curvature.

This is an explicit *specified-geometry* calculation answering the
primordial redistribution analogy, not a Planck-to-radiation initial
condition prediction or full radiation-era relativistic gravity solver.

Use a three-torus with periodic length L in each coordinate, the same
total material inventory in every arrangement, and Gaussian concentrations
as assumed spatial profiles. 'M_inferred' names the redistributed
material, not a new gravitational species.

At a chosen post-handoff classical slice, assume the *restricted locked
scalar response* m=ell=Q_ij=0, q_s=-Phi_J and
    gamma_ij=a_s² exp(2q_s) delta_ij.

Assume a weak subhorizon nonrelativistic Jordan Poisson source
    Laplacian(Phi_J) = 4pi G a_J² [rho_m-rhobar].
The source has zero average and hence has a unique zero-mean periodic
potential. This is NOT the exact Planck boundary or full radiation bath.

The exact intrinsic structural Ricci scalar is
 R3=-exp(-2q_s)/a_s²*(4 Laplacian q_s+2|grad q_s|²).
Its proper-volume average with periodic boundary conditions equals
 <R3>_V = (2/a_s²)
     [int exp(q_s)|grad q_s|² d³x]
      / [int exp(3q_s)d³x] >= 0.

The linear curvature is
 delta R3=+16pi G a_J²/a_s² [rho_m-rhobar],
which has zero coordinate mean; its fluctuations and the nonlinear
proper-volume mean are not the FLRW homogeneous K.

Crucial: conserved M_u by itself does not mandate suppression.
Redistribution into sufficiently *narrow* lumps can INCREASE average
curvature even though it reduces the number of large-scale coherent
material modes. Comparisons must retain the complete density profile.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import json
import numpy as np


@dataclass(frozen=True)
class MassGeometry:
    grid: int = 64
    box_length: float = 2 * np.pi
    G: float = 0.05
    aJ: float = 1.
    aS: float = 1.

    def __post_init__(self):
        if self.grid < 24 or self.grid % 2 or self.box_length <= 0 or \
           self.G <= 0 or self.aJ <= 0 or self.aS <= 0:
            raise ValueError("positive, resolved, even cubic periodic geometry")

    @property
    def cell_volume(self):
        return (self.box_length / self.grid)**3


def periodic_gaussian_density(centers, masses, width, p=MassGeometry()):
    """Return rho in a periodic box, with exact discrete masses.

    Periodic minimal-image Gaussian profiles; the discretized integral of
    each mass cloud is normalized separately, so mass is not duplicated.
    """
    pos=np.asarray(centers, dtype=float)
    weights=np.asarray(masses, dtype=float)
    if pos.ndim!=2 or pos.shape[1]!=3 or weights.shape!=(pos.shape[0],) or \
       len(weights)<1 or np.any(weights<0) or weights.sum()<=0 or \
       width<=0 or not np.isfinite(width) or \
       not np.isfinite(pos).all() or not np.isfinite(weights).all():
        raise ValueError("n centers Nx3, positive total mass, nonnegative masses and width")
    n=p.grid; L=p.box_length
    axis=np.arange(n,dtype=float)*L/n
    xyz=np.meshgrid(axis,axis,axis,indexing="ij")
    rho=np.zeros((n,n,n),dtype=float)
    for center,mass in zip(pos,weights):
        delta=[(coord-x+L/2)%L-L/2 for coord,x in zip(xyz,center)]
        cloud=np.exp(-.5*sum(d*d for d in delta)/width**2)
        cloud*=mass/(cloud.sum()*p.cell_volume)
        rho+=cloud
    return rho


def solve_structural_curvature(density,p=MassGeometry()):
    """Solve periodic potential and exact conformal scalar 3-curvature.

    Gauge choice: mean Phi=0; restricted locked structural field q=-Phi.
    k=0 density is subtracted so homogeneous curvature is not computed.
    """
    rho=np.asarray(density,dtype=float)
    if rho.shape!=(p.grid,)*3 or not np.isfinite(rho).all() or \
          rho.min() < -1e-13 or rho.sum()<=0:
        raise ValueError("finite positive cubic density on this grid")
    k=np.fft.fftfreq(p.grid,d=p.box_length/p.grid)*2*np.pi
    kx,ky,kz=np.meshgrid(k,k,k,indexing="ij")
    ksq=kx*kx+ky*ky+kz*kz
    zero_mask=ksq>0
    density_fluctuation=rho-rho.mean()
    rho_ft=np.fft.fftn(density_fluctuation)
    phi_ft=np.zeros_like(rho_ft)
    phi_ft[zero_mask]=(-4*np.pi*p.G*p.aJ**2 *
                       rho_ft[zero_mask]/ksq[zero_mask])
    # The mean reference potential is fixed to zero.
    q_ft=-phi_ft
    q=np.fft.ifftn(q_ft).real
    gradients=[np.fft.ifftn(1j*kp*q_ft).real for kp in (kx,ky,kz)]
    g2=sum(g*g for g in gradients)
    lap=np.fft.ifftn(-ksq*q_ft).real
    linear_curvature=-4*lap/p.aS**2
    nonlinear_curvature=-np.exp(-2*q)*(4*lap+2*g2)/p.aS**2
    volume_weight=np.exp(3*q)
    vmean=float(np.sum(volume_weight*nonlinear_curvature)/
                np.sum(volume_weight))
    integrated_gradient_mean=float(
        2*np.sum(np.exp(q)*g2)/(p.aS**2*np.sum(volume_weight)))
    variance=float(np.sum(volume_weight *
             (nonlinear_curvature-vmean)**2)/np.sum(volume_weight))
    density_predicted=16*np.pi*p.G*(p.aJ/p.aS)**2 *density_fluctuation
    return dict(
        M_u=float(rho.sum()*p.cell_volume),
        rho_mean=float(rho.mean()),
        rho_peak=float(rho.max()),
        q_max=float(q.max()),q_min=float(q.min()),
        physical_mean_R3=vmean,
        gradient_identity_R3=integrated_gradient_mean,
        maximum_local_R3=float(nonlinear_curvature.max()),
        minimum_local_R3=float(nonlinear_curvature.min()),
        RMS_local_linear_R3=float(np.sqrt(np.mean(linear_curvature**2))),
        linear_coordinate_mean_R3=float(linear_curvature.mean()),
        nonlinear_R3_variance=variance,
        poisson_linear_max_error=float(np.max(
          np.abs(linear_curvature-density_predicted))),
        primordial_K_J_predicted=False,
        numerical_model="post-handoff static periodic Poisson plus exact 3D conformal scalar curvature")


def equal_spaced_centers(p=MassGeometry()):
    L=p.box_length
    return [(i*L,j*L,k*L) for i in (.25,.75)
                        for j in (.25,.75) for k in (.25,.75)]


def redistribution_comparison(p=MassGeometry(), initial_width=.6,
                              after_width=.3):
    initial=solve_structural_curvature(
        periodic_gaussian_density([(p.box_length/2,)*3],
                         [1.],initial_width,p),p)
    after=solve_structural_curvature(
        periodic_gaussian_density(equal_spaced_centers(p),
                         [1/8]*8,after_width,p),p)
    return dict(
        initial_width=float(initial_width),after_width=float(after_width),
        M_u_before=initial["M_u"],M_u_after=after["M_u"],
        before=initial,after=after,
        average_curvature_ratio=float(
            after["physical_mean_R3"]/initial["physical_mean_R3"]),
        average_curvature_reduction_percent=float(100*(
            1-after["physical_mean_R3"]/initial["physical_mean_R3"])),
        maximum_local_curvature_ratio=float(
            after["maximum_local_R3"]/initial["maximum_local_R3"]),
        actual_planck_to_radiation_physics_solved=False)


def diagnostic():
    p=MassGeometry()
    moderate=redistribution_comparison(p,.6,.3)
    very_compact=redistribution_comparison(p,.6,.1)
    return dict(scope="specified post-handoff material redistribution curvature, not derived early SDMC",
       geometry=asdict(p),
       moderate_redistribution=moderate,
       overly_compact_redistribution=very_compact,
       conclusions=[
        "M_u conserved while mean R3 can fall with distributed spatial mass",
        "peak local R3 need not fall when the physical mean does",
        "arbitrarily narrow surviving wells can increase averaged R3",
        "the Poisson zero-mode subtracts homogeneous density and cannot select FLRW K",
        "neither radiation pressure nor scalar/UV source constraints have been solved",
        "curvature dynamics across the prethermal Planck boundary not derived"])
if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
