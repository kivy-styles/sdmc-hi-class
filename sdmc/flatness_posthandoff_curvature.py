"""Post-handoff structural 3-curvature from conserved M_u redistribution.

THIS IS A GEOMETRIC POST-HANDOFF SNAPSHOT, NOT A PLANCK-TO-RADIATION
DYNAMICAL SIMULATION OR AN INDEPENDENT SDMC CURVATURE PREDICTION.

Manuscript B, Part III, §102:
    gamma_ij=M² h_ik(e^{2L})^k_j.
In scalar-only conformally flat restriction:
    gamma_ij=a_s² exp(2q_s) delta_ij,
    q_s=m+ell/3-Phi_J.
Its exact intrinsic spatial Ricci scalar is
    R_s=-(exp(-2q_s)/a_s²)(4 ∇² q_s+2|∇q_s|²).
The physical volume weight is a_s³ exp(3q_s).

For periodic boundaries the exact physical spatial average is
    <R_s>_V = 2/a_s² [<exp(q_s)|∇q_s|²>_0
                        /<exp(3q_s)>_0] >= 0.

The first-order coordinate-average of -4 ∇²q/a_s² is 0, but
the physical-volume average is generally positive at O(q²).
This is NOT the homogeneous FLRW curvature constant K, and
neither the volume mean nor its variance is automatically
an observable effective Omega_k.

For a curved reference metric gamma=a_s² exp(2q) bar_gamma:
    R_s=a_s^-2 exp(-2q)[6K_R-4 bar_Δq-2|bar_∇q|²].
For a compact, boundaryless constant-curvature reference,
    <R_s>_V= a_s^-2
      [6K_R∫e^q +2∫e^q|bar_∇q|²]/∫e^{3q},
with reference measure understood. Thus starting with
K_R=0 is an ASSUMPTION, not a result of fragmenting M_u;
for K_R>0 in this conformal class the mean stays >0.

The optional weak-field Poisson source used below is
    ∇² q_s = -alpha_delta * delta_m,
a subhorizon matter-motivated toy when q_s≈-Phi_J and
mapper/independent trace are ignored. At the real early
radiation transition the fully relativistic Jordan+structural
constraint must replace this source. In particular neither
alpha_delta nor q_s's primordial spectrum is fixed here.

The total material inventory M_u is conserved in all blob
examples. M_inferred is actual redistributed material
from the retained core, not a newly created particle species.
"""
from __future__ import annotations
import json
import numpy as np
from scipy.special import iv

def exact_conformal_curvature_1d(q,grad,lap,a_s=1.):
    """Three-dimensional R_s for a field depending on one periodic x."""
    q=np.asarray(q)
    grad=np.asarray(grad)
    lap=np.asarray(lap)
    if not (q.shape==grad.shape==lap.shape) or a_s<=0:
        raise ValueError("matching finite fields and positive a_s needed")
    R=-np.exp(-2*q)*(4*lap+2*grad*grad)/a_s**2
    vol=np.exp(3*q)
    physmean=np.mean(R*vol)/np.mean(vol)
    identity=2*np.mean(np.exp(q)*grad*grad)/(a_s*a_s*np.mean(vol))
    variance=np.mean((R-physmean)**2*vol)/np.mean(vol)
    return dict(mean=float(physmean),positive_identity=float(identity),
                rms=float(np.sqrt(variance)),
                coordinate_linear_mean=float(np.mean(-4*lap/a_s**2)),
                pointwise_min=float(R.min()),
                pointwise_max=float(R.max()))

def harmonic_example(k=4,C=.02,a_s=1.,resolution=8192):
    """Fixed curvature / density-source mode amplitude C=A*k²."""
    if int(k)!=k or k<1 or resolution<32 or k>=resolution//4:
        raise ValueError("resolved nonzero integral k required")
    if C<=0 or a_s<=0:raise ValueError("positive source and ruler")
    A=C/k**2
    x=2*np.pi*np.arange(resolution)/resolution
    q=A*np.cos(k*x)
    grad=-A*k*np.sin(k*x)
    lap=-A*k*k*np.cos(k*x)
    result=exact_conformal_curvature_1d(q,grad,lap,a_s)
    bessel=A*A*k*k*(iv(0,A)-iv(2,A))/iv(0,3*A)/a_s**2
    return dict(k=int(k),A=float(A),source_amplitude=float(C),
                mean_bessel=float(bessel),**result)

def conserved_gaussian_material(positions,sigma,N=64,
                               L=2*np.pi,core_fraction=0.):
    """Construct a normalized 3D periodic matter-density field.

    Positions are distributed material centers; each gets equal
    fraction (1-core_fraction)/n of the total conserved M_u.
    A retained core at (L/2,L/2,L/2), when requested, is an
    additional material portion of the SAME total mass.
    The density is normalized so mean rho=1 and δmean=0.
    """
    pos=np.asarray(positions,dtype=float)
    if pos.ndim!=2 or pos.shape[1]!=3 or len(pos)<1:
        raise ValueError("n by 3 material coordinates expected")
    if not np.isfinite(pos).all() or not (0<=core_fraction<1) or \
       N<16 or sigma<=0 or L<=0:
        raise ValueError("regular periodic source required")
    x=np.arange(N)*L/N
    densities=np.zeros((N,N,N),dtype=float)

    def add_blob(position,weight):
        factors=[]
        for axis in range(3):
            dx=(x-position[axis]+L/2)%L-L/2
            factors.append(np.exp(-.5*(dx/sigma)**2))
        gauss=(factors[0][:,None,None]*
            factors[1][None,:,None]*factors[2][None,None,:])
        densities[:,:,:]+=weight*gauss/gauss.mean()

    for p in pos:
        add_blob(p,(1-core_fraction)/len(pos))
    if core_fraction:
        add_blob(np.full(3,L/2),core_fraction)
    return densities

def structural_snapshot_from_delta(delta,alpha_delta=.002,
                                   a_s=1.,L=2*np.pi):
    """Conserved periodic δ -> q by weak-field Poisson, then exact R.

    All k are Fourier COMOVING wavenumbers, not homogeneous K_FLRW.
    The source is illustrative only. Physical scalar curvature of
    its specified conformally flat metric is nevertheless exact.
    """
    delta=np.asarray(delta,dtype=float)
    if delta.ndim!=3 or delta.shape[0]<16 or \
       any(n!=delta.shape[0] for n in delta.shape) or \
       not np.isfinite(delta).all():
        raise ValueError("finite cubic 3D density-contrast grid expected")
    if alpha_delta<=0 or a_s<=0 or L<=0:
        raise ValueError("positive spatial ruler and coupling expected")
    # The mean density is conserved and its contrast must have zero
    # Fourier zero mode; eliminate floating point roundoff only.
    d=delta-delta.mean()
    n=delta.shape[0]
    kn=np.fft.fftfreq(n,d=L/n)*2*np.pi
    kx,ky,kz=np.meshgrid(kn,kn,kn,indexing="ij")
    k2=kx*kx+ky*ky+kz*kz
    d_k=np.fft.fftn(d)
    q_k=np.zeros(d_k.shape,dtype=complex)
    q_k[k2>0]=alpha_delta*d_k[k2>0]/k2[k2>0]
    q=np.fft.ifftn(q_k).real
    dq=(np.fft.ifftn(1j*kx*q_k).real,
        np.fft.ifftn(1j*ky*q_k).real,
        np.fft.ifftn(1j*kz*q_k).real)
    lap=np.fft.ifftn(-k2*q_k).real
    g2=sum(a*a for a in dq)
    R=-np.exp(-2*q)*(4*lap+2*g2)/a_s**2
    dV=np.exp(3*q)
    mean=np.mean(R*dV)/np.mean(dV)
    exact_integrated_identity=(2*np.mean(np.exp(q)*g2)
                                 /(a_s**2*np.mean(dV)))
    variance=np.mean((R-mean)**2*dV)/np.mean(dV)
    raw_linear_mean=np.mean(-4*lap/a_s**2)
    return dict(
      sources="conditional Newtonian Poisson source; exact curvature of specified post-handoff structural metric",
      mean_curvature=float(mean),
      integrated_identity=float(exact_integrated_identity),
      rms_curvature=float(np.sqrt(variance)),
      curvature_rms_about_zero=float(np.sqrt(np.mean(R*R*dV)/np.mean(dV))),
      min_curvature=float(R.min()),max_curvature=float(R.max()),
      max_potential=float(q.max()),min_potential=float(q.min()),
      coordinate_mean_linear_curvature=float(raw_linear_mean),
      max_poisson_residual=float(np.max(abs(lap+alpha_delta*d))),
      mean_density_contrast=float(d.mean()),
      density_contrast_rms=float(np.sqrt(np.mean(d*d))),
      dimensions="dimensionless curvature times chosen coordinate length squared",
      background_K_R_determined=False,
      jordan_Omega_k_determined=False)

def material_fragmentation_case(n=8,seed=7,grid=64,
                                sigma_at_8=.35,alpha_delta=.002):
    if n<1 or sigma_at_8<=0:raise ValueError("n>=1")
    rng=np.random.default_rng(seed)
    L=2*np.pi
    x=rng.uniform(0,L,(n,3))
    sigma=sigma_at_8*(8/n)**(1/3)
    density=conserved_gaussian_material(x,sigma,N=grid,L=L)
    r=structural_snapshot_from_delta(
        density-1,alpha_delta=alpha_delta,L=L)
    return dict(n=int(n),seed=int(seed),
       individual_mass_fraction=float(1/n),
       redistributed_fraction=1.,
       total_mass_fraction=1.,
       sigma=float(sigma),
       mean_density=float(density.mean()),
       **r)

def ensemble_report(ns=(8,32,128),seeds=range(1,7),grid=64):
    groups=[]
    for n in ns:
        samples=[material_fragmentation_case(n=n,seed=s,grid=grid)
                 for s in seeds]
        groups.append(dict(n=int(n),sample_count=len(samples),
             mean_structural_curvature=float(np.mean(
                  [v["mean_curvature"] for v in samples])),
             mean_local_rms_curvature=float(np.mean(
                  [v["rms_curvature"] for v in samples])),
             scatter_structural_mean=float(np.std(
                  [v["mean_curvature"] for v in samples])),
             max_density_mass_error=float(max(
                  abs(v["mean_density"]-1) for v in samples)),
             max_integrated_identity_residual=float(max(
                  abs(v["mean_curvature"]-v["integrated_identity"])
                  for v in samples))))
    return dict(scope="conditional post-handoff conserved matter snapshots on flat conformal reference, NOT primordial K prediction",
       groups=groups,loss_of_mass=False,
       main_limitation="The homogeneous structural spatial curvature and full early thermal dynamics are not specified.")

def diagnostic():
    return dict(one_mode=[harmonic_example(k) for k in (1,2,4,8,16)],
       gaussian_ensemble=ensemble_report(),
       no_flatness_claim="Zero reference K_R was imposed in conformal post-handoff snapshot; does not compute origin of K_R or K_J.")

if __name__=="__main__":
    print(json.dumps(diagnostic(),indent=2))
