# SDMC Manuscript B — Part III flatness continuation: structural curvature immediately after material redistribution

## 1. Purpose and scope

This study retains the completed second-route structural-gravity
interpretation and addresses only the *flatness investigation*.
The conceptual step is to calculate the intrinsic structural
3-curvature of a nonhomogeneous post-handoff material distribution,
rather than assuming one deep primordial gravitational depression
simply converts into exactly the same sum of smaller curvatures.

We distinguish:
  (i)  local geometric curvature R_s(x,t_h),
  (ii) volume-averaged structural scalar curvature <R_s>_D,
  (iii) structural curvature fluctuation variance and
  (iv) physical homogeneous Jordan FLRW curvature Omega_k,J.
These four quantities are NOT interchangeable.

Manuscript B §100 defines actual redistributed material
M_inferred, not a new source mass, with conserved inventory
  M_u=M_remaining+M_inferred.
The sample post-handoff matter distributions conserve M_u
exactly while changing the spatial arrangement.

The exact prethermal Planck boundary S_P=1 has
rho_v,P=rho_P and rho_X,P=0; there is no assumed classical
Jordan matter Poisson equation at that boundary.
The following calculations are conditioned on a meaningful
post-emergence structural spatial slice t_h.

## 2. Exact nonlinear structural scalar curvature

Manuscript B §102:
  gamma_ij=M² h_ik(exp(2L))^k_j.

In the conformally flat scalar-only restriction, neglecting
traceless shape and allowing an independent local mapper and
trace deformation,

  gamma_ij=a_s² exp(2q_s) delta_ij,
  q_s=m+ell/3-Phi_J,
  a_s=S/S0.

The exact scalar 3-curvature is

  R_s^(3)(x)=
   -(exp(-2q_s)/a_s²)[4 Δq_s+2|∇q_s|²].

At linear order R_s=-4 Δq_s/a_s²,
matching Manuscript B (102.5). The exact nonlinear expression
includes gradient-square and spatial volume-weighting terms.

For a periodic or compact boundaryless domain,
  dV_s=a_s³ exp(3q_s)d³x,
  <R_s>_V=∫R_s dV_s/∫dV_s.

Using integration by parts,

  <R_s>_V=
     (2/a_s²)[∫ exp(q_s)|∇q_s|² d³x]
                  /[∫exp(3q_s)d³x] >= 0.

At first order, <R_s>_coordinate=0 under periodic boundaries.
At next-to-leading order,
  <R_s>_V≈ 2/a_s² <|∇q_s|²>_coordinate >0.

This effect is NOT a new dark fluid; it is the nonlinear
geometrical average of the specified metric. The physical
mean generally does not vanish exactly even though the
linear mean does. This second-order mean cannot be equated
to the global background curvature K_FLRW without the
Hamiltonian/extrinsic-curvature and averaging dynamics.

## 3. Curved reference background: avoid assuming the answer

A more general conformal 3-metric is

  gamma_ij=a_s² exp(2q_s) bar_gamma_ij,

with reference curvature
  R[bar_gamma]=6 K_R.

Then

  R_s=a_s^-2 exp(-2q_s)
     [6K_R-4 bar_Delta q_s-2|bar_grad q_s|²].

For a compact boundaryless reference domain,

  <R_s>_V =
    [6 K_R ∫exp(q_s) dV_bar
      +2∫exp(q_s)|bar_grad q_s|² dV_bar]
    /[a_s²∫exp(3q_s)dV_bar].

The conformal-gradient term is nonnegative. For an initially
POSITIVE K_R, the mean retains a positive contribution and
cannot vanish exactly within this restricted class.
For negative K_R, partial cancellation is possible but is
not guaranteed or dynamically selected. On a finite
nonperiodic domain a boundary term also enters.

Thus converting a single gravitational depression into
many localized deformations need not erase a pre-existing
homogeneous FLRW curvature parameter, and the flat
periodic illustrative simulations below assume K_R=0
as a reference rather than predict that selection.

## 4. Exact harmonic test at fixed material density-source amplitude

Take a single spatially varying conformal structural mode

  q_s(x)=A cos(k x), L=2pi,

representing curvature concentrations varying along x
(they are not independent bound 3D galaxies).

At physical-volume weight exp(3q), the exact mean is

  <R_s>_V = A² k²/a_s² *
       [I0(A)-I2(A)]/I0(3A),

with modified Bessel functions I_n.

At leading order, the local curvature RMS is
  R_rms≈2sqrt(2) A k²/a_s².

Suppose a fixed source density contrast is related to q
by an illustrative weak-field Poisson equation
  Δq=-alpha_delta delta_m.
For fixed delta amplitude, A*k² is constant, not A.
Thus set A=0.02/k², a_s=1. Then

 k=1:  mean R=3.996602e-4; RMS~.05655
 k=2:  mean R=9.999469e-5; RMS~.05657
 k=4:  mean R=2.499992e-5; RMS~.05657
 k=8:  mean R=6.249999e-6; RMS~.05657
 k=16: mean R=1.562500e-6; RMS~.05657.

At constant density contrast and smaller length scales,
the mean curvature decreases about k^-2 while local
curvature fluctuations remain fixed.

However if A (the metric deformation amplitude) is held
fixed rather than the density source, the mean grows ~k².
Thus "more localized wells" alone is NOT a universal
flattening mechanism: the amplitude–width scaling law matters.

## 5. Three-dimensional conserved-mass redistribution experiment

We also construct genuinely localized Gaussian material
concentrations in a three-dimensional PERIODIC comoving
box of side L=2pi.

At each n in {8,32,128}, we distribute the SAME total
material mass M_u among n equal-mass concentrations.
Each normalized Gaussian satisfies
  <rho_i>_coordinate=mass fraction m_i/M_u,
so <rho_total>=1, exactly up to floating-point precision.

To keep their nominal total occupied volume approximately
constant while increasing n, set
  sigma_n=.35*(8/n)^(1/3).

Position centers are randomly sampled. Six deterministic
seeds {1,2,3,4,5,6} provide a modest ensemble.

For a subhorizon, weak-field Jordan-matter-inspired model
after emergence, choose an ILLUSTRATIVE relation
  Δq_s=-alpha_delta delta_m,
  delta_m=rho_total/<rho_total>-1,
  alpha_delta=.002.

This model keeps m=ell=0 in the relation
q_s=m+ell/3-Phi_J and uses q_s≈-Phi_J. It is NOT
the full relativistic photon-baryon constraint during
the Planck-to-radiation handoff, nor a calculation of
surviving independent structural response.

Solve the mean-zero periodic Poisson equation by FFT,
  q_k=alpha_delta delta_k/k², for nonzero Fourier k,
and then calculate exact intrinsic curvature of
gamma_ij=a_s² exp(2q_s)delta_ij on a 64³ grid.

Average over six seeds in each n:

 n=8:
   <R_s>_V ≈2.17578e-5,
   local curvature RMS ≈0.031535,
 n=32:
   <R_s>_V ≈9.96248e-6,
   local curvature RMS ≈0.031999,
 n=128:
   <R_s>_V ≈4.47988e-6,
   local curvature RMS ≈0.032318.

The mean has decreased by nearly a factor 4.86 between
8 and 128 concentrations, while local fluctuation RMS
is nearly unchanged. Total matter mass is exactly
conserved. Each snapshot independently satisfies
the Poisson relation and the exact positive
curvature averaging identity to numerical accuracy.

This numerically REALIZES the *conditional geometrical*
part of the initial analogy, but does not show that the
actual stiff-to-radiation microphysics realizes the
chosen mass profiles, Poisson source, or density contrasts.

## 6. Consequences for actual near-flatness

The observed cosmic near-flatness is a Jordan-space
geometric property. Neither a vanishing linear curvature
average nor a small nonlinear averaged R_s directly
fixes the homogeneous K_J or Omega_k,J of an inhomogeneous
realization. Averaged dynamics require expansion/shear
and other gravitational constraints; in SDMC the structural
3-metric does not yet define a unique full 4D structural
spacetime metric.

The result relevant to the current flatness investigation is:

  Material redistribution with an appropriate
  amplitude–length scaling CAN make <R_s>_V small
  compared with local curvature fluctuations.

It is mathematically possible for the post-handoff
spatial geometry to contain many local curvature
concentrations while having a small average curvature.

It is NOT established that a previously order-unity
homogeneous k_R is transformed to a small K_J during
the Planck-to-radiation transition. The Gaussian/periodic
tests began with a flat reference geometry K_R=0 and
used externally specified distributions. No fragmentation
conservation law exists for the scalar curvature,
curvature-squared action energy or gravitational binding
energy akin to the conserved M_u partition.

The physical next step is to solve the actual spatial
Hamiltonian + relative-metric equations ON A NONFLAT
REFERENCE background through the thermal handoff, derive
its initial spatial geometry and post-handoff matter
power spectrum without assuming their amplitude, and
then evaluate <R_s>, variance, and the Jordan curvature
in the same solution.

## 7. Reproducibility

  sdmc/flatness_posthandoff_curvature.py
  tests/test_flatness_posthandoff_curvature.py
  .github/workflows/sdmc_flatness_posthandoff_curvature.yml

The workflow verifies mass conservation, source
normalization, Fourier-Poisson residual, exact nonlinear
curvature identity, the Bessel analytic one-mode average,
the 3D randomized localization comparison, and that
background curvature and observable Jordan Omega remain
UNDETERMINED in these illustrative snapshots.

This work does not reopen the CDM/second-route investigation.
