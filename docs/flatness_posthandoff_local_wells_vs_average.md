# SDMC Part III — Post-handoff local structural curvature and averaged near-flatness

## Scope

This continuation follows the Planck-to-radiation flatness investigation,
not the earlier CDM comparison. The particular hypothesis is that
an initially concentrated configuration of conserved material M_u
redistributes nonuniformly during the structural-to-Jordan handoff.
The resulting geometry contains local structural wells even if its
large-scale background is nearly spatially flat.

Manuscript B Part III §102 gives the *scalar-sector* structural
spatial metric and its linear intrinsic curvature. The equations
below extend those expressions through the full nonlinear
conformal metric identity. This does NOT derive the microscopic
post-handoff mode amplitudes from the Planck boundary.

## 1. Conservation and geometric source

M_u = M_remaining + M_inferred, with M_inferred distributed material,
not a new mass species. Its spatial transport must satisfy the
ordinary conserved-material continuity equation. The material
partition does not determine a unique metric.

The structural spatial metric is

 gamma_ij=M² h_ik (e^{2L})^k_j,
 L^i_j=(ell/3)delta^i_j+Q^i_j.

For the scalar-only conformally flat sector set Q=0 and write

 gamma_ij=a_s² e^{2q_s(x)}delta_ij,
 q_s=m+ell/3-Phi_J,

where m is a local mapper perturbation. The actual physical
post-handoff spatial metric requires m,ell,Phi_J and their
constraints; conservation of M_u alone does not determine them.

## 2. Exact intrinsic spatial Ricci scalar

For a conformal metric in three spatial dimensions,

 R3=-a_s^-2 exp(-2q_s)
             [4 Laplacian(q_s)+2|grad q_s|²].

To first order, delta R3=-4/a_s² Laplacian q_s.
This reproduces Manuscript B Eq (102.5).

For a compact periodic three-dimensional domain, linear coordinate
averaging gives

 <delta R3>_coordinate=0

because the Laplacian is a total derivative.

Crucially that is not the full *physical* volume average.
The proper structural volume is dV_s=a_s³ e^{3q_s}d³x.
Integration by parts yields

 <R3>_physical =
   (2/a_s²) [int exp(q_s)|grad q_s|² d³x]
                  / [int exp(3q_s)d³x] >=0.

A boundary term is present on a finite nonperiodic domain;
with boundaries or traceless Q_ij this special sign result
need not hold. It does NOT imply the global FLRW curvature
index K is positive: the metric is an explicitly flat-reference
periodic conformal example, not the averaged FLRW solution.

For |q|<<1, the leading mean is
  <R3>_physical≈2/a_s² <|grad q_s|²>_coordinate
  =2/a_s² integral k² P_q(k)d³k/(2pi)³.

The leading local linear-curvature variance is
  <(delta R3)^2>≈16/a_s^4
      integral k⁴ P_q(k)d³k/(2pi)³.

These power formulas require appropriate Fourier normalizations
and a finite smoothing/UV cutoff if the short-wavelength
spectrum diverges. Averaging exact curvature and deriving an
effective FLRW backreaction source are *different tasks*.

## 3. Exact solvable single-mode benchmark

Let q_s(x)=epsilon*cos(kx), extended uniformly in y,z on
a periodic three-torus. Then

 R3(x)=a_s^-2 exp[-2epsilon cos(kx)]
       [4 epsilon k² cos(kx)-2epsilon² k² sin²(kx)].

Its *physical volume* mean evaluates exactly to

  <R3>=epsilon² k²/a_s² *
    [I_0(epsilon)-I_2(epsilon)]/I_0(3epsilon),

where I_n is the modified Bessel function of first kind.
For small epsilon:

 <R3>≈epsilon² k²/a_s²,
 RMS(delta R3_linear)=sqrt(8)|epsilon|k²/a_s².

At epsilon=.1, k=1, a_s=1:
  R3 at mode crest =4*.1 exp(-.2)=0.327492...
  R3 at trough     =-4*.1 exp(.2)=-0.488561...
  physical mean R3 =0.0097909662...
  first-order coordinate mean=0,
  linear RMS R3=0.2828427125...

The physical mean is much smaller than the characteristic
local curvature, but it is not exactly zero.
It is suppressed by one extra power of epsilon for
weak wells. This illustrates the coexistence of strongly
varying *local* curvature with small *averaged* curvature.
It does NOT establish an inflation-free solution to the
observed flatness problem or the actual SDMC epsilon.

## 4. Multi-mode well network

For q_s(x)=sum_j epsilon_j cos(k_j x+phi_j)
with nondegenerate harmonics,
the weak second-order mean is

  <R3>≈a_s^-2 sum_j epsilon_j² k_j²,

and the linear curvature variance is

  Var(R3_lin)=8 a_s^-4 sum_j epsilon_j² k_j^4.

The exact result still follows from the proper-volume
identity, not from the naive coordinate average.

Random phase/position changes the local curvature
distribution; a power spectrum controls statistical means
when given. The number of Planck-mass reference units
does not determine epsilon_j or k_j.

## 5. A critical distinction about the original deep well

One "deep primordial gravitational well" refers to a
gravitational potential or material configuration.
It cannot be identified with the background FLRW intrinsic
curvature without solving the Hamiltonian and momentum
constraints. A spatially constant potential has no
Laplacian and no local intrinsic 3-curvature in this
conformal model, although it may affect clocks or
other 4D geometric quantities.

There is no general conservation law
 "R3_primordial = sum_i R3,well_i".
Likewise matter conservation M_u=constant does not
imply conservation of gravitational binding energy.
Radiation pressure can oppose collapse, and the
post-thermal early state may contain perturbations
without mature gravitationally bound wells.

## 6. Exact next piece needed from the Planck-to-radiation dynamics

At the end of the transition t_h, extract

  a_s(t_h), a_J(t_h), h_ij(x,t_h),
  m(x,t_h), ell(x,t_h), Q^i_j(x,t_h)
  and the physical material density/current.

Then use

  gamma_ij(x,t_h)=M²h_ik exp(2L)^k_j

to compute the full spatial Riemann tensor,
the volume mean <R3> and its variance without
assuming conformal Q=0, and the ADM
expansion/shear/energy constraint.

Only with the resulting physically averaged
expansion and 3-curvature can the actual
post-transition Jordan Omega_k be determined
and propagated to today.

The present calculation is an exact geometrical
response test for a chosen structural configuration,
not a prediction of that configuration.

## 7. Implementation

sdmc/flatness_posthandoff_curvature.py
tests/test_flatness_posthandoff_curvature.py

Numerical quadrature validates the exact conformal
identity, periodic integration by parts, modified
Bessel closed form, structural scaling, and
multi-mode weak-spectrum expressions. Tests explicitly
do NOT claim that the homogeneous FLRW curvature
index has been calculated.
