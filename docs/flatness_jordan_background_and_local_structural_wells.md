# SDMC flatness: homogeneous Jordan background and local structural wells

## Meaning of the two geometries

The flatness condition for cosmological background is K_J=0 for
the *homogeneous FLRW background*. The observable Jordan metric
is not required to be locally homogeneous: galaxies, lenses and
radiation perturbations produce local Jordan potentials.

In the candidate independent supporting-space description
(Manuscript B Part III §102), the structural 3-metric is

 gamma_ij=M(x,t)² h_ik (exp(2L))^k_j,
 L^i_j=(ell/3)delta^i_j+Q^i_j.

For scalar perturbations only, choose spatial coordinates
 h_ij=a_J² exp(-2 Phi_J(x,t))delta_ij,
 M(x,t)=Mbar(t) exp(m(x,t)),
 a_s=Mbar a_J.

Then

 gamma_ij=a_s² exp(2q_s)delta_ij,
 q_s=m+ell/3-Phi_J.

It is therefore CONSISTENT as a kinematic ansatz to let
the homogeneous background hbar_ij=a_J²delta_ij be
spatially flat, while gamma_ij has local curvature.

There is a precise locking limitation:

- If Phi_J=0, m=0 but ell!=0, Jordan is exactly locally
  homogeneous while structural gamma has local curvature.
  This is possible geometrically but means L!=0 and so
  has not reached exact independent-relative-mode locking.

- If L=0 and m=0 but Phi_J!=0, both full Jordan h and
  gamma have local perturbations; the Jordan FLRW
  background can remain exactly flat. This is the
  most direct interpretation of physical local gravity
  after homogeneous structural–Jordan locking.

- If L=0, Phi_J=0, m!=0, the structural gamma is
  nonuniform despite homogeneous Jordan, but the local
  mapper field must have an allowed independent
  constrained spatial perturbation. This mechanism
  is not established merely by writing m(x).

- If L=0, m=0, Phi_J=0, gamma=a_s²delta_ij.
  No local structural curvature exists.

Thus the statement "Jordan metric equals global homogeneous
average of local structural metric" is NOT in general an
identity of the manuscript's pointwise metric map. It can
be a candidate **coarse-graining hypothesis**, but its
map and physical averaging dynamics need derivation.

## Exact demonstration that metric averaging is not curvature averaging

Use a periodic structural scalar mode

  q_s(x)=epsilon cos(kx),
  gamma_ij=a_s² exp(2q_s(x))delta_ij.

Then the conformal spatial scalar curvature is

  R3[gamma]= -a_s^-2 exp(-2q_s)
    [4 Laplacian(q_s)+2|grad q_s|²].

Its physical volume mean is

  <R3[gamma]>_V =
   epsilon² k²/a_s²
      [I0(epsilon)-I2(epsilon)]/I0(3epsilon) >0

for nonzero epsilon and periodic boundaries.

But the coordinate average metric is

  <gamma_ij>_coord =
       a_s² I0(2epsilon)delta_ij.

This constant-coefficient spatial metric has

  R3[<gamma>_coord]=0,

even though <R3[gamma]>_V>0.

The operations "average metric, then take curvature"
and "take curvature, then physically average" do not commute.

Numerically for epsilon=.1,k=1,a_s=1:

  R3[<gamma>]=0;
  <R3[gamma]>_V=0.0097909661985093;
  local RMS linear R3=0.282842712474619;
  coordinate mean metric factor I0(.2)=1.010025027795;
  coordinate mean scale factor sqrt[I0(.2)]=1.005000013828;
  proper-volume scale [I0(.3)]^(1/3)=1.007486111408.

Therefore, merely declaring the Jordan metric to be an
average of structural gamma does NOT uniquely give
a scale factor, physical spatial curvature, or a
gravitational evolution equation. The definition
of average matters.

The special positive sign of <R3[gamma]> holds only
for this periodic conformal geometry, and does not
establish a positive global FLRW K. Traceless strain,
boundary conditions, perturbation gauge, backreaction
and extrinsic curvature have not been solved.

## Jordan metric physical constraint

For local Jordan density perturbations, setting
Phi_J=Psi_J=0 everywhere while maintaining arbitrary
delta rho_m generally violates the local gravitational
constraint, unless other physical fields cancel their
source exactly. The fact that K_J=0 on the homogeneous
background is compatible with Phi_J, Psi_J nonzero.

The action-derived structural 3-metric may contain
additional deformations, but observable matter and
photon trajectories depend on the full Jordan
gravitational field. An arbitrary spatial average
of gamma cannot replace its dynamical 4D equations.

In particular the previous Planck-to-radiation
structural suppression 0.0858219 and the averaged
inhomogeneous gamma 3-curvature do not by themselves
predict the actual K_J or the observed Omega_k0.

## Physical flatness interpretation and next test

The strongest retained formulation is:

 "The Jordan FLRW metric describes the large-scale
  homogeneous background; the structural metric
  can carry independent small-scale geometric
  inhomogeneities. The two remain dynamically linked
  by the mapper and relative strain."

The next genuine dynamical test must independently
generate the post-handoff gamma(x), m(x), L(x),
h_J(x), extrinsic geometry and stress from the
Planck-to-radiation action, calculate an action-
consistent large-scale Jordan background, and
show whether its present |Omega_k| is small.

A flat coarse-grained coordinate metric alone is
not an observational flatness solution.
No CDM comparison or source normalization is reopened.

## Test implementation

sdmc/flatness_jordan_background_vs_structural_wells.py
tests/test_flatness_jordan_background_vs_structural_wells.py

The code constructs exact geometric examples for
homogeneous Jordan / inhomogeneous supporting metric,
exact locked local Jordan perturbations, mapper-only
variations, and fully homogeneous exact lock.
It compares the averaged structural metric against
the physically averaged local 3-curvature and tests
both candidate coarse-grained spatial rulers.
