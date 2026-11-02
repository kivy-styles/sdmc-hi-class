# SDMC flatness — Post-handoff curvature from an explicitly redistributed material inventory

## Status and objective

The question is whether the user's analogy —
an initially concentrated M_u configuration forming many distinct,
nonuniform gravitational wells during Planck-to-radiation
structural–Jordan handoff — can be *used to calculate*
the structural spatial curvature immediately afterward, and
whether its large-scale mean can decrease while local wells persist.

This continuation does not reopen CDM. It also makes no
claim that an actual initial Planck boundary had one classical
material "lump"; that initial condition requires a physical
primordial-material origin and matching calculation.

Part III §102 gives the structural–Jordan spatial metric:

  gamma_ij=M² h_ik (exp(2L))^k_j,
  L^i_j=(ell/3) delta^i_j+Q^i_j.

For a restricted scalar mode with Q=0,

  gamma_ij=a_s² exp(2q_s)delta_ij,
  q_s=m+ell/3-Phi_J.

Here m=delta ln M is a perturbation of the mapper, not
a Planck or ordinary material mass.

## 1. The direct local curvature equation

The exact three-dimensional intrinsic scalar curvature is

 R3_s=-exp(-2q_s)/a_s²
        [4∇²q_s+2(∇q_s)²].

At first order,

 delta R3_s=-4∇²q_s/a_s².

For the baseline locked scalar case m=ell=Q=0, let q_s=-Phi_J.
If post-handoff matter is nonrelativistic, weakly gravitating,
subhorizon, and the scalar/metric constraints reduce to
Newtonian Poisson, then

 ∇² Phi_J=4pi G a_J² (rho_m-bar rho_m).

Therefore

  delta R3_s=16pi G (a_J/a_s)² (rho_m-bar rho_m).

This is an action-compatible kinematic *test baseline*,
not a full SDMC early radiation or UV closure. The exact
prethermal Planck boundary does not admit the assumed
classical Poisson equation.

## 2. Conserved Gaussian fragment construction

Use a cubic periodic spatial domain of side 2pi and grid 64³,
normalized total M_u=1, G=.05, a_J=a_s=1 in arbitrary
dimensionless units.

The density is a sum of normalized periodic Gaussian lumps
such that the integral of each is exactly m_i,
and the sum remains M_u=1 in every configuration.

Before: 1 Gaussian concentration of mass 1 at box center,
comoving width sigma=.6.

After (moderately concentrated): 8 Gaussian concentrations
of mass 1/8 each at 8 spatially separated octants,
width sigma=.3 each.

The density contrast rho-bar rho has zero volume mean.
Solve its Poisson potential by Fourier inversion with the
zero mode removed. Compute q_s=-Phi_J and R3_s using spectral
spatial derivatives. Average with the physical metric volume
dV_s=a_s³ exp(3q_s)d³x.

The exact periodic integration-by-parts identity is

 <R3_s>_V=(2/a_s²)
    [int exp(q_s)|∇q_s|² d³x]/
    [int exp(3q_s)d³x]>=0.

The comparison is at fixed a_s and a_J to isolate the
mass-redistribution effect. Actual structural expansion
would provide an additional a_s^-2 factor for fixed q,
but the sourced q must be solved as the universe expands.
The number of concentrations does not equal the number
of Planck masses or necessarily gravitational wells.

## 3. Numerical results

Before (one concentration, width .6):
 M_u=1
 physical mean R3_s=0.00013018751559455962
 largest positive local R3_s=0.6667786749516493
 RMS linear local R3_s=0.0504365039231869
 max |q_s|≈0.044368

After (eight equally spaced concentrations, width .3):
 M_u=1
 physical mean R3_s=0.00003221578677563654
 largest positive local R3_s=0.7126635508344074
 RMS linear local R3_s=0.05043647647414826
 max |q_s|≈0.011092

Hence

 mean-curvature ratio=0.24745680588886512
 reduction in positive mean R3=75.2543194111%
 maximum local R3 after / before≈1.06882.

The modeled spatial mean is smaller even though the
maximum local intrinsic scalar curvature is slightly
larger: localized strong wells and small average
curvature coexist at fixed conserved M_u.

Linear coordinate-averaged R3=0 by periodicity,
but the exact proper-volume mean is positive and
thus needs the second-order calculation.

## 4. Counterexample: overly narrow fragments

Keep precisely the same eight centers and total M_u=1,
but use width sigma=.1 for each after concentration.

Then the numerics give
 physical mean R3_s≈0.0001525595160064761,
 largest positive local R3_s≈18.2487125.

The physical average is now approximately 17.18%
LARGER than the one-lump initial state. Compact cores
may preserve/increase gradient energy while reducing
long-wavelength density coherence. The result is
sensitive to the spatial density profile and scale.

Conclusion: fragmentation itself is NOT an automatic
flatness-selection or curvature-suppression law.

## 5. Why average spatial curvature is not observable FLRW K

The above proper-volume average is the intrinsic
3-curvature of the selected *structural spatial slice*.
It is not necessarily the global FLRW constant-curvature
index K_J nor the present Omega_k,J.

In GR's ADM convention (unit normal, c=1),

  R3+K_ADM²-K_ij K^ij=16piG rho+2Lambda.

With K_ij=-H_loc gamma_ij-sigma_ij,
  K_ADM²-K_ij K^ij=6H_loc²-sigma_ij sigma^ij,

so

  R3=16piG rho+2Lambda-6H_loc²+sigma_ij sigma^ij.

This shows why spatial density alone does not
determine the full gravitational geometry: the
extrinsic curvature, shear, and lapse/shift/matter
constraints are essential.

In actual SDMC scalar–tensor gravity the Hamiltonian
constraint further includes the nonminimal coupling F,
braiding G3, structural spatial fields and mapper
gradients. None has been solved in this simple
Poisson model.

For periodic conformal q, a positive second-order
mean R3 can arise from local perturbations even while
the reference background is globally flat. It must not
be relabeled an FLRW K or an independently conserved
gravitational-energy inventory.

## 6. What the actual Planck-to-radiation solution must supply

At the first classical Jordan-compatible slice after
the prethermal handoff, derive rather than assign:

 - the physically selected material concentrations or
   stochastic spatial spectrum;
 - radiation and structural kinetic energy densities and
   energy-transfer rate J from the actual microphysical
   action;
 - lapse N, mapper M, structural tracking p and
   boundary-matched scalar coefficient functions;
 - the spatial relative strain ell and Q_ij plus its
   sources and initial data;
 - the Jordan metric potentials from the full ADM
   Hamiltonian and momentum constraints;
 - the 3D structural spatial metric gamma and full
   evolution of intrinsic and extrinsic curvature.

Once these are independently available,
   R3_s(x,t_h), <R3_s>_physical, Var(R3_s),
and the Jordan post-handoff K_J / Omega_k,J can be
calculated on the same physical initial data.

The mass redistribution clue is therefore supported
as a mathematically possible *smoothing of the mean
structural curvature with surviving local wells*.
It does not, without a transition-produced density
distribution and full constraints, prove the
mechanism dynamically solved observational flatness.

## 7. Implementation

Module: sdmc/flatness_mass_poisson_curvature.py
Tests: tests/test_flatness_mass_poisson_curvature.py

Tests explicitly check exact mass conservation,
spectral linear Ricci-to-density identity,
periodic vanishing of linear average,
nonlinear proper-volume gradient identity,
a strong-local/weak-average configuration,
and a more compact distribution that
INCREASES rather than decreases the mean.


## 8. Remove the flat-reference assumption: exact non-flat S³ slice

The periodic Poisson calculation assumed a flat reference
Euclidean three-torus. This is appropriate for computing
post-handoff *local inhomogeneous curvature* but cannot
by itself establish how an originally non-flat universe
became nearly flat.

We therefore performed a second, independent, exact
geometric calculation with an explicitly positively
curved three-dimensional reference metric:

  dbar s² = dchi²+sin²chi dOmega_2²,
  Rbar3=+6,
  0<=chi<=pi.

Define the structural metric

  gamma_ij=a_s² exp(2q(chi)) bar h_ij,
  q(chi)=eps coschi.

This is not a primordial SDMC initial-condition prediction;
it is a non-flat metric family with explicit finite
spatial deformation. The resulting exact curvature is

  R3_s(chi)=exp(-2eps coschi)/a_s² *
   [6+12eps coschi-2eps² sin²chi].

Correct physical averaging using
  dV_s=a_s³ exp(3q) sin²chi dchi dOmega_2
gives

  <R3_s>_V =
      (1/a_s²)
      [6 int exp(q) dVbar
          + 2 int exp(q)|grad_bar q|² dVbar]
       / [int exp(3q) dVbar] > 0.

The gradient identity is obtained from integration
by parts on the closed S³ manifold. Local R3 can
change sign even while its physical mean stays
strictly positive.

For a_s=1, calculate:

  eps=0: mean=6.000000, 0% reduction;
  eps=.1: mean=5.955267, 0.75% reduction;
  eps=.5: mean=5.021311, 16.31% reduction;
  eps=1: mean=3.191283, 46.81% reduction;
  eps=2: mean=.871078, 85.48% reduction;
  eps=3: mean=.186632, 96.89% reduction;
  eps=4: mean=.035172, 99.41% reduction;
  eps=5: mean=.006137, 99.8977% reduction.

This proves a **geometric existence** statement:
a non-flat spatial manifold can have an extremely
small physical *average intrinsic scalar curvature*
compared with its original homogeneous reference
while retaining a non-flat topology and high local
curvature. It does not mean the effective homogeneous
FLRW curvature label K has changed or that observable
Omega_k has been made acceptably small.

The substantial suppression at eps~3–5 requires
large spatial metric variations e^(2q), beyond the
small-perturbation regime. No action-derived dynamics,
physical matter distribution, radiation transport,
stability or global averaging constraints have been
provided that generate these profiles.

This is nevertheless an important control:
we are no longer relying entirely on an initially flat
coordinate background to check the mathematical
possibility of small averaged curvature.

The two tests must not be combined as one actual
transition:
 (i) conserved M_u redistributed through a periodic
     Newtonian matter source on a flat reference;
 (ii) non-flat positively curved S³ geometry with
     prescribed structural deformation.

The missing step is to solve the *same* physical
Planck-to-radiation equations to derive both the
material redistribution and the spatial deformation
from common initial conditions, then evaluate its
Jordan FLRW distance/expansion observables.

Module: sdmc/flatness_nonflat_structural_slice.py
Tests: tests/test_flatness_nonflat_structural_slice.py
