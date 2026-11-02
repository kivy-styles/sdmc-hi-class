# SDMC Manuscript B — Flatness: probability of near-flatness after material redistribution

## Purpose and limits

Question: Does the Planck-to-radiation handoff make near-flatness a
likely rather than fantastically improbable result of initial material
configurations? Can we calculate the *fraction of physically allowed*
initial configurations that become observationally near-flat?

The physical answer cannot yet be given as a unique number: the
current SDMC primordial action/transition has no derived normalized
probability measure over the admissible prethermal configurations,
no distribution of their material density correlations, and no
self-consistent mapping of each such initial state through the
nonflat structural-to-Jordan ADM evolution to the *observable*
present Jordan Omega_k. Conservation of M_u and positivity of matter
density constrain configurations but do not assign relative chances.

This study therefore reports two SEPARATE conditional experiments.
It NEVER equates suppression of averaged structural Ricci curvature
with the observational Jordan FLRW curvature fraction.

The historical estimate "initial density adjusted to 1 part
in 10^60" is a fine-tuning/SENSITIVITY argument in a particular
classical standard cosmological parameterization and time, not
itself a probability of 10^-60 unless a justified measure is
introduced. The chosen parameterization and measure matter.

For an observational scale, Planck+BAO gave
 Omega_k,0 = 0.0007 +/- 0.0019 (68% CL)
(Planck 2018 results VI, A&A 641 A6 2020).
Use |Omega_k,0|<0.004 as a simple approximate near-flatness
benchmark, not a strict observational probability statement.
An additional tighter target |Omega_k,0|<0.001 is reported.

## Experiment I — Conditional observable Jordan curvature

Use the published manuscript conditional FLAT-BACKGROUND reference
propagation law (not a self-consistent nonflat solution):

  |Omega_k,J,0|=A |C_k,0 Omega_k,R,P|,
  A=0.0917953.

Take C_k,0=1 in a common fixed-curvature FLRW coordinate
normalization, and assume by fiat

  X = Omega_k,R,P ~ Uniform(-B,B),

where X is the provisional *structural* Planck curvature diagnostic.
It is NOT necessarily conventional primordial Jordan Omega_k,
and uniformity over it is not derived from the SDMC Planck boundary.

Then the measure fraction is exactly

  P(|Omega_k,J,0|<epsilon | Uniform(-B,B), Ck=1)
       =min(1, epsilon/(A B)).

For epsilon=.004:
  B=1:    P=0.04357521572 (4.3575%);
  B=.1:   P=0.43575215725 (43.5752%);
  B=.01:  P=1 (100%).

For epsilon=.001:
  B=1:    P=0.01089380393 (1.0894%);
  B=.1:   P=0.10893803931 (10.8938%);
  B=.01:  P=1 (100%).

Changing the assumed prior range alone changes "probability"
by almost two orders of magnitude. These are not SDMC's
predicted probabilities. A prior that concentrates near X=0
naturally yields a much higher fraction but would need a
microscopic justification, rather than being assumed.

The present code leaves Ck as an explicit constant, not an
arbitrary evolving time-dependent curvature ratio. A
time-only spatial mapper cannot arbitrarily change Ck when
both slices have fixed homogeneous FLRW curvature labels.

## Experiment II — Material configurations and structural averaging

We also calculate a fraction over positive, finite,
mass-conserving post-handoff matter profiles, keeping
the same M_u=1 in every realization.

This is a deliberately restricted *plane-symmetric*
periodic 3-D geometric model (q=q(x)), rather than
eight truly separated 3D spherical wells. The x-domain
is periodic with L=2pi; density is uniform in y,z and
inhomogeneous in x. The initial source is a unit mass
Gaussian sheet at a fixed location of width sigma0=.6.
After redistribution, eight equal-mass Gaussian sheets
each contain 1/8 of M_u.

The normalized nonzero Fourier coefficients are

  rho_k ∝ exp(-k²sigma²/2) *
             (1/8)sum_i exp(-ik x_i).

In the weak-field subhorizon Jordan Poisson approximation,
with no independent mapper or L perturbation,
q=-Phi_J, and the leading second-order proper-volume
average structural spatial Ricci curvature is

  <R3_s> ≈ (2/a_s²)<(grad q)²>
          ∝ sum_{k>=1}|rho_k|²/k².

The ratio to the prescribed initial single broad Gaussian is

  S(x_i,sigma)=
      sum_{k=1}^100
         exp(-k²sigma²) |(1/8)sum_i exp(-ikx_i)|² / k²
      -------------------------------------------------
      sum_{k=1}^100 exp(-k²sigma0²)/k².

The normalization is independent of total mass and G since
the same M_u, a_s, a_J and linear Poisson coefficient appear
on both sides. This ratio is the weak-geometry analog of
earlier exact structural-curvature calculations, NOT
present observable Omega_k,J.

The sample criterion S<0.25 means 75% *mean structural*
curvature suppression. S<0.10 means 90%; S<0.01 means 99%.

The Monte Carlo makes completely explicit, arbitrary
statistical priors on x_i and sigma:

  A. Independent uniform positions: all eight sheet
     positions independently uniform on [0,2pi).

  B. One clustered group: a uniform common center, and
     positions displaced by independent N(0,.2²) offsets.
     Marginal positions remain uniform, but pair
     correlations enhance collective large-scale density.

  C. Almost evenly spaced: eight sheet positions
     around 2pi*i/8 with random collective translation
     and independent N(0,.1²) jitter.

In every model, the single width sigma per configuration
is independently Uniform(.1,.6). The initial width=.6
is held fixed. Every distributed material mass equals1/8.
Thus total M_u=1 exactly in all configurations.

The numerical experiment samples 20,000 independent
configurations per prior, using a fixed pseudorandom seed
4102 for reproducibility. It does not filter by a full
SDMC Hamiltonian constraint, thermal birth of matter,
genuine 3D strain, radiation pressure, N,p or stability.
These are "allowed" only in the restricted positive-density
Poisson/conformal toy, not certified physically allowed
Planck-to-radiation initial conditions.

### Numerical Monte Carlo results

                                          probability of specified
                                          structural mean reduction
  assumed distribution                   >=75%    >=90%    >=99%
  independent uniform positions           76.225%  34.370%   0.505%
  grouped near one common center            0.000%   0.000%  0.000%
  near evenly spaced with .1 jitter       100.000% 100.000% 87.405%

The Monte Carlo's observed means of the ratio S were
  uniform: 0.18050,
  clustered: 1.33973 (average curvature often *increases*),
  evenly spaced: 0.00475.

The regular arrangement is highly special; it has
systematic cancellations of low-k Fourier harmonics.
It is not necessarily favored by stochastic primordial
physics. Likewise the one-cluster prior concentrates
material more strongly and gives *less*, not more,
mean-curvature suppression in this toy.

This directly contradicts claims that simply packing
material into small tight gravitational wells inevitably
makes the universe globally flatter.

The 90% sample proportion for independent uniform
positions is about .3437, with ordinary binomial
sampling uncertainty ≈0.00336 (1σ); its 95% Wilson
interval is roughly .337–.350. This quantifies numerical
sampling noise only, NOT uncertainty in physical priors.

For a fixed sigma with independent uniform positions,
the expected weak-curvature ratio has a closed form:
  E[S] = (1/8)sum_k exp(-k²sigma²)/k²
               /[sum_k exp(-k²sigma0²)/k²].

There is no requirement for the number of gravitational
wells to equal a number of Planck-mass units.

## What can and cannot be concluded

The computation verifies that, *within an explicitly
chosen restricted ensemble*, many nonuniform
distributions of the same conserved matter lead to
relatively small averaged intrinsic structural curvature.
It also shows the fraction ranges from approximately
0 to 100% when the assumed spatial correlation prior
is changed. Such wide sensitivity makes a unique
"physical probability of near-flatness" impossible
without the primordial ensemble measure and the
ADM/structural-to-Jordan observable mapping.

Formally the desired answer is

 P_nearflat =
  (∫ dmu[I] 1_physically_allowed(I)
     1_{|Omega_k,J,0(I)|<epsilon})/
   (∫ dmu[I] 1_physically_allowed(I)).

To evaluate it physically the full SDMC genesis
must specify:
 - a gauge-reduced and normalized measure dmu[I]
   at or after the allowed earliest classical slice;
 - constraints at S_P=1: rho_v,P=rho_P, chi_P=0,
   and allowed material creation/redistribution history;
 - initial scalar/metric and matter correlation spectra
   selected by the thermal transition;
 - positive density, proper fluid evolution, causality
   and stable F/G2/G3/structural kinetic dynamics;
 - Hamiltonian/momentum constraints, evolving shear and
   spatial curvature from the SAME microphysical solution;
 - a nonflat Jordan observable curvature prediction
   through late300 and CMB/BAO distance likelihoods.

Conclusion: the screenshot's "1 in 10^60" is not a physical
probability. SDMC's material-redistribution interpretation
offers many mathematical near-flat examples, but physical
near-flatness probability is not yet identified by the
available action and input distribution.

Code: sdmc/flatness_probability_ensembles.py
Tests: tests/test_flatness_probability_ensembles.py
