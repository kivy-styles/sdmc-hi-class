# SDMC Manuscript B continuation — Joint flatness and Part III geometric-gravity prerequisites

## Central question

What fraction of *physically allowed* Planck-to-radiation initial configurations
produce BOTH observational near-flatness and the CDM-like geometric response
investigated in Manuscript B, Part III?

The actual physical fraction is **not yet calculable** from the available
model because (i) the normalized initial-state measure for admissible
prethermal SDMC solutions is absent, (ii) no self-consistent non-flat
primordial 4D action solution produces the post-handoff material/metric
spectrum, and (iii) a full observable Jordan pressure/growth/lensing
response of the supporting-space stress is not derived.

This is not a reopening of the CDM question: it applies the conditions
already defined in Part III as a filter on the flatness ensemble.

Two different metrics and criteria must not be conflated:

  near-flatness: |Omega_k,J,0| < eps_obs, for observational Jordan
                 FLRW curvature reconstructed from full distances and H;
  cold-like: rho_g>0 with effectively P_g≈0, rho_g~a_J^-3,
             small effective sound speed, correct growth and lensing,
             and a physical abundance consistent with the accepted
             baseline, without a second CDM inventory.

Neither is equivalent to a density profile being narrowly peaked.

## 1. Reanalyze the SAME ensemble of configurations

Use the last probability study's exact priors, not independent
probabilities multiplied across different experiments:

 - 8 positive equal-mass concentrations, M_u conserved;
 - positions independently uniform, near one common center, or
   close to evenly spaced 2pi j/8 with jitter;
 - every configuration draws sigma from Uniform(.1,.6);
 - initial reference is one concentration of width .6;
 - first 100 Fourier modes, 20,000 draws per class, seed 4102.
 - plane-symmetric periodic spatial field; not a physical
   isotropic 3D Planck-era structure ensemble.

Each material distribution is represented by

  delta_k = exp(-k²sigma²/2) * 1/8 sum_i exp(-ik x_i),
  k>=1.

This is the normalized density contrast for the test toy;
mass k=0 is fixed and never produces a new material inventory.

### Structural near-flatness proxy

On the flat-reference locked conformal scalar branch with weak
subhorizon Poisson sourcing,

  <R3_s>_V^(2) ∝ sum_{k>=1}|delta_k|²/k².

Define S as this weighted sum relative to the initial
single-reference Gaussian width .6.

  criterion F := S<0.1

is ≥90% *averaged structural Ricci-curvature* reduction.

**F is NOT |Omega_k,J,0|<.004** and is not an observational
constraint by itself. For an initial nonflat reference,
additional background spatial curvature terms must be included.

### Local-well geometric prerequisite

For the same configuration, define

  D_rms=[2 sum_{k=1}^{100}|delta_k|²]^{1/2},

and the short-scale density power fraction

  f_short =
      [sum_{k=8}^{100}|delta_k|²]/
       [sum_{k=1}^{100}|delta_k|²].

The proxy W requires

  D_rms >= 1, f_short >= .05.

This requires substantial density contrast and at least
5% power at wavelengths represented by k>=8.

The thresholds are deliberately disclosed and arbitrary,
NOT derived physical cold-matter/lensing requirements.
For our Poisson trial, linear local structural Ricci
curvature follows density contrast, but this does NOT
show the correct action-derived Jordan gravitational force.

The *joint geometric proxy* is J=F and W.

## 2. Numerical overlapping fractions

Exactly 20,000 samples per prior, 60,000 in total.
The data below count the same samples in F and W.

                                P(F)      P(F and W)
  independent uniform          34.37%       1.60%
  one clustered center          0.00%       0.00%
  near evenly spaced          100.00%       2.16%

Conditional on F and an independently random position prior,
only approx 1.60/34.37 = 4.655% also meet the local-well proxy.

The finite-sampling 95% Wilson interval for 1.60% is about
[1.435%,1.783%]; for 2.16% it is [1.968%,2.371%].
These intervals quantify only random Monte Carlo uncertainty,
not physical prior uncertainties.

The joint prerequisite is sensitive to local-density threshold:

  D_rms threshold         uniform      clustered     near-even
       .5                    3.885%        0%         17.44%
      1.0                    1.600%        0%          2.16%
      2.0                    0%            0%          0%
with f_short>=.05 and S<.1 throughout.

In the exactly equally spaced eight-concentration limit,
  sum_{j=0}^7 exp(-2pi i k j/8)=0, k mod 8 != 0,
but modes k=8,16,... survive, especially for small sigma.

Thus large-scale homogenization and local inhomogeneity
can coexist mathematically. Broad wells, however, can
have negligible surviving short-scale contrast even when
low-k curvature vanishes, and strong short-scale density
contrast does not itself prove pressureless geometric gravity.

## 3. Check against the ACTUAL Part III gravitational dust branch

Manuscript B Part III §107 gives the GR Buchert irrotational
dust comparison

  rho_g,av = -(R_D+Q_D)/(16piG),
  P_g,av   = (R_D-3Q_D)/(48piG).

For that *specific* averaged-GR realization,
P_g,av=0 requires R_D=3Q_D. Positive rho_g,av then needs

  R_D=3Q_D<0,
  rho_g,av=-R_D/(12piG)>0.

For the restricted periodic conformal spatial geometry used
in the earlier structural flatness test,

  gamma_ij=a_s² exp(2q)delta_ij,

the proper-volume averaged R is

  <R3_s>_V=2/a_s² *
     [int exp(q)|grad q|²]/
     [int exp(3q)] >=0.

If we additionally impose **exact spatial conformal locking**
gamma=M² h_J on the same averaged spatial foliation,
the physical Jordan spatial Ricci scalar has the same sign.
Hence the positive averaged Ricci curvature of those
nontrivial scalar-only periodic profiles cannot realize
**positive-energy Buchert geometric dust** under R_D=3Q_D.
The fraction satisfying both the positivity condition
and that special dust branch is exactly zero for this
restricted conformally locked scalar-periodic family.

THIS IS NOT A GENERAL SDMC NO-GO. Manuscript B explicitly
distinguishes the supporting spatial metric and its
Jordan effective stress from a GR Buchert average of the
Jordan physical geometry. Nonconformal Q^i_j, nonlinear
shear/correlations, signed large-scale curvature, different
slicings, and the direct action-derived geometric source
are not captured by our flat-reference scalar-only toy.
The manuscript's §107 itself cautions that periodic Newtonian
backreaction can cancel to Q_D=0 at leading order.

The Part III direct-action route instead needs a physical
rho_g>0, pressure, energy-transfer law, linear clustering
and Jordan metric potential modifications determined by
the full F/G2/G3+UV ADM action. A surviving localized
material-density mode does not certify these criteria.

## 4. Interpretation and what remains to be calculated

P(J|prior) as computed above is not the desired physical
joint probability. The latter is

 P_both =
  [integral over I dmu(I) 1_admissible(I)
     1_{|Omega_k,J,0(I)|<eps_obs}
     1_{positive ColdLikeGravity(I)}]/
  [integral over I dmu(I) 1_admissible(I)].

The required cold-like indicator must include
  (a) physical rho_g>0 with a justified normalization;
  (b) dustlike pressure and a^-3 density dilution;
  (c) low physical sound speed/adequate growth through radiation
      and matter epochs;
  (d) correct Jordan Weyl lensing, CMB and BAO;
  (e) exact conserved M_u with no independently inserted CDM.

None can be inferred solely from the density Fourier spectrum.

A future actual joint test must use one *shared* primordial
Planck-to-radiation state measure, solve the fully constrained
nonflat handoff for each configuration, calculate both Jordan
curvature and structural response on that same solution,
then count those meeting the observational requirements.

Consequently the mathematical answer is:
  - joint geometric prerequisite can be quantified;
  - special positive-energy Buchert dust branch is excluded
    in the stated locked conformal scalar periodic toy;
  - physical simultaneous flatness and CDM-like success
    probability is presently **undetermined**.

## Reproducibility

Code: sdmc/flatness_joint_cold_necessary.py
Tests: tests/test_flatness_joint_cold_necessary.py
Comparison: sdmc/flatness_probability_ensembles.py
Physics reference: Manuscript B Part III Sections 102, 106–111, 116.
