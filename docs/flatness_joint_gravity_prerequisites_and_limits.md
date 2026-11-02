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


## 5. A further restriction: the same toy initial data must satisfy GR constraints

The earlier 1.60% / 2.16% geometric intersection is still only a
necessary proxy. A genuine allowed classical initial state must
satisfy a Hamiltonian and momentum constraint *and* energy positivity.

As a controlled comparison (not an SDMC UV solution), equip the
SAME 60,000 mass-position-and-width seed configurations with
plane-symmetric conformal spatial metric

  h_ij=e^{2q(x)}delta_ij,

where q=-Phi is calculated from the seed periodic Poisson density
and is treated as an initial geometric template. Use unit a_J=a_s
in an arbitrary length convention.

Impose irrotational GR dust constant-mean-curvature initial data

  K_ij=-H h_ij, H constant in each realization,
  sigma_ij=0, J_i=0.

Then the momentum constraint D_j(K^ij-h^ij K)=8piG J^i
is exactly satisfied. The GR Hamiltonian is

  R3+6H²=16piG rho_phys.

The density rho_phys is reconstructed from the nonlinear geometry.
It need not equal the Gaussian seed density; the seed fixed q only.

The average physical volume is V_phys=Vcoord<e^{3q}>, and for
prescribed identical conserved material mass M_u=1, solve

  H² = [16piG M_u/V_phys - <R3>_phys]/6.

Then

  rho_phys(x)=[R3(x)+6H²]/(16piG),

and check H²>0, rho_phys(x)>=0, total integral rho_phys
sqrt(det h) d³x=M_u with numerical residual near machine precision.

This builds genuine GR dust initial slices satisfying both initial
constraints, given the additional freely selected conformal
metric and CMC assumptions. It is not the action of SDMC,
not a Planck-boundary initial state, not an early radiation
fluid and not a dynamical transition/evolution.

The same Fourier ensembles, seed 4102, widths Uniform(.1,.6)
and original geometric prerequisites were retained. Across
20,000 samples of each:

                           P(F and W)   P(CMC positive)
  independent uniform        1.60%        60.545%
  one common cluster         0.00%         0.000%
  almost evenly spaced       2.16%        96.860%

Adding the GR CMC constraints to F and W:

  independent uniform: 2 / 20,000 =0.0100%;
  one common cluster:  0 / 20,000 =0%;
  almost evenly spaced:80 / 20,000 =0.4000%.

Among the GR-CMC-admissible subsets these become
  independent uniform: 2/12,109 ≈0.0165%;
  clustered: no allowed samples in this restricted construction,
             so conditional probability is undefined;
  almost evenly spaced:80/19,372≈0.4130%.

These are again strictly **toy prior fractions**, NOT physical
Planck SDMC near-flatness+CDM-like success rates. Requiring
a CMC dust slice itself is an unmotivated extra restriction
for the radiation-era early universe. It is informative
precisely because a constraint filter drastically changes
the joint proxy fraction and because the Gaussian seed
density is not the same as the reconstructed exact GR
density.

The special Buchert geometric source on these CMC slices has
  Q_D=0 (uniform expansion, zero shear);
  rho_g,av=-R_D/(16piG)<=0;
  P_g,av=R_D/(48piG).

Indeed, beyond CMC, on ANY exactly locked periodic conformal
scalar spatial slice with R_D>=0, requiring a positive
Buchert geometric density demands Q_D<-R_D. Then

  w_g,av=P_g,av/rho_g,av
       =1+4R_D/[3(-R_D-Q_D)] >= 1,

strictly >1 for R_D>0. Thus under these locked-metric,
periodic boundary, dust Buchert effective-source assumptions,
positive effective geometric energy is incompatible with
pressureless behavior. One cannot rescue positive
pressureless Buchert dust by merely adding enough shear
while retaining the same nonnegative average intrinsic R_D.

This is a strong restriction on the PARTICULAR Buchert
effective-dust route. It does NOT exclude a distinct
positive action-derived SDMC supporting-space rho_g,
an independent traceless metric Q_i^j or nonconformal map,
a negative reference spatial curvature, nonperiodic
boundary effects, or an evolving non-dust/radiation
structural–Jordan geometry.

In particular "CDM-like success" is still NOT demonstrated
for any of these GR-CMC profiles: physical scalar/structural
stress, growth, sound speed and lensing remain absent.

Code: sdmc/flatness_joint_cmc_constraints.py
Tests: tests/test_flatness_joint_cmc_constraints.py

## 6. Final status of the requested fully simultaneous calculation

We now have a layered same-sample test:
  1. fixed total material inventory M_u;
  2. specified weak-reference averaged structural curvature proxy;
  3. nonzero small-scale material power;
  4. exact GR comparator initial Hamiltonian/momentum constraints;
  5. exact restriction on the special positive dust Buchert branch.

The full physically requested joint event
 |Omega_k,J,0|<eps_obs AND physical SDMC cold-like gravity
remains UNCOMPUTED, since neither indicator follows from
the reduced configuration ensemble and no normalized
physical primordial measure or complete coupled 4D
SDMC nonflat primordial trajectory has been supplied.

This establishes a quantitative **constraint-based limitation**,
rather than any measured primordial probability.
