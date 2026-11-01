# SDMC Manuscript B — Part III: Primordial-Well Occupation and Scalar–Geometric Response

## 500. The surviving curvature-to-scalar-to-strain pathway

The preceding source audit established that material, minimally
coupled to the Jordan metric, does not directly source the independent
relative strain L in the *displayed* action. The curvature-selector
mixing term is proportional to W and vanishes when W=0. On the
unexcited locked-flat branch the positive quadratic L operator
admits no sourced nonzero static lump.

This does not settle the response about an *occupied* structural
state. The accepted F(sigma)R_J/2 part of the rolling scalar action
receives a Jordan scalar-curvature source even when W=0. The
supporting-space kinetic/restoring coefficient B f_L² can depend
on sigma. An already active relative-strain configuration then
responds to changes of sigma through the same action.

No extra material-to-L force and no CDM fluid is introduced
by this calculation. It probes an existing interaction that
is absent to first order only about L=0.

## 501. Reduced Jordan-time mode action and exact background equation

With a real spatial cosine ell=A(t_c)cos(kx), W=0, and
the existing relative-trace operator, write

 S_ell=int dt_c a_J³ (K_ell(sigma)/4)
          [dot A²-omega_L²(a_J)A²],

 K_ell=B f_L²/3>0,
 omega_L²=m_L²+c_L²k²/a_J².

This section treats the mass and gradient coefficient as
independent of sigma to isolate the K_sigma interaction.
If those coefficients also depend on sigma the final
source includes a term proportional to their derivatives.

Euler–Lagrange variation produces

 ddot A+(3H_J+dot K_ell/K_ell)dot A+omega_L² A=0.

This equation uses the Jordan time t_c; a_s and t_s
are not silently substituted for observable Jordan clocks.

## 502. Linearization around an occupied structural mode

Let
 sigma=sigma_bar+epsilon delta sigma,
 A=A0+epsilon delta A.

At fixed background Jordan spatial geometry,

 delta ddot A
 +(3H_J+dot K0/K0)delta dot A
 +omega_L² delta A
 =
 -delta(dot K/K) dot A0
 -delta(omega_L²) A0.

For the conditional coefficient history
 K_ell(sigma)=K0 exp(beta sigma)
and sigma-independent omega_L²,

 delta(dot K/K)=beta delta dot sigma,

so

 delta ddot A
 +(3H_J+beta dot sigma_bar)delta dot A
 +omega_L² delta A
 =-beta dot A0 delta dot sigma.

This is a source linear in the *preexisting mode velocity*.
If A0=dot A0=0, delta A remains zero from zero
perturbation initial data. If beta=0 or delta sigma is
constant, this particular interaction does not drive
delta A. If A0 is occupied and rolling, scalar changes
can parametrically modulate it.

The exact earlier conclusion remains valid:
no stationary new relative strain arises on the
unexcited W=0 branch through the displayed terms.

## 503. Scalar curvature-loading is present but not yet physically solved

The covariant scalar–tensor sector has
 S_J ⊃ int sqrt(-g_J) F(sigma)R_J/2,

so the sigma equation includes F_sigma R_J/2.
For illustration, the numerical subcase uses a
canonical scalar disturbance s with a prescribed
small Jordan spacetime Ricci impulse delta R_J:

 ddot s+3H_J dot s
      +[m_sigma²+c_sigma²k_sigma²/a_J²]s
       =(F_sigma/(2K_sigma))delta R_J.

The fully accepted late300 F/G2/G3 scalar dynamics
are not identical to this reduced scalar equation.
Furthermore delta R_J must be computed from the
physical coupled Jordan lapse/shift, matter and
metric equations. The numerical impulse is an
EXTERNAL CONTROL INPUT, not a prediction from
the primordial material inventory.

In e-fold x=ln a_J with
 H_J=H0 exp(-h x), sigma_bar=v x,

 s_xx+(3-h)s_x
     +(omega_sigma²/H_J²)s
      =F_sigma delta R_J/(2K_sigma H_J²),

 A0_xx+(3-h+beta v)A0_x
     +(omega_L²/H_J²) A0=0,

 delta A_xx+(3-h+beta v)delta A_x
     +(omega_L²/H_J²)delta A
       =-beta s_x A0_x.

This permits a clean linear-response experiment with
a distinct causal role for each interaction.

## 504. Direct amplitude/energy tangent checks

The numerical work solves all three equations
simultaneously and independently evolves the exact
nonlinear exponential-K amplitude equation at
sigma=sigma_bar±epsilon s.

For the quadratic mode-energy observable

 E_L=K_ell/4[H_J² A_x²+omega_L² A²],

its first variation is

 delta E_L =
   beta s E_L0
   +K0/2[H_J² A0_x delta A_x+
                    omega_L² A0 delta A].

The first term is the direct coefficient variation;
the second comes from the actual excited mode's
linear dynamics.

These two contributions can interfere; a nonzero
instantaneous delta E is not automatically a
new net cosmological material-density source.

A centered ±epsilon finite-difference derivative
of the full exponential-K action's mode equation
is compared to the analytically linearized delta A
and delta E.

Essential negative controls:
- If A0=A0_x=0, the induced geometric trace and
  its energy response vanish even though the
  canonical scalar curvature disturbance exists.
- If beta=0, the scalar does not drive this L mode.
- If F_sigma=0 or deltaR_J=0, this canonical
  curvature-loading route produces no delta sigma
  and no resulting delta A under zero initial
  scalar disturbance.
- If an occupied A0 is present, the response scales
  linearly with A0 while E0 scales quadratically.

These checks are not full physical metric-potential
tests and are not late300 no-CDM likelihood results.

## 505. An exact continuity identity for the variable-K mode

An additional consistency result follows directly
from the W=0 quadratic trace action.

Define G=c_L²k²/a_J² and omega_L²=m_L²+G.
The reduced mode energy and isotropic pressure are

 E=K/4[dot A²+(m_L²+G)A²],
 P=K/4[dot A²-m_L² A²-G A²/3].

An independent auxiliary pressure-like quadratic
factor in the scalar energy exchange is

 Pstar=K/4[dot A²-(m_L²+G)A²].

Using the *exact on-shell* A equation,

 dot E+3H_J(E+P)
       =-(dot K/K)Pstar.

Equivalently, in x=ln a_J,

 E_x+3(E+P)
       =-(K_x/K)Pstar.

This is an action-derived reduced energy-transfer
identity, not a phenomenological adjustment.

The gradient contribution to the actual isotropic
spatial pressure is -G A²/3, not -G A². This exact
distinction makes the scalar energy transfer
identity cancel the ordinary gradient redshifting
when K is constant.

For a massive high-frequency adiabatic oscillator,
phase averaging gives
 <dot A²>≈m_L²<A²>
and negligible G fraction, so <Pstar>≈0.
Consequently an adiabatic scalar-dependent K does
not necessarily create a large averaged geometric
energy transfer even if it modifies the
instantaneous oscillating field.

This is a new bottleneck for claims that K_sigma
alone creates the missing cosmological
matter-like abundance from an empty state.

## 506. Relation to primordial gravity, rho_v, rho_X and rho_g

At the exact Planck boundary the manuscript states

 S_P=1, rho_v,P=rho_P, rho_X,P=0, chi_P=0.

Classical F/G2/G3 plus spatial-strain equations
apply on the one-sided post-Planck branch,
S_A approximately 1; no pre-Planck ordinary
FLRW extrapolation is made here.

The primordial original material-inventory identity

 M_u=M_remaining+M_inferred

cannot determine A0 or its velocity, the
geometric occupation number, the support-space
spectral correlations, or the physical scalar
coupling coefficient beta.

The accepted bare/active split
rho_X=chi rho_v
must not be confused with the separate
supporting-space rho_g, and rho_v must not be
counted as a third classical fluid.

The present calculation establishes a conditional
existing-field response:
  primordial occupied L
     -> scalar-dependent amplitude/energy modulation.

It does not demonstrate
  primordial M_u -> nonzero occupied L
or
  occupied L -> observed Phi_J/Psi_J
at the required amplitudes.

## 507. Physical closure tests and falsifiability

The physical no-CDM hypothesis requires:
 (i) a valid post-Planck canonical state that
     fixes the preexisting structural mode
     occupation and correlations;
 (ii) a self-consistent rolling F/G2/G3
     scalar trajectory and curvature source;
 (iii) one ADM lapse/shift/metric variation
     including L background energy and exchange;
 (iv) positive kinetic and stable gradients,
     isotropic averaged stress and causal
     perturbative propagation;
 (v) the correct Jordan Poisson/Weyl potentials,
     equality and CMB transfer from the same
     geometric state without a reference CDM
     energy density added in the equations.

The present source experiment satisfies none of
these observational closure requirements. It
tests the formal existence and the dynamical
selection rule of one existing interaction.

A central emerging result is that a source-free
empty relative geometry stays empty at W=0,
while a preoccupied geometry can respond to
the scalar even without a direct L-matter
interaction. This is a conditional
occupation-dependent response, not a
derivation of initial occupation.

Core code:
  sdmc/pg_occupied_scalar_curvature_transfer.py
  tests/test_pg_occupied_scalar_curvature_transfer.py
  .github/workflows/sdmc_pg_occupied_well_scalar_transfer.yml

Source references: SDMC Manuscript B, Part III
Sections 101–106, 111–115; supporting-space
ADM and rolling G3 follow-up calculations.


## 508. Can an initially empty classical mode become quantum-occupied during activation?

Section 115 of Manuscript B identifies the missing physical
canonical initial data and points out that B=0 makes the isolated
supporting-space kinetic rank degenerate. The existing technical
appendix considers a canonical Bogoliubov mode occupation.

A strictly CLASSICAL zero initial condition A=dot A=0 remains zero
under the homogeneous quadratic W=0 equation no matter how
K(t) varies. This is a uniqueness result for the ordinary
classical initial-value problem.

If the structural mode admits an independent quantum
description on a valid positive-rank post-Planck branch,
a different mechanism is possible: a time-dependent
canonical normalization can mix positive and negative
frequency vacuum modes, leading to a nonzero expected
particle/mode occupation in the final basis.

This is not a classical spontaneous source and does not
remove the requirement for an initial quantum state.

## 509. Exact canonically normalized oscillator

For one real mode,
 S_A = integral dt_c (a_J³ K_ell/4)
                 [dot A²-omega_L² A²].

Define the positive normalization
  b(t_c)=a_J³ K_ell/2,
  v=sqrt(b) A,
  Gamma=dot b/b=3H_J+dot K_ell/K_ell.

After an integration by parts the quadratic action becomes

  S_v = integral dt_c 1/2 [
          dot v² - Omega_v² v²],

  Omega_v² =
    omega_L² - (1/2) dot Gamma - (1/4)Gamma².

This is an exact result of the reduced action.
It is NOT evidence for a particular quantum state.

When B→0 the variable v becomes degenerate and
the quadratic canonical description loses validity.
Therefore a positive initial B_floor is required
for the concrete numerical diagnostic.
No pre-Planck FLRW equation is assumed.

## 510. Controlled local kinetic-ramp occupation experiment

To separate kinetic-ramp production from uncomputed
cosmological expansion, the experimental benchmark
takes locally frozen Jordan a_J=1 and H_J=0,
W=0, positive m_L and c_L², with

  omega_k²=m_L²+c_L² k²,

  K_ell(t) =
   K_i exp[ln(K_f/K_i)u(t)],

  u(t)=[1+exp(-t/tau)]^-1.

The normalized canonical initial positive-frequency
mode is

  v_i=1/sqrt(2 omega_k),
  dot v_i=-i sqrt(omega_k/2).

The canonical equation is integrated through a
finite smooth transition,

  ddot v + Omega_v²(t)v=0.

When dot K/K vanishes at the late endpoint,
the mode can be decomposed into final
positive/negative frequencies:

  alpha_k=sqrt(omega_k/2)v_f
             +i dot v_f/sqrt(2omega_k),

  beta_k=sqrt(omega_k/2)v_f
             -i dot v_f/sqrt(2omega_k).

The mode occupation is

  n_k=|beta_k|².

The code verifies the canonical identity

  |alpha_k|²-|beta_k|²=1

and the independent late energy identity

  E_final - omega_k/2 = omega_k n_k.

For an unchanging K_i=K_f, n_k≈0, as required.

For the same finite ratio K_f/K_i=5,
a rapid but positive-frequency transition
produces a nonzero mode occupation while a
slow transition remains nearly adiabatic.

These numbers measure the response of the
chosen quantum oscillator to an input
coefficient profile. They are not predicted
relic abundances or proof of SDMC dark matter.

## 511. Consequence for the primordial gravitational well

The classical empty-branch source audit and
quantum-occupation diagnostic are compatible:

  Classical zero initial strain + homogeneous
    positive-rank quadratic action -> zero.

  Assumed initial positive-frequency vacuum
    + nonadiabatic post-Planck K(t) -> possible
    final nonzero canonical occupation.

Thus the absence of a direct late material-to-L
force does not logically exclude production of
geometric quanta while the kinetic sector
activates. But this conditional route is
only meaningful if the canonical field remains
physical after all ADM gauge and secondary
constraints and if a valid initial vacuum
can be specified near S_A≈1.

The bare structural density rho_v,P=rho_P at S=1
does not fix that vacuum, the transition rate,
the physical UV cutoff or the spectrum n_k.

A phenomenological CDM fraction may not be
substituted for any missing normalization.

## 512. Cold-like dilution and abundance remain separate

A nonrelativistic structural excitation with
physical mode energy approximately m_L could
contribute a mean energy density of the form

 rho_g(a_J) ~ a_J^-3
      integral d³k/(2pi)³ omega_k n_k,

only after the physical gauge-invariant modes,
measure and renormalized stress have been
established.

The numerical single-mode occupation test does
not evaluate this integral and does not determine
the observed additional gravitational abundance.

Even if an initial population exists, the
effective perturbative sound speed and the
Jordan Poisson/Weyl potentials must be calculated
from the same full tensor/scalar action.
The action-level oscillation condition is not
a gravitational lensing result.

## 513. The revised physical bottleneck

The source investigation has isolated two
routes with distinct unresolved requirements:

Route A: occupied classical supporting geometry.
Curvature and sigma can modulate preexisting
geometric modes, but the primordial occupation
and its normalization remain uncomputed.

Route B: quantum creation during a finite
post-Planck kinetic handoff.
A positive-rank canonical EFT permits a
nonzero Bogoliubov mode occupation from a
nonadiabatic coefficient evolution, but its
physical initial state, coefficient history,
gauge constraints and UV cutoff remain to
be derived from SDMC.

Neither route establishes a permanent static
relative-strain well sourced by M_u after W=0.
Neither route establishes that the scalar or
geometric fields reproduce CDM gravitational
observables.

The decisive next calculation requires a
complete 4D supporting-space parent action
and valid quantum/classical matching on the
post-Planck S_A≈1 branch, followed by a
self-consistent matter+radiation+scalar
Einstein–Boltzmann evolution.

Additional code:
  sdmc/pg_canonical_handoff_occupation.py
  tests/test_pg_canonical_handoff_occupation.py


## 514. Structural energy bookkeeping and the bare/active split

A time-dependent K_ell is a background source of canonical
oscillator work, not a spontaneous creation of physical
energy from nothing. The full scalar-structural action
must include the backreaction of the produced mode
population on the field that controls B f_L², including
the scalar equation and the Jordan Hamiltonian constraint.

If the occupied geometry ultimately forms a separately
conserved cold-like source in a physical FLRW solution,

  rho_g=Gamma_g rho_v,
  rho_v=rho_P S^-2,
  p_v=d ln S/d ln a_J,

then exact dustlike dilution imposes

  rho_g'=-3rho_g,
  Gamma_g'=(2p_v-3)Gamma_g,

where prime means d/dln a_J.
During ideal radiation structural tracking p_v=2,
Gamma_g proportional to a_J; during matter tracking
p_v=3/2, Gamma_g is constant; during mature coasting
p_v=1, Gamma_g proportional to a_J^-1.

This is an IDENTITY, not a primordial abundance
prediction. A kinetic quench can establish a
nonzero single-mode n_k under assumed quantum initial
conditions, but cannot by itself supply the
normalization of Gamma_g, the full renormalized
three-dimensional spectrum or its physical pressure
and metric response.

To achieve ordinary matter/radiation equality near
the comparison z~3466 instead of a baryon-only
z~527 under the same reference radiation history,
the geometric population must be present and
effectively cold sufficiently early. The physical
acoustic peaks additionally require the correct
potential evolution and photon-baryon forcing.

No full Planck or lensing likelihood is claimed
for the illustrative source and quench calculations.


## 515. Fourier/gauge restriction on the reduced modulation experiment

The source pulse is a TIME-ONLY envelope in the local
mode experiment. The mode wavenumber k_sigma appears
in the illustrative scalar oscillator frequency, but
the numerical system does not solve a spatial
convolution of an inhomogeneous delta sigma_k(x)
with an independent ell_q(x). In a full cosmological
perturbation theory, their product generally couples
different spatial Fourier modes and the gradients of
K(sigma) also contribute to the spatial field equation.

Consequently the correct status is a local
parametric-response demonstration, not a
gauge-invariant k-by-k Poisson or lensing kernel.
The scalar response also omits feedback of
the excited strain stress into sigma, lapse,
shift and the Jordan metric. This feedback is
essential for a physical energy-transfer claim.

The numerical exact Ward identity checks the
specified reduced mode action, not the completed
rolling Horndeski fluid or spatially varying
tensor-metric stress.
