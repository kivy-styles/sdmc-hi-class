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
