# SDMC Manuscript B, Part III — Flatness Selection and the Survival of Geometric Gravity

## 476. The two physically separate tasks of the spatial sector

Manuscript B, Part I (Sections 43–49 and 53–56) retains a flat
Jordan cosmological baseline with K_J=0. Its candidate microscopic
spatial completion adds a relative structural strain L and a
transient intrinsic-curvature selector.

A selected flat background and local curved gravitational wells
are *not contradictory*:
- K_J=0 is the homogeneous/isotropic Jordan FLRW curvature index.
- The inhomogeneous structural metric gamma_ij may have nonzero
  Riemann[gamma] even when the Jordan background is flat.
- A nonzero inhomogeneous mode stress has to backreact on the
  Jordan lapse, shift and metric; otherwise the background
  is not a solution of the fully coupled action.

The question is whether the same existing supporting-space
action can (i) favor the flat Jordan branch during its formation
and (ii) retain an excited geometric gravitational sector
after the curvature selector switches off.

This study is an analytically controlled and numerically tested
**conditional handoff**. It is not a complete Planck-to-present
field-equation solution.

## 477. Homogeneous curvature selector

The adopted candidate SDMC intrinsic curvature selector has

  V_K[gamma]=eta_K/2 * R_ijkl[gamma]R^ijkl[gamma].

On a locked homogeneous constant-curvature FLRW structural
slice with gamma_ij=(a_s²) q_ij(K_J),

  R_ijkl[gamma]R^ijkl[gamma]=12 K_J²/a_s⁴,

so the curvature cost is

  V_K^FLRW=6 W eta_K K_J²/a_s⁴.

For W eta_K>0, the minimal energetic value at a fixed
structural scale is K_J=0. The source material establishes
this as the candidate curvature-selection potential.

But a constant-curvature FLRW index is not automatically a
time-dependent propagating scalar degree of freedom. With
a mature fixed spatial metric topology, merely adding an
energy cost proportional to K_J² does not make K_J
continuously evolve through -1, 0, +1.

Therefore the selection is physically a PRE-FLRW /
Jordan-geometry-formation process unless a complete
nonflat spatial covariant parent action derives a lawful
inter-branch transition or selection probability.

The post-handoff kinematic FLRW identity remains

  Omega_K=-K_J/(a_J² H_J²),
  d Omega_K/d ln a_J=2q_J Omega_K,

where q_J=-1-dot H_J/H_J². Thus K_J=0 is an exact
invariant branch; it is not necessarily a stable attractor
from nonzero K_J. During decelerating matter-like expansion
q_J≈1/2, the *relative* magnitude of a nonzero
Omega_K tends to grow approximately as a_J.

The corrected manuscript's flatness-status sections
explicitly distinguish accepted K_J=0 calculations
from hypothetical K_J !=0 curved extensions.

## 478. Exact locking versus occupied geometric excitations

The accepted candidate spatial action contains

  S_support=int dt_c d³x N sqrt(h) [
      B f_L²/2 Tr(D_nL D_nL)
     -B f_L² c_L²/2 Tr(D_iL D^iL)
     -B f_L² m_L²/2 Tr(L²)
     -W eta_K/2 R_ijkl[gamma]R^ijkl[gamma] ].

The homogeneous energy minimum is L=0,
R_ijkl[gamma]=0 on the selected flat branch.

If *everywhere* L=0 and its derivatives vanish, then
rho_support=0 and no late gravitational abundance follows
from this relative strain sector.

However, *near-locking* L !=0 with small amplitude and
nonzero occupation can carry positive quadratic energy
even if the selector W ->0 after the early handoff:

  rho_support≈(B f_L²/2)[ Tr(dot L²)
                 +m_L² Tr(L²)
                 +c_L²/a_J² Tr((grad L)²)].

The surviving energy is not the early selector energy,
which tends to zero when W->0 at bounded curvature.

Consequently a 'flat branch' and an 'occupied
structural geometry' are logically compatible.
They are NOT the exact L=0 vacuum branch treated by
strict IR integration-out calculations.

## 479. Trace mode equation with activating lock B and transient selector W

For one real spatial cosine trace
  ell(x,t_c)=A(t_c) cos(kx)

and f_L² constant over the illustrative handoff, define

  K_trace(x)=B(x) f_L²/3,
  a_s=S/S0, x=ln a_J,
  H_J(x)=H0 exp(-h x), a_s(x)=exp(p_v x),

  C_k=4 eta_K W(x)k⁴/a_s⁴,
  Lambda_k=2 C_k/(9 K_trace),

  Omega²=m_L²+c_L² k²/a_J²+Lambda_k.

The post-Planck Jordan-mode Euler–Lagrange equation is

  ddot A+(3H_J+dot K_trace/K_trace)dot A
      +Omega² A=0.

In e-fold time this becomes

  A_xx+[3-h +dln K_trace/dx] A_x
     +Omega²/H_J² A=0.

The module deliberately starts at B>0, because at B=0 the
quadratic kinetic term degenerates and the classical mode
equation cannot be extrapolated across that point without
a UV completion. It does not extend classical FLRW to
times earlier than the exact Planck boundary S=1.

An illustrative smooth handoff is

  u=1/(1+exp[-(x-x_h)/Delta]),
  B=B_floor+(1-B_floor)u,
  W=4 W_peak u(1-u).

The manuscript gives a similar schematic pulse.
Neither the shape nor its time placement is an accepted
SDMC prediction.

## 480. Adiabatic survival despite strong locking

For a massive structural oscillator with

  m_L >> H_J,  m_L² >> c_L²k²/a_J²+Lambda_k,

and slowly varying B, the mode action invariant is

  I= a_J³ <rho_trace>/Omega ≈constant.

The reduced structural-mode energy per Jordan physical
volume scales

  <rho_trace>≈a_J^-3 Omega.

For constant m_L this approaches

  <rho_trace> ~ a_J^-3.

The oscillatory amplitude obeys

  A_envelope ∝
    [a_J³ K_trace Omega]^-1/2.

Thus increasing B damps the *relative geometric strain
amplitude* but does not necessarily erase the already
present occupation energy. This is a physical difference
between a heavy field in its vacuum (integrate out) and
a heavy field carrying a coherently oscillating population.

The module numerically integrates two orthogonal phases
of the same existing structural trace mode and compares
their average energy, envelope and adiabatic invariant.
The W=0 mass-restoring baseline, an illustrative transient
W pulse in a matter-like scaling, and a second pulse on the
radiation tracker p_v=2, H_J proportional to a_J^-2,
are checked. The radiation tracker better matches the
manuscript's placement of thermal/Jordan branch formation.
None of these coefficient profiles is an accepted background
history; the simulations test mathematical compatibility only.

None of this predicts the actual primordial
L_ij occupation population.

## 481. An abundance-amplitude consistency bound (not a prediction)

For a single trace mode with

  ell=A_peak cos(m_L t_c) cos(kx),

in a mass-dominated adiabatic regime, the space-time
averaged energy is

  <rho_trace>=K_trace m_L² A_peak²/4,

where K_trace=Bf_L²/3.

Let the *assumed comparison fraction* of the total Jordan
critical density be

  f=<rho_trace>/(3F H_J²).

Then the amplitude needed to produce f is

  A_peak=sqrt(12 f F/K_trace)*H_J/m_L.

For the illustrative choice f=.25, F=K_trace=1,
m_L/H_J=1000, one obtains A_peak≈0.001732.

This demonstrates only that positive energy can coexist
with numerically small relative geometric strain when
the restoring frequency is large. It does not explain
the initial population, show that m_L is below the
effective-theory cutoff, prove the correct lensing or
structure growth, or predict the assumed 25% fraction.

Do not import a CDM value into this expression as a
supposedly predicted SDMC geometric abundance.

## 482. The remaining flatness-versus-gravity consistency tests

A completed SDMC background must solve, from the single
full physical action,

  delta S_total/delta N = 0,
  delta S_total/delta N^i = 0,
  delta S_total/delta h_ij =0,
  delta S_total/delta sigma=0,
  delta S_total/delta L^i_j=0,

with a nonzero physical structural occupancy if the
surviving gravitational mode is claimed to source
the observed extra attraction.

A flat Jordan FLRW average can exist alongside localized
curved wells only when the averaged stress is isotropic
and the total Hamiltonian and momentum constraints
are simultaneously satisfied.

The initial flatness selection also requires a
physical branch-selection mechanism: a lower
selector potential by itself is not a proof of a
time-evolution from curved FLRW to flat FLRW.

The heavy-relative-strain *vacuum* decoupling
(L=0) and an *occupied* heavy structural sector
must not be conflated. An occupied sector can be
mathematically equivalent to an additional cold
massive field in the low-energy description,
even if the fundamental field originates in
the geometry. That physical distinction must be
tested, not merely renamed.

If W switches off while L remains excited, a
consistent late pg calculation must use L kinetic
and restoring/gradient sources, not continue to
claim late gravitational abundance from the
curvature-squared selector.

For a massive oscillating field, small density-
perturbation effective sound speed and CDM-like
clustering require physical wavelengths satisfying
c_L k/(a_J m_L)<<1 and the full coupled field
perturbation calculation; background dustlike
dilution alone is insufficient.

References: SDMC Manuscript B, Part I Sections
34–42, 43–49, 53–56; Manuscript B, Part III
supporting-space ADM action; Cembranos et al.
(JHEP 2016) on averaged oscillating-field
perturbations, DOI:10.1007/JHEP03(2016)013.

Implementation:
  sdmc/pg_flatness_handoff.py
  tests/test_pg_flatness_handoff.py
  .github/workflows/sdmc_pg_flatness_handoff.yml

The numerical output labels all handoff profiles,
mode occupation and energy fractions as conditional
internal model diagnostics, not physical cosmological
likelihood tests.
