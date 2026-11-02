# SDMC — Nonzero-curvature calculation from the Planck-to-radiation handoff

## 1. Scope and physical distinction

The observable flatness target is small |Omega_{k,J}|, not necessarily exact
K_J=0. This continuation tests whether the early thermal handoff's physical
equations compute that magnitude, rather than merely applying archived
kinematic suppression factors to a presumed initial curvature.

No cold-dark-matter fluid, normalization or perturbation fit is reopened.

The exact Planck boundary S_P=1 is prethermal, with
  rho_v,P=rho_P, rho_X,P=0, chi_P=0.
The classical Jordan Hamiltonian does not apply below or exactly at
that boundary unless a microscopic matching prescription establishes its
validity. The following Jordan homogeneous constraint is used only after
the earliest classical start S_A approximately 1.

## 2. Direct nonflat extension of the canonical Jordan constraint

Manuscript B canonical Section 63 uses G2=k1 Z+k2 Z²-V,
G3=g(sigma) Z, F(sigma) R_J/2 and homogeneous
  Z=v²/2, v=dot sigma/N, H_J=dot a_J/(N a_J).

Its scalar energy contribution in the adopted rearrangement is

 rho_X = k1 v²/2 + 3 k2 v⁴/4 + V
         -3 H_J F_sigma v + 3 g H_J v³
         -(g_sigma v⁴)/2.

The exactly homogeneous Jordan FLRW Hamiltonian, allowing
k_J not identically zero, is

  3 F (H_J² + k_J c²/a_J²)
   = rho_b + rho_r + rho_X + rho_UV.

The full nonflat scalar and ADM equations must remain consistent
with this lapse constraint and scalar/thermal conservation equations.

In particular rho_UV is NOT determined merely by renaming
structural density rho_v; it follows from the lapse variation
of the actual independent spatial sector.

## 3. A restricted ultraviolet curvature-energy calculation

Manuscript B Section 54's candidate 3-spatial operator has

  S_UV contains - integral dt N a_J³ W V_K[gamma],
  V_K=eta_K/2 Riem[gamma]².

On an *exactly locked spatial FLRW* slice,
 gamma_ij=M(t)² h_ij,
 a_s=M a_J,
 V_K=6 eta_K k_J²/a_s^4.

If W, eta_K, a_s and the curvature-sector metric are held
independent of the lapse during this restricted variation, the
corresponding UV Hamiltonian density is

  rho_UV,K=6 W eta_K k_J²/a_s^4.

This is only the restricted fixed-lapse form; if the full
structural lapse, mapper and independent 3-metric varied,
additional terms or constraints can arise.

Define the physical curvature contribution (c=1 normalization)
  x=k_J/a_J²,
  lambda=6 W eta_K/M^4,
  E0=rho_b+rho_r+rho_X+rho_UV,other.

Then

  3F(H_J²+x)=E0+lambda x²,

which is the algebraic equation

  lambda x² -3F x + (E0-3F H_J²)=0.     (*)

If lambda>0 and its discriminant is nonnegative, (*) has
two real formal branches:

 x_small=2(E0-3F H²)/
           [3F+sqrt((3F)²-4lambda(E0-3F H²))],

 x_large=[3F+sqrt((3F)²-4lambda(E0-3F H²))]/(2lambda).

The small branch reduces continuously to the ordinary
 x=(E0-3F H²)/(3F)
as W approaches zero. The large branch generally goes as
3F/lambda and may lie far beyond the validity of the spatial EFT.

For an illustrative independent flat-background residual
  3F=1, H=1, E0=1, lambda=.06,
the equation has BOTH
  x=0 and x=16.666667.

This shows that merely adding a positive curvature-squared
potential to the lapse constraint does not uniquely select
a small curvature. A stable physical branch and its initial
data must emerge from the full dynamics.

At E0=1.002 with those other coefficients, the two roots
are approximately 0.00200024 and 16.66466643.
When W=0, the only root is 0.002.

Moreover, for any physical exact-FLRW evolution
  k_J=a_J² x
must remain constant through the established Jordan era.
Solving the constraint independently at each time with
an archived flat background would fail to enforce this
unless the time-dependent scalar and thermal equations
are solved simultaneously.

## 4. A genuine curved FLRW evolution with entropy-producing transfer

To establish what thermal transfer can accomplish by itself,
consider a separate conservative two-component fixed-F GR
classical toy, with NO inflation and NO cold species:

  X'=-3(1+w_X)X - Gamma X,
  R'=-4R+Gamma X,
  Gamma=J/(H X)>=0.

The source cancels from the total density:
  (X+R)'=-3(1+w_X)X-4R.

The reduced equation-of-state test uses exactly Manuscript B's
G2 positive-kinetic form

  w_X=(1+r)/(1+3r),
  r(x)=exp[6(x-1)],
  x=ln a_J.

A specified positive pulse
  Gamma(x)=5 exp[-(x-1.2)²/(2*0.5²)]
transfers X energy to radiation. In a positive-temperature bath,
J>0 corresponds to positive comoving thermal entropy production,
but the pulse is NOT dynamically predicted by the SDMC action.

With 3F=1 and X(0)=1, R(0)=1e-8, the curved
Friedmann constraint is integrated self-consistently:

  H_J²(X,R,K;a_J)=(X+R)-K/a_J².

Different constant K values yield different nonflat expansion
histories, while obeying the same prescribed toy w_X and transfer
per e-fold. K is not assigned the same value in all solutions:
that is exactly the unselected primordial datum being tested.

The exact Friedmann kinematic identity for fixed K is

  dln |Omega_k,J| / dln a_J = 2q_J.

For fixed F without vacuum,
  q_J=[X(1+3w_X)+2R]/(6F H_J²)>0,

so |Omega_k,J| always grows while the classical toy expands,
despite decreasing positive q and radiation entropy production.

In the normalized example with open K=-1e-6:
  at x=0: q≈1.9927, |Omega_k,J|≈9.99999e-7;
  at x=.8: |Omega_k,J|≈1.27668e-5;
  at x=2: |Omega_k,J|≈1.43863e-4;
  at x=3: |Omega_k,J|≈.0010620;
  at x=4: q≈1.0, |Omega_k,J|≈.00779453.

The increase is by about 7794.5, not suppression.
The same qualitative behavior occurs for sufficiently small
closed K>0 before any recollapse; very large positive K
eventually makes H² negative, marking the recollapse/turnaround.

This is a *nonflat classical toy calculation* with real curvature
feedback to H_J. It is NOT a full SDMC nonflat transition: it
does not include the full dynamical F(sigma), G3 kinetic braiding,
the bare structural normalization rho_v, independently solved
N(t), p(t), mapper and curvature-selector UV stress.

Thus it proves a restricted obstruction, not a no-go theorem
for early prethermal structural-to-Jordan curvature suppression.

## 5. Direct interpretation for the actual SDMC transition

Manuscript B's provisional structural curvature measure satisfies

  Omega_k,R ∝ N^-2 for fixed k_R:
  N: 1 -> 2/3 -> 3.41350846,
  |Omega_k,R,rad|/|Omega_k,R,P|≈.0858219.

The Jordan mapping, on corresponding classical FLRW slices,

  Omega_k,J=(k_J/k_R) p² M² Omega_k,R.

This can give very small early Jordan *fractional curvature*
when the early mapper M is tiny, but for fixed k_J the normal
later growth law still applies. Neither relation independently
calculates the primordial k_J amplitude or its probability.

The exact geometric locking condition on a shared FLRW slice
gamma=M(t)²h implies k_R=k_J in a common spatial normalization.
A new variable C_k(t)=k_J/k_R cannot be obtained by simply
replacing it with an arbitrary function while both are exact
fixed-curvature FLRW metrics.

To obtain actual small nonzero final curvature, SDMC needs
a physical handoff derivation that either:
  (a) computes a small nonzero Jordan k_J/curvature radius as
      geometry emerges from the prethermal structural sector;
  (b) demonstrates a nonstandard effective curvature process
      that changes the relevant averaged spatial curvature
      before the final constant-k_J FLRW stage; or
  (c) shows its admissible primordial boundary conditions
      naturally give sufficiently small curvature once
      propagated through the same action.

The candidate positive curvature potential admits a
homogeneous flat minimum, but its effect on physical
nonzero-curvature solutions cannot be read from that
minimum alone. Need full UV lapse+spatial variation,
unambiguous B,W and matter-bath functions, and independent
initial data or boundary selection.

## 6. Test and reproducibility

Module: sdmc/flatness_nonflat_handoff.py
Suite:  tests/test_flatness_nonflat_handoff.py

The tests verify non-flat ADM curvature roots, multiplicity,
positivity and discriminant, two independent fixed-k FLRW
solutions with the same transfer, the exact classical
curvature-growth identity, correct total transfer, lack of
CDM input, and the prethermal limitation.

Scientific verdict: the actual curvature *can be calculated*
only after SDMC provides enough independent physical input
to close and integrate the non-flat early system. The current
Planck-to-radiation background and lapse alone do not
provide a unique absolute curvature. The suppression route
remains open as a spatial-emergence mechanism, but an
ordinary radiation-producing decelerating FLRW handoff
by itself cannot accomplish observable curvature suppression.
