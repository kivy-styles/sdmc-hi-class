# SDMC Manuscript B, Part III — Primordial Structural-Well Source Closure

## 492. The actual direct material-to-strain variation

Manuscript B, Part III Eq (106.1) displays a supporting-space
action with positive trace/traceless kinetic, gradient and mass
terms, and a transient curvature-selector term. The ordinary
material and radiation actions are minimally coupled to the
physical Jordan metric, as in Manuscript B and the older
SDMC matter-loaded structural-release investigation.

Write the displayed classical action schematically as

 S_total =
   S_{F,G2,G3}[g_J,sigma]
   +S_b[g_J,psi_b]+S_r[g_J,psi_r]
   +S_support[g_J,sigma,M,L].

Holding the independent Jordan metric and the material
configuration fixed, minimal matter coupling means

 delta S_b/delta L = delta S_r/delta L = 0.

The F/G2/G3 sector contains sigma and the Jordan metric but no
explicit independent relative strain L, hence its direct
partial derivative with respect to L is also zero.

This is a DIRECT off-shell variation. Material is still
coupled indirectly through Jordan geometry when the
complete metric/mapper/strain constraints are solved.

The only explicitly displayed linear metric-to-relative-trace
coupling in Section 106 is the W-dependent curvature mismatch.
No additional direct matter-to-L interaction may be quietly
introduced to make the primordial well persist.

## 493. Exact source-channel Hessian from accepted quadratic action

For one real cosine mode define
  u=m-Phi_J,
  K=B f_L²/3,
  omega0²=m_L²+c_L²k²/a_J²,
  C=4W eta_K k⁴/a_s⁴.

Manuscript B Eq (106.3) gives the positive reduced
static supporting energy

 E_support = K omega0² ell²/4
         + C (u+ell/3)²/2.

To audit the bilinear coupling, add a purely diagnostic
positive Jordan metric stiffness Q and an abstract
matter-to-metric source j_u:

 E_audit=Q u²/2-j_u u+E_support.

The Q>0 energy is NOT the actual physical signed
ADM metric Hessian; it is only a reproducible channel
test. Its exact Hessian in (u,ell) is

 H =
  [ Q+C,               C/3 ]
  [ C/3,  K omega0²/2+C/9 ].

Thus the direct metric-to-trace mixing is

  H_(u,ell)=C/3=4W eta_K k⁴/(3a_s⁴).

It is proportional to W and vanishes exactly at
W=0 or k=0. The material source j_u is present
only in the metric channel, not in the ell channel.

Solve H [u,ell]^T=[j_u,0]^T. Eliminating ell gives

  ell_eq=-3 Lambda u/(omega0²+Lambda),
  Lambda=2C/(9K),

  (Q+K_geom) u=j_u,

  K_geom=C omega0²/(omega0²+Lambda).

The geometric metric stiffness, as well as the static
source-supported ell, disappears as W→0.

The diagnostic positive Q channel typically screens,
rather than enhances, u; the physical gravitational
Jordan Hessian is SIGNED and its potential/slip
response cannot be inferred from this positive toy.

As a physical subcase cross-check, the earlier
Jordan F(sigma)R plus canonical scalar quasistatic
constraint Hessian is exactly recovered when W=0.
That earlier subcase does not include rolling G3
and cannot stand in for the complete late300 model.

## 494. A static positive-energy integral obstruction at W=0

When W=0, K>0, m_L²>0 and c_L²>=0, and no
additional source acts on ell, its static Euler-Lagrange
equation on a fixed Jordan slice is

  m_L² ell - (c_L²/a_J²) Delta ell=0.

Multiply by ell and spatially integrate over a
periodic domain, or a localized domain with
vanishing surface flux:

  integral [m_L² ell²
    +(c_L²/a_J²)|grad ell|²] d³x=0.

Both terms are nonnegative; for m_L²>0,
the only regular finite-energy solution is

  ell=0.

The same positive-coercive argument applies to
each component of the traceless strain Q when the
spatial operator and mass are positive and there
is no source or nonlinear soliton interaction.

This is stronger than checking the local constant
mode's equilibrium: there is no STATIC localized,
unforced relative-strain lump in the healthy
positive quadratic W=0 theory under these boundary
conditions.

It does NOT rule out:
 - ordinary Jordan gravitational wells with curved h_ij;
 - excited TIME-DEPENDENT structural L modes with rho~L²;
 - nonzero late direct coupling from a separately
   specified and verified parent action;
 - nonlinear bound states or topological sectors;
 - degenerate massless zero modes, which fall outside
   the strictly positive-mass lemma.

The code independently checks the integrated field
times Euler-Lagrange residual against the positive
mass-plus-gradient norm on a periodic 3D Fourier grid.

## 495. What remains of the primordial well when L=0

A locked structural metric obeys gamma_ij=M² h_ij.
It may retain nonzero intrinsic curvature if the
physical Jordan metric is inhomogeneously curved.

For scalar modes,
  q=m-Phi_J+ell/3,
  delta R_s=4 k² q/a_s².

When W=0 and ell_eq=0,

  delta R_s=4 k²(m-Phi_J)/a_s²

may be nonzero although the selector energy
W eta_K Riem²/2 is zero.

Therefore:
 - LOCAL gravitational curvature survives in the
   ordinary locked Jordan system.
 - The separate STATIC relative-strain distortion
   is absent from the strictly positive quadratic
   W=0 source-free branch.
 - Nonzero Jordan curvature is not an additional
   independent gravitational energy or mass by itself.

## 496. A genuine surviving scalar curvature-loading channel

The accepted covariant scalar-tensor action contains

 S_J ⊃ integral sqrt(-g) [F(sigma)/2] R_J.

Its scalar Euler-Lagrange equation contains the
direct nonminimal curvature contribution

  delta S_J/delta sigma ⊃ F_sigma R_J/2.

This term can be nonzero even when W=0, since
matter can make the Jordan curvature nonzero.

It is a sigma-channel, not a direct L-channel.
It must be combined with the full G2/G3 scalar
equation and the Jordan gravitational constraints.

The archived late300 No-Slip scalar–tensor sector
has its own constrained gravitational response:
the mere presence of F_sigma R does not prove
an additional CDM-sized gravitational contribution,
a new independent material source, or required
early clustering.

The code verifies the F_sigma R/2 contribution
by independent central differences of F(sigma)R/2.

## 497. Existing excitation versus matter-driven static response

The coefficient B(sigma) or f_L²(sigma) can vary
with the rolling scalar. Around an excited
L_0 !=0, even the mass/restoring energy

 E_L=K(sigma) omega0² L²/4

has a nonzero mixed sigma–L Hessian:

  d² E_L/(d sigma dL)
     =K_sigma omega0² L/2.

At exactly L=0, this vanishes.

Consequently a nonzero INITIAL structural mode
can dynamically couple to sigma and Jordan gravity
after W shuts off, but the existing action's
quadratic kinetic/mass terms do not CREATE such
a mode from an exactly unexcited state.

The supported energy is positive and quadratic,
rho_L∝dot L²+m_L²L²+c_L²|grad L|²/a_J².
Its abundance and primordial correlations are
not supplied by the bare normalization
rho_v,P=rho_P at S_P=1.

The material ledger
  M_u=M_remaining+M_inferred
does not furnish an independent canonical
occupation number for L or a fixed geometric
energy fraction.

## 498. Why adding a new direct L source is a new physical assumption

A thought experiment with an extra static term
 -g_new j_m ell
would yield ell_eq=2 g_new j_m/(K omega0²) at W=0.

That elementary identity shows which type of
additional operator is required to maintain a
source-supported relative-strain depression.
It is NOT contained in the displayed action
and has NOT been introduced into the physical
SDMC field equations or interpreted as an
observational result.

A proposed direct matter/relative-strain coupling
would have to be constructed from the complete
diffeomorphism-invariant parent action, respect
the universal Jordan matter coupling and matter
conservation, and pass strong-field, equivalence
and cosmological constraints. It cannot simply
be tuned to supply the missing CDM density.

## 499. What is actually determined now

Established on the stated classical branch:
  1. Minimal Jordan material does not directly
     source L in the accepted displayed action.
  2. The linear curvature-mismatch coupling is
     exactly proportional to W k^4/a_s^4.
  3. The static Schur response disappears when
     W=0, and the signed canonical Jordan
     quasistatic subcase returns to its
     scalar-only gravitational response.
  4. The healthy unforced positive massive
     spatial operator has no static nonzero
     finite-energy relative-strain lump.
  5. A locked Jordan primordial gravitational
     well may still have local metric curvature.
  6. F_sigma R_J/2 can continue to source the
     scalar even when the relative strain no
     longer has direct material sourcing.
  7. A PRE-EXCITED L can carry positive energy
     and interact quadratically with sigma,
     but cannot be equated with the original
     source-supported primordial depression.

Not established:
  - full 4D mapper/triad action and any linear
    matter–structural operator that might survive W=0;
  - independent, healthy primordial initial
    structural excitation state;
  - a separate geometric background abundance
    with conserved dustlike stress;
  - full rolling G3 plus excited structural ADM
    Poisson, Weyl, growth and CMB transfer functions.

The theory's dark-matter replacement question
therefore remains open. But the source-channel
audit eliminates a specific unsupported inference:
retaining an original Jordan gravitational well
does not automatically retain an independent
static structural gravitational source.

Main code:
  sdmc/pg_primordial_source_action_audit.py
  tests/test_pg_primordial_source_action_audit.py
  .github/workflows/sdmc_pg_primordial_source_action_audit.yml

Manuscript sources: Part III Sections 98,101,106,
111,112,114,115, and the earlier minimally coupled
Jordan matter sector. All toy positive-Q and
new-operator examples are expressly auxiliary,
not physical late300 predictions.
