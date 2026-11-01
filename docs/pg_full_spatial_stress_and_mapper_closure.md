# SDMC Manuscript B — Part III: Excited Three-Dimensional Spatial Stress, Structural–Jordan Mapper, and Pressureless Closure

## 466. Objective and scope

The preceding parts derived the covariant rolling F/G2/G3 ADM kinetic
constraints, coupled them to the supporting-space scalar trace,
calculated reduced energy and momentum sources from excited wells,
and verified a cubic spatial-metric source from the exact
conformally flat structural Riemann invariant.

This continuation calculates off-shell ADM *spatial metric source
projections* for the full symmetric 3D relative strain, including
the trace ell and traceless tensor Q_ij. It also resolves why
different restricted pressure calculations of the same curvature
term give different values.

The accepted supporting action is (Manuscript B Part III, Section 106)

  S_g = ∫ dt_c d³x N_ADM sqrt(h)
    { (B f_L²/2)[(D_n L)²-c_L²(DL)²-m_L²L²]
       -(W eta_K/2) R_ijkl[gamma]R^ijkl[gamma] }.

The structural scale a_s=S/S0 belongs to the structural clock t_s;
ADM velocity and observable H_J use Jordan time t_c. The mapper
M=a_s/a_J connects metrics; its complete 4D transformation and
equation of motion are not yet established.

The displayed operator class specifies gamma in a locked-flat
reference frame but not a unique general nonlinear 4D tensor
completion for all off-shell h_ij and Q_ij configurations.
Consequently the metric-congruence continuation below is a
*defined positive-definite spatial completion*, not a uniqueness
theorem for the accepted physical SDMC action.

## 467. Positive-definite trace/traceless metric completion

Decompose the mixed relative strain into

  L^i_j=ell delta^i_j/3 + Q^i_j,
  Tr Q=0,  Q_ij=Q_ji on the reference spatial slice.

Let epsilon_ij(t_c,x) denote an arbitrary symmetric Jordan
spatial metric strain, with

  h_ij = a_J² [exp(2epsilon)]_ij.

For this spatial completion choose

  gamma = a_s² exp(2ell/3)
       exp(epsilon) exp(2Q) exp(epsilon).

For symmetric epsilon and Q, gamma remains symmetric positive
definite, and at epsilon=Q=0 it reduces to
gamma_ij=a_s² exp(2ell/3)delta_ij.
For Q=0 and epsilon=zeta I,
gamma_ij=a_s² exp[2(ell/3+zeta)]delta_ij,
exactly the conformal metric of the previous calculation.

This completion keeps the structural ruler, Jordan metric,
and structural relative shape logically separate. It is
the action to be varied, not an effective CDM density ansatz.

## 468. Exact 3D curvature of the traceless-excited state

For any symmetric positive spatial gamma_ij(x), calculate

  Gamma^k_ij[gamma]
    =gamma^kl/2 (
       partial_i gamma_jl+partial_j gamma_il
       -partial_l gamma_ij),

  R_ij =
    partial_k Gamma^k_ij
    -partial_j Gamma^k_ik
    +Gamma^k_ij Gamma^l_kl
    -Gamma^k_il Gamma^l_jk.

In THREE spatial dimensions the Riemann tensor contains no
independent Weyl tensor. The exact identity is

  R_ijkl R^ijkl
    =4R_ij R^ij - (R^(3))².

The numerical module evaluates the Ricci tensor and the invariant
using periodic three-dimensional Fourier derivatives, including
nonlinear gamma and its inverse. A separate test confirms that
the zero-traceless-Q limit reproduces the earlier explicit
conformal Riemann squared formula.

The exact spatial curvature operator contains fourth derivatives
of gamma and produces higher-gradient gravitational effects.
These do not by themselves define an extra particle species.

## 469. Spatial metric stress from finite off-shell variation

At N_ADM=1, N^i=0, with the background spatial metric held fixed
during the single-epoch action evaluation, define the exact
spatially averaged reduced action per Jordan background volume
a_J³ as

  Lbar(epsilon) = <
    exp(Tr epsilon) [
       (B f_L²/2) Tr(dot L²)
      -(B f_L² m_L²/2) Tr(L²)
      -(B f_L² c_L²/2) h^ij Tr(partial_i L partial_j L)
      -(W eta_K/2) Riem[gamma]² ] >.

The spatial-gradient contraction is written for the fixed
reference-frame matrix L; connections associated with a
fully time-dependent orthonormal frame and the Lie transport
of its traceless part remain to be derived from the parent action.
This is therefore a partial metric stress, with a clearly stated
reference-frame convention.

An arbitrary off-shell metric test strain epsilon_ij(x) gives

  Delta Lbar =
      (delta S_g/(a_J³ V delta epsilon)) · delta epsilon.

At the reference h=a_J²delta,

  delta Lbar/delta epsilon[basis]
   = <T^i_j basis^j_i>

within this supporting-space metric completion, with T
the source defined by Jordan spatial-metric variation.
No explicit formula for a cosmological \rho_g or P_g is
assigned before the variation.

The code computes all SIX independent spatial tensor
projections by central differences of the exact action:
 three diagonal projections T^x_x, T^y_y, T^z_z,
 and three symmetric off-diagonal projections T^x_y,
 T^x_z, T^y_z.
It also accepts spatially varying test strains, giving
inhomogeneous shear-source projections without assuming a
homogeneous isotropic stress.

The trace gives a restricted mean isotropic stress; the
diagonal spread and off-diagonal projections quantify the
anisotropic stress. Vanishing average anisotropic projection
for a symmetric ensemble does not imply local isotropy.

## 470. A decisive pressure comparison for the curvature operator

For an arbitrary structural trace/traceless state, the reduced
action energy components on the reference slice are

  Ekin  = <B f_L²/2 Tr(dot L²)>,
  Emass = <B f_L² m_L²/2 Tr(L²)>,
  Egrad = <B f_L² c_L²/2 h^ij Tr(partial_iL partial_jL)>,
  Ecurv = <W eta_K/2 Riem[gamma]²>.

Under an ISOTROPIC Jordan strain epsilon_ij=epsilon delta_ij,
and holding the mapper M and structural L fixed:

  sqrt(h) -> exp(3epsilon)sqrt(h),
  h^ij -> exp(-2epsilon)h^ij,
  gamma_ij -> exp(2epsilon)gamma_ij,
  Riem[gamma]² -> exp(-4epsilon)Riem[gamma]².

Thus the four action terms scale as

  Kinetic:  exp(3epsilon) Ekin,
  Mass:     exp(3epsilon) Emass,
  Gradient: exp(epsilon) Egrad,
  Curvature:exp(-epsilon) Ecurv.

Spatial ADM variation yields, EXACTLY for this fixed-mapper
completion,

  P_mapper = Ekin - Emass - Egrad/3 + Ecurv/3.   (470.1)

In particular the curvature operator alone has the
*isotropic spatial metric stress* contribution

  P_curv,mapper = +Ecurv/3.                     (470.2)

This is a direct variation of the nonlinear spatial Riemann
squared action, independently verified by computing exact
inhomogeneous 3D curvature for varied metrics.

The earlier local oscillator study used a DIFFERENT
restricted Jordan scale comparison: varying a_J at
fixed structural ruler a_s, so that M=a_s/a_J
compensated the metric scale and gamma remained fixed.
That diagnostic produced

  P_ruler=Ekin-Emass-Egrad/3-Ecurv.             (470.3)

The difference is

  P_mapper-P_ruler=(4/3)Ecurv.                 (470.4)

The two calculations hold different structural variables
fixed. One must NOT identify either one on its own with
the final gauge-invariant total physical SDMC P_g.

The complete mapper/structural scalar equation, rolling
Horndeski lapse/shift, and the structural curvature/strain
operator need to be varied together to settle the
physical pressure splitting and exchange.

## 471. Consequence for continuity and the geometric dust question

For the same evolving local reduced energy E, define two
bookkeeping exchanges

  J_ruler = dot E+3H_J(E+P_ruler),
  J_mapper= dot E+3H_J(E+P_mapper).

Using (470.4),

  J_mapper=J_ruler+4H_J Ecurv.                 (471.1)

For a single reduced trace oscillator in an expanding Jordan
era, with

  Omega²=m_L²+c_L²k²/a_J²+Lambda,
  Lambda=2C/(9K), C=4W eta_Kk⁴/a_s⁴,

the action-derived on-shell Ward identity verifies this
relation for arbitrary initial trace amplitude, velocity
and the prescribed coefficient evolution.

At high frequency and slow coefficient variation
(phase-averaged WKB), define G=c_L²k²/a_J² and
Omega²=m_L²+G+Lambda. The fixed-mapper restricted
pressure and exchange satisfy

  w_mapper =
       (G+2Lambda)/(3Omega²),

  J_mapper/(H_J E)=
    [(w_run-4p_v-k_run+4)Lambda]/(2Omega²),

where

  p_v=dln S/dln a_J,
  k_run=dln K/dln a_J,
  w_run=dln(W eta_K)/dln a_J.

The logarithmic dilution identity is

  dln E/dln a_J =
      -3(1+w_mapper)+J_mapper/(H_J E)

    = -3 +
       [-2G+(w_run-4p_v-k_run)Lambda]/
       (2Omega²).

Therefore the energy dilution prediction of the previous
adiabatic investigation remains intact; only the allocation
between a restricted stress and an explicit exchange term
depends on which quantities are held fixed in the off-shell
variation.

For constant couplings and CURVATURE-DOMINATED WKB:

  Ekin≈Ecurv≈E/2,
  w_mapper≈2/3,
  J_mapper/(H_J E)≈-2(p_v-1),
  E~a_J^(-3-2p_v).

Thus:
  matter tracking p_v=3/2: w_mapper=2/3, J/(HE)=-1,
                           dln E/dln a_J=-6;
  radiation tracking p_v=2: w_mapper=2/3, J/(HE)=-2,
                            dln E/dln a_J=-7.

The reduced curvature-dominated mode therefore is NOT
an independently conserved pressureless source in this
specific fixed-mapper spatial variation.

Even imposing the conditional coefficient evolution
w_run=4p_v+k_run to maintain Omega constant gives

  E~a_J^-3,
  w_mapper≈2/3,
  J_mapper/(H_J E)≈+2.

That is a cancellation between positive spatial pressure
and energy transfer, NOT separately conserved dust.

The mass-restoring adiabatic trace limit, with negligible
gradient and curvature fractions, instead yields
w≈0, J≈0 and E~a_J^-3. However, neither its primordial
abundance nor its observed gravity/lensing response is
derived by this reduced oscillator calculation.

## 472. Inhomogeneous shear and the FLRW isotropy condition

An excited traceless Q_ij contributes energy even when the
trace ell has zero spatial average. Its quadratic kinetic,
mass and gradient contributions, as well as its exact
supporting curvature, enter the mean ADM Hamiltonian source.

An anisotropic structural state can give
<T^x_x> != <T^y_y> or an off-diagonal <T^x_y>,
as confirmed by the six stress projections.

For a valid *homogeneous and isotropic averaged Jordan FLRW*
background, the average traceless stress must be negligible
or canceled dynamically:

  <pi^i_j> = <T^i_j -delta^i_j T^k_k/3> ≈0.

This does not require every local gravitational well to
have zero shear, but it does constrain the ensemble
or averaging operation. Equal-axis configurations
can have isotropic averages while local shear remains.

The ADM off-shell stress calculation alone does NOT
produce a physical causal closure for such an ensemble.
The initial structural power spectrum, covariance,
evolution, and transfer through the scalar–Jordan
constraints remain required.

## 473. Formal background and primordial normalization closure

To obtain one physical no-CDM geometry:

  3F H_J² = rho_b+rho_r+rho_X+bar rho_g,

must arise from the combined ADM Hamiltonian constraint
*with the same F/G2/G3 scalar and supporting-space action*.

An excited structural well gives a nonzero reduced source,
but there is no unique prediction for its amplitude from
the quadratic action coefficients alone. The energy changes
quadratically when the initial excitation amplitude changes.
The structural radius and M_u shell-volume identity do not
set the required canonical initial mode amplitudes.

One must determine the physical background state and
initial spectrum without importing omega_c or the
conventional CDM matter transfer function as a target.

Baryon/material bookkeeping and geometric excitation energy
must remain distinct to avoid double counting:
  M_u = M_remaining + M_inferred
is a conserved material inventory identity, not itself
an extra source of geometric excitation energy.

The exact Planck boundary has rho_v,P=rho_P at S=1,
rho_X,P=0 and no classical pg energy spectrum specified.
The classical action begins only on the post-Planck
classical branch; no unreliable t<t_P extrapolation is used.

## 474. Tested implementation and limitations

Modules:
  sdmc/pg_full_spatial_stress.py
  sdmc/pg_mapper_pressure_accounting.py

Tests:
  tests/test_pg_full_spatial_stress.py
  tests/test_pg_mapper_pressure_accounting.py

New checks include:
- exact general 3D metric Ricci/Riemann scalar identity;
- analytic conformal 3D curvature comparison;
- positive metric and nonzero traceless-strain curvature;
- six independent homogeneous spatial stress projections;
- an arbitrary INHOMOGENEOUS shear projection;
- direct metric variation of all four energy terms;
- +Ecurv/3 fixed-mapper isotropic curvature pressure;
- exact difference between fixed-mapper and fixed-ruler stress;
- exact on-shell Ward exchange identity;
- matter/radiation WKB curvature pressure and transfer;
- compensated a_J^-3 energy with nonzero pressure/transfer;
- arbitrary positive mixture and free amplitude degeneracy;
- prior 185 structural and scalar-tensor regression tests.

The numerical 3D Ricci/Riemann implementation uses periodic
spectral differentiation; arbitrary nonlinear spatial fields
need resolution/convergence tests to control Fourier aliasing.

The present result has NOT derived the complete parent
spacetime tensor metric action, its full covariant connection
for time-dependent traceless Q, the coupled inhomogeneous
F/G2/G3 lapse/shift constraints, or a physical FLRW
no-CDM cosmological abundance.

Therefore it establishes *action-level spatial metric
stress projections of a clearly defined nonlinear
supporting-space completion* and the pressure/exchange
obstruction of the curvature WKB branch, not the final
gauge-invariant total T^mu_nu or observational victory
over the standard cosmology.

The mathematically remaining task is to specify the
spacetime parent action for gamma_ij, M and L_ij,
derive all Euler-Lagrange/ADM constraints on a nonzero
inhomogeneous background, and then compute physical
mu(k,a), Sigma(k,a) and cosmological likelihoods
without imported CDM normalization.
