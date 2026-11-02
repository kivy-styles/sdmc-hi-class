# SDMC flatness — suppression branch quantitative audit (2026-11-02)

## Status and precise question

Flatness selection at K=0 is NOT a required observable target:
CMB+distance measurements admit sufficiently small nonzero
present Jordan curvature. This continuation tests the
suppression mechanism instead and does not reopen CDM.

The parameters below come from *Manuscript B Part I, §28
and the flatness continuation §31–33*, not from a newly
integrated non-flat late300 solution.

The Planck boundary is prethermal; the curvature fraction at
the exact boundary is used as a counterfactual structural
diagnostic, not a classical measurable Jordan FLRW fraction.

## 1. The reduced two-stage structural suppression

Assuming a fixed nonzero structural curvature label k_R,
the provisional structural fraction depends on the lapse:

  Omega_k,R=-k_R R0²/N².

The archived lapse sequence N: 1 -> 2/3 -> 3.41350846
therefore yields

  Planck to crossover: (1/(2/3))²=2.25  (grows);
  Crossover to radiation: ((2/3)/3.41350846)²
                        ~=0.0381431 (suppressed);
  Full sequence: (1/3.41350846)²
                        ~=0.0858219 (suppressed).

The result is 91.4178% reduction in the magnitude of this
*provisional structural curvature fraction*. It is NOT
an 8.58% curvature of space and NOT the Jordan observable.

## 2. The archived Jordan-frame mapping

From homogeneous background kinematics,

  Omega_k,J=-k_J/(a_J² H_J²),
  Omega_k,R=-k_R/(a_s² H_R²),

and the bridge variables are

  a_s=M a_J, H_R=p H_J,
  p=dln a_s/dln a_J.

Thus

  Omega_k,J=C_k p² M² Omega_k,R,
  C_k=k_J/k_R.

This relation is algebraically exact given the definitions,
but does not alone compute k_J for an actually non-flat
structural-to-Jordan handoff.

The earlier flat late300 reference history gave, for C_k=1,

  freeze-out Jordan/initial structural: 6.35e-17,
  recombination Jordan/initial:          2.10696e-4,
  present Jordan/initial:                0.0917953.

These are CONDITIONAL diagnostic factors using the
accepted background H and mapping; curvature backreaction
for a non-flat solution has not been solved.

The implied recombination-to-present relative growth is
0.0917953/0.000210696=435.67652.
The recombination-era inferred reference Hstar/H0
=1091 sqrt(435.67652)=22772.2964.

Early small Omega_k can therefore later be substantially
larger when aH changes. The small mapped early number
is not by itself evidence of a lasting new physical
curvature-damping mechanism.

## 3. What near-flatness requires without exact K=0

Using illustrative magnitude target |Omega_k,J,0|<0.001,

  0.0917953 |C_k,0 Omega_k,R,P| <0.001,

or

  |C_k,0 Omega_k,R,P| <0.01089380393.

For C_k,0=1, the initial provisional structural
curvature fraction must have magnitude <~0.01089.

For |Omega_k,R,P|=1, C_k,0 must have magnitude <~0.01089
(a further factor 91.7953). Neither requirement
can be declared dynamically natural or unnatural until
the primordial initial-condition statistics and actual
curved handoff equations are solved.

For |Omega_k,R,P|=.01, C_k,0=1 gives
|Omega_k,J,0|=.000917953 (meets example target).
No exact zero is necessary.

Do not elevate the illustrative 0.001 numerical target
to an exact experimental detection or absolute threshold.

## 4. New consistency result: exact locking does NOT evolve C_k

The manuscript's candidate microscopic branch has
  gamma_ij=M(t)² h_ij  (L=0).

For corresponding constant-curvature spatial FLRW
metrics and a common normalized coordinate chart,
  (3)R[gamma]=M^-2 (3)R[h].

Writing the two spatial Ricci scalars gives
  6 k_R/a_s²=M^-2 6 k_J/a_J².

Since a_s=M a_J, this implies
  k_R=k_J and hence C_k=1, *at all epochs of exact
  homogeneous conformal locking in the shared chart*.

More generally, if k_R and k_J are each independently
constant FLRW curvature labels, C_k=k_J/k_R is also
constant even if it is not one due to differing chart
normalizations. Therefore

  A_k=dln|C_k|/dln a_J = 0.

The manuscript's expression
  dln|Omega_k,J|/dln a_J=2q_J+A_k
is an algebraic extension. A nonzero A_k requires
a physically specified spatial transition that leaves
a model with two fixed-k FLRW metrics, or changes the
definition of effective averaged spatial curvature.
It is not produced by the time-dependent mapper M(t)
or the structural clock N(t) alone.

This is a diagnostic obstruction to invoking a freely
chosen late-time C_k~0.01 while retaining exact locking
of two fixed-k FLRW spatial geometries. A *microscopic
non-FLRW or nonconformal phase before/during emergence*
could conceivably select the small nonzero k_J that
freezes after formation. It must be DERIVED, not fitted.

The positive curvature-selector action
  V_K=eta_K Riem[gamma]²/2
has a quadratic minimum at k=0; having that minimum does
not compute a continuous evolution of a global FLRW
curvature label, which is a fixed parameter of standard
FLRW spatial geometry.

## 5. No-CDM scope and next physical calculation

Keep second-route geometric-gravity/extra-mass case closed.
This exercise asks solely whether an early admissible
structural curvature maps into small |Omega_k,J,0|.

The actual missing calculation is a non-flat solution of
the *same* 4D covariant scalar–tensor + supporting-space
field equations across early handoff, with K!=0 allowed
but not prescribed to vanish. It must supply:
  - material/radiation/active-scalar transfer with
    self-consistent curvature backreaction;
  - the actual inhomogeneous or nonconformal spatial map
    determining k_J and the physical curvature parameter;
  - endpoint propagation to today for a family of
    initial structural curvature values;
  - stable local perturbations and observational metric
    constraints without using the FLRW flat solution as
    an exact non-flat solution.

There is no honest new dynamical proof of sufficient
suppression until those equations and boundary conditions
are supplied.

## 6. Code and tests

Module: sdmc/flatness_suppression_audit.py
Test:   tests/test_flatness_suppression_audit.py

Tests cover:
 - full two-stage 2.25 * 0.0381431 factor;
 - very small early Jordan curvature and later growth;
 - conditional target |C_k0 Omega_kR,P|<0.0108938;
 - exact bridge identity with nonzero kJ;
 - C_k fixed for two exact FLRW metrics with fixed k;
 - constant C_k=1 under pure time-only conformal locking;
 - safe refusal to claim curved background self-consistency.
