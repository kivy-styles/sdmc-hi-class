# SDMC Part III — Reduced structural trace and Jordan constraint Schur complement

## 1. The precise question

The quadratic trace action has now produced a positive spatially averaged
energy and a Jordan-time reduced energy-exchange identity. Those results do
not yet establish an independently conserved, pressureless geometric source
in the observable Jordan Friedmann equation.

The next algebraic step is to determine what the *actual* quadratic
structural trace term contributes to the Jordan scalar metric constraint
when ell is eliminated. The other gravitational coefficients are those
of the already-present F, G2, G3 scalar–tensor sector, the ADM lapse and
shift, and the material sources. They MUST NOT be invented or calibrated
from omega_c to interpret this calculation.

The physical structural ruler a_s=S/S0 uses R_s=c t_s, whereas cosmic
measurements use Jordan a_J, H_J and time t_c.

## 2. Reduced static trace action

For an individual scalar Fourier mode k the quadratic static energy is

    E= (K omega0²/2)ell² + C (m-Phi_J+ell/3)²,

    K=B f_L²/3,
    omega0²=m_L²+c_L² k²/a_J²,
    C=4 W eta_K k⁴/a_s⁴,
    m=delta ln M, M=(S/S0)/a_J.

Here m is a perturbation of the existing structural mapping, NOT an
additional matter particle. It is related to the existing structural
scalar perturbation and gauge variables by the full theory.

The stationary trace variation gives

    [K omega0²+2C/9]ell = -(2C/3)(m-Phi_J),

    Lambda=2C/(9K),
    ell*= -3 Lambda/(omega0²+Lambda)(m-Phi_J).

Substituting ell* into the quadratic energy gives the exact restricted
Schur complement:

    E_eff = Kgeom(k) (m-Phi_J)²,

    Kgeom = C omega0²/(omega0²+Lambda).

The trace denominator corresponds to a reduced scalar restoring
frequency omega_eff²=omega0²+Lambda. Static elimination of an evolving
mode is a justified approximation only when the mode evolves much
faster than H_J AND faster than the time variations of its drivers.
The code provides the necessary frequency/Hubble ratio, not a
sufficient proof of full quasi-static validity.

## 3. Infrared and ultraviolet structural response

For positive K and positive W eta_K, and assuming c_L²>0:

  if k -> 0 and omega0² remains positive (including a nonzero mass),

    Kgeom(k) = [4 W eta_K/a_s⁴] k⁴ + O(k^6).

  if k -> infinity,

    Kgeom(k) = (9K/2) c_L² k²/a_J² + subleading corrections.

The long-wavelength metric response from this *restricted* static
relative-strain term is therefore derivative suppressed (k⁴),
whereas the high-k response approaches a k² stiffness. Neither
limit behaves by itself like a k-independent CDM density source.

If the mass and gradient terms both vanish, degeneracy must be treated
with care: a zero-frequency trace constraint can cancel the relative
strain and the static elimination is not an autonomous propagating
mode. k=0 also has no curvature-squared response.

## 4. Embedding in a general Jordan scalar constraint

Let x denote the remaining Jordan scalar variables, e.g.

    x=(Phi_J, Psi_J, delta sigma, delta N_ADM, beta_ADM, ...).

The full action would give a (possibly indefinite or singular,
pre-gauge-fixing) quadratic scalar kernel H_Jordan(k,t).
The structural trace adds

    Delta E = Kgeom(k) (m-Phi_J)².

On a slice where u=m-Phi_J=v^T x is an appropriate linear combination,
the updated Hessian is

    H_full = H_Jordan + 2 Kgeom v v^T.

With the full action derived and after a valid choice of physical
variables/gauge, the remaining source equation is

    H_full x = J_matter + J_other.

The code implements a generic Schur complement on a supplied real
symmetric nonsingular eliminated block. It does NOT generate
H_Jordan or claim to fix gauge degeneracies.

This is how the supporting-space action would enter the Jordan
constraints once the actual ADM and F/G2/G3 kernels are supplied.

**Important:** the absence of a first-order independent geometric
density perturbation about an exactly unexcited UV background
does not imply the quadratic action cannot alter *linear metric
equations*. Varying a quadratic action gives linear equations.
It does mean this calculation has not supplied the nonzero cold
background density and cold fluid perturbation of a separate
physical rho_g sector.

## 5. Fully coupled but schematic positive-stiffness control

To verify the algebra without pretending to know H_Jordan, define a
positive quadratic diagnostic

    E_toy= (A/2) Phi²+(D/2)m²
          +(K omega0²/2)ell²+C(m-Phi+ell/3)² -J Phi,

    A>0, D>=0.

A is a positive *toy* potential stiffness, not the physical signed
gravitational constraint. D is a toy stiffness of the existing
structural mapping/scalar, not an introduced particle.

After eliminating ell:

    E_toy,red = (A/2)Phi²+(D/2)m²+
                    Kgeom(m-Phi)²-J Phi.

The algebraic mapper equation gives

    m= [2Kgeom/(D+2Kgeom)] Phi.

Then the observable-potential toy equation reads

    [ A + 2Kgeom D/(D+2Kgeom) ] Phi = J.

The toy relative response to its unmodified version is

    R_toy = A/[A+2Kgeom D/(D+2Kgeom)] <= 1.

In the D->0 limit, the mapper freely follows Phi, u=0, and
the extra stiffness vanishes. In the D->infinity limit,
m is stiff, giving additional stiffness 2Kgeom.

Therefore this *positive-stiffness toy block* can screen response
but cannot create the ~6.57 amplitude of missing gravitational
attraction by itself. It would be INVALID to generalize this result
to the physical, indefinite and constrained Jordan gravitational
system without deriving its actual Hessian.

## 6. Checks implemented

sdmc/pg_constraint_schur.py includes:

- reduced scalar-trace kernel and static frequency ratio;
- full toy 3x3 static action Hessian [Phi,m,ell];
- exact analytic and numerical elimination of ell;
- direct simultaneous solution for all three variables;
- generic Schur elimination for a supplied full symmetric Hessian;
- asymptotic k⁴ and k² checks;
- mapper-free and stiff-mapper limiting responses.

Tests are in tests/test_pg_constraint_schur.py.

Run:

    python3 -m unittest discover -s tests -p 'test_pg_constraint_schur.py' -v
    python3 -m sdmc.pg_constraint_schur

The dimensionful interpretation of k, W eta, and K must be
consistent with a chosen physical unit system before comparison
to observed modes in h/Mpc.

## 7. Exact missing inputs for physical closure

To turn the algebraic Schur complement into the SDMC no-CDM physical
metric response, the following must be derived from the *complete*
accepted action, not from an effective reference-CDM Poisson source:

(1) the quadratic ADM lapse/shift and mixed metric Hessian;
(2) the active scalar F/G2/G3 perturbation/mixing terms;
(3) how the structural mapping m follows delta S and the metric;
(4) the curvature-activation W and coefficient history;
(5) primordial inhomogeneous structural initial conditions and
    averaged homogeneous geometric stress;
(6) gauge-invariant propagating modes, stability, sound speeds;
(7) the observable potentials Phi_J and Psi_J from the closed
    gravitational constraints.

The previous quadratic 3D averaged structural source and energy
exchange remain algebraic/internal-action results; this branch
does NOT yet simulate a physically predicted rho_g, a nonlinear
halo mass function, CMB anisotropies or weak-lensing spectra.
