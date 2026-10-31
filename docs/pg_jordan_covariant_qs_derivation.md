# SDMC Part III — Covariant Jordan quasistatic constraint checkpoint

## Scope: a real subcase of the accepted covariant operator class

The manuscript's classical covariant action is

  S = integral sqrt(-g)[G2(sigma,Z)-G3(sigma,Z) Box sigma+
                         F(sigma) R/2]+S_b+S_r,

  G2=k1(sigma)Z+k2(sigma)Z²-V(sigma),
  G3=g(sigma)Z, Z=-g^ab(d_a sigma)(d_b sigma)/2.

The **physical late300 No-Slip sector has G3 and rolling sigma**.
The calculation here is an independently verifiable **restricted
quasi-static checkpoint** on the g=0, frozen-background-gradient,
linear-in-Z subfamily. Its static quadratic action is derived from
the actual covariant F R/2 + k1 Z - V subclass, supplemented by
the previously derived *reduced* structural trace-eliminated
spatial term. It does not silently assign a gravitational strength
or cold energy density to the supporting-space geometry.

Structural quantities obey R_s=c t_s, sigma=S/S0, a_s=S/S0.
Jordan observables use a_J, t_c and the potentials Phi_J/Psi_J.
With Newtonian-gauge background a_J fixed, the perturbation of
ln M (M=sigma/a_J) is m=delta sigma/sigma. This is a fixed-gauge
identity and not a new physical scalar species.

## 1. Derive the signed Jordan action, not a positive-gravity toy

At one Jordan epoch define Q=k_comoving²/a_J² and
s=delta sigma. Use the static Newtonian-gauge metric

  ds_J²=-(1+2 Psi_J)dt_c²+
            a_J²(1-2 Phi_J)delta_ij dx^i dx^j.

In the high-k weak-field frozen-coefficient limit of
F(sigma)R/2 + k1(sigma)Z -V(sigma), and after integration by
parts, the quadratic Jordan metric and scalar spatial action
per physical volume is

 L_J,2 = F Q(Phi_J² -2Phi_J Psi_J)
         +F_sigma Q s(Psi_J -2Phi_J)
         -(k1 Q+m_eff²)s²/2
         -delta_rho_b Psi_J.

m_eff² denotes the local effective scalar mass including
potential/curvature background contributions; it is **not**
identically V_sigma_sigma on a general rolling background.
The sign of the gravitational Hessian is indefinite because
Psi_J is a constraint, not a freely propagating positive
energy mode. This is why the earlier positive-definite toy
cannot be used as the physical gravitational coupling.

The structural spatial action adds the exact reduced Schur
complement term

 L_support,2 = -Kgeom(k,a) (s/sigma-Phi_J)²,

 Kgeom=C omega0²/(omega0²+Lambda),
 C=4W eta_K k⁴/a_s⁴,
 omega0²=m_L²+c_L² k²/a_J²,
 Lambda=2C/(9K), K=B f_L²/3.

It does not contribute an independent CDM density. The
remaining kinetic braiding of G3, background time derivatives,
ADM lapse/shift and tensor/trace structural sectors must be
derived separately for the accepted late300 solution.

## 2. Variation gives three signed Jordan constraints

Vary the combined restricted L2 with respect to
x=(Psi_J, Phi_J, s). Then H x=(delta_rho_b,0,0), with

  H = [
    [0, -2F Q, F_sigma Q],
    [-2F Q, 2F Q-2Kgeom, -2 F_sigma Q+2Kgeom/sigma],
    [F_sigma Q, -2 F_sigma Q+2Kgeom/sigma,
                      -(k1 Q+m_eff²)-2Kgeom/sigma²]
  ].

The explicit equations are

  -2F Q Phi_J+F_sigma Q s=delta_rho_b,                (C1)

  2F Q(Phi_J-Psi_J)-2 F_sigma Q s
       +2Kgeom(s/sigma-Phi_J)=0,                       (C2)

  F_sigma Q(Psi_J-2Phi_J)
       -(k1 Q+m_eff²)s
       -2Kgeom(s/sigma-Phi_J)/sigma=0.                (C3)

This matrix is **physically signed** for the restricted
covariant-action subclass, unlike the earlier positive
stiffness comparison. An arbitrary overall action-sign
change would also change the source term and is not made.

The supporting-space term changes both (C2) and (C3).
Because the curvature-squared spatial action is quadratic
around the unexcited trace branch, and carries no linear
Psi_J factor at this order, it does NOT itself introduce
a separate first-order source in (C1). This should not be
confused with saying that a quadratic structural action cannot
affect linear gravitational equations.

The modified slip condition in this restricted action is

  Phi_J-Psi_J =
      (F_sigma/F)s
      -[Kgeom/(F Q)](s/sigma-Phi_J).

In particular the measured Weyl potential can change even
though no CDM fluid was added. Whether its change can match
real lensing depends on the complete rolling G3 theory.

## 3. Known scalar-tensor benchmark, no structural kernel

For Kgeom=0, solve (C1-C3). Define
Psi_GR=-delta_rho_b/(2F_ref Q) for the same baryon
density perturbation, and

  mu=Psi_J/Psi_GR,
  Sigma=(Phi_J+Psi_J)/(2Psi_GR).

Then for the restricted subclass,

  mu_ST = F_ref/F [1+
    F_sigma² Q/{2F(k1 Q+m_eff²)+3F_sigma² Q}],

  Sigma_ST=F_ref/F,

  Phi_J-Psi_J=F_sigma s/F.

This is the standard scalar-tensor static linear benchmark.
For healthy k1>=0, mass2>=0, 1<=mu_ST/(Fref/F)<=4/3.
For k1=mass2=0 and F_sigma!=0, the familiar massless
Brans-Dicke/f(R) scalar force gives exactly 4/3.
For F_sigma=0 it reduces to GR with effective F.

This bound applies ONLY to the non-braided, positive-kinetic
subclass. It does NOT apply to the full late300 action with
rolling g(sigma)Z, nor to the physical averaged geometric source.

## 4. Actual accepted covariant EFT coefficients also implemented

The accepted action has, on a homogeneous rolling background,

  alpha_M = sigma_dot F_sigma/(H_J F),

  alpha_B = sigma_dot(2Z g-F_sigma)/(H_J F),

  alpha_K = [2Z/(H_J² F)]
       [k1+6k2 Z-4g_sigma Z]
       +12 sigma_dot Z g/(H_J F),

  D = alpha_K + (3/2) alpha_B².

These are the Bellini-Sawicki variables for this restricted
luminal Horndeski operator class, with G3=g(sigma)Z and
G4=F(sigma)/2. On a no-slip background,

  alpha_B=-2alpha_M
  => Z g = -F_sigma/2
     (provided sigma_dot is nonzero).

The module computes and tests these identities independently.
They determine the requisite *time-dependent coefficients*,
but do not, by themselves, supply the full ADM scalar
constraint Hessian of the rolling action. In the non-braided
G3=0 subcase, alpha_B=-alpha_M, NOT the accepted No-Slip
condition, except when alpha_M=0.


### Exact on-shell No-Slip covariant kinetic identity

For a rolling background with Z>0, the accepted No-Slip condition gives

  Zg=-F_sigma/2,

and therefore the G3 contribution to alpha_K becomes

  12 sigma_dot Z g/(H_J F)=-6 alpha_M.

Define the explicit microscopic kinetic combination

  K_cov = [2Z/(H_J² F)](k1+6k2 Z-4g_sigma Z).

Then

  alpha_K = K_cov - 6 alpha_M,

  D = alpha_K + 3 alpha_B²/2
    = K_cov - 6 alpha_M + 6 alpha_M²,

  Q_s=2F D/(2-alpha_B)².

This is an exact *necessary background identity* of the accepted rolling
G3 Horndeski subclass, not a derived No-CDM source. If D<=0 or Q_s<=0,
that branch is scalar-ghost unstable even if the No-Slip identity holds.
The sound-speed numerator and tensor/structural mixing require separate
gradient/constraint equations.

The function no_slip_kinetic_closure checks the identity and the
necessary kinetic-positivity condition over arbitrary supplied
background coefficient sets. It does not import an archived CDM density.

## 5. Numerical interpretation

The module tests the full signed H via numerical differentiation
of the **same** quadratic Jordan action, solves (C1-C3)
simultaneously, verifies the F-only analytic mu and Sigma,
the 4/3 massless scalar-tensor limit, and the scale-dependent
change after adding Kgeom.

Do not interpret the numerical sample, including its
matter source, as a late300 prediction. The F and scalar
coefficients are example values chosen for unit-regression
tests, and the structural source is not action-derived on
a cosmic inhomogeneous background.

A scan is useful because it demonstrates two distinct facts:
*the physical signed Jordan static constraints do not obey
the earlier positive-definite toy-screening theorem*, and
*the actual reduced structural kernel does not automatically
supply the order-five additional force or a cold-matter
perturbation history*.

## 6. Still missing for the accepted SDMC no-CDM action

The next genuinely physical step is to obtain the complete
ADM quadratic action for the accepted rolling G3, G2 and F
coefficients, including their time derivatives and
lapse/shift constraints; add the full supporting-space scalar,
trace and traceless operators with proper covariance and
structural clock mapping; and solve for Phi_J, Psi_J and
delta sigma together with the primordial source.

Then test physical mu(k,a), Sigma(k,a), kinetic positivity,
sound speeds and causal initial conditions in the full
Boltzmann hierarchy.

Current file: sdmc/pg_jordan_qs_covariant.py
Tests: tests/test_pg_jordan_qs_covariant.py
Workflow: .github/workflows/sdmc_pg_jordan_qs_covariant.yml

External cross-checks:
- Clifton et al., Living Reviews in Relativity (2012), scalar-tensor
  Jordan field equations and force limits.
- Baker et al., Testing General Relativity in Cosmology,
  Living Reviews in Relativity (2018), Jordan-frame
  Brans-Dicke perturbations and Weyl potential.
- Bellini & Sawicki (2014), EFT of Horndeski linear perturbations.
