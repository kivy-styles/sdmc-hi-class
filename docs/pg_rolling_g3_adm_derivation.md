# SDMC Manuscript B, Part III — Rolling G3 ADM and Jordan lapse–shift kinetic constraints

## 442. Complete covariant scalar–tensor operator to unitary-gauge ADM

The accepted classical SDMC scalar–tensor action (before adding the
supporting-space sector) contains

    S_J = ∫d^4x sqrt(-g_J) [
        F(sigma)R_J/2 + k1(sigma)Z+k2(sigma)Z²-V(sigma)
        -g(sigma)Z Box(sigma) ],

    Z=-g_J^ab ∂_a sigma ∂_b sigma/2.

This calculation begins *after* the exact Planck quantum boundary,
where rho_v(S=1)=rho_P and rho_X=0. The classical rolling action is not
extrapolated to t<t_P.

On a monotonic sigma(t_c) branch, choose unitary gauge delta sigma=0.
Let N be the Jordan ADM lapse, h_ij the Jordan spatial metric,
K_ij its extrinsic curvature, and v=dot sigma>0. Then

    Z=v²/(2N²).

Covariant integrations by parts, retaining lapse dependence, give the
exact scalar–tensor gravitational ADM operator

    S_J ≡ ∫ dt_c d³x N sqrt(h) [
      A2(N,t_c)+A3(N,t_c)K+
      (F/2)( R^(3)+KijK^ij-K² ) ]       (boundary terms omitted),

    A2 = k1 Z + [k2-(2/3)g_sigma]Z² -V,

    A3 = -(F_sigma v)/N+(g/3)(v/N)³.

This identity holds for the selected G2/G3/G4 family with G5=0
and monotonic sigma in unitary gauge. It is the full ADM rewriting
of that *covariant scalar–tensor sector*, not of the separate
supporting-space L_ij, curvature W, structural clock or M_u sector.

The term -F_sigma v K/N comes from integrating F R by parts,
and +(g/3)(v/N)³K from integrating -g Z Box(sigma).
The g_sigma shift of the Z² coefficient must not be omitted.
The full lapse dependence is essential for the braiding terms.

## 443. Explicit lapse derivatives

Write X=v²/2 and p=k2-(2/3)g_sigma. Evaluated at N=1,

    A2_N  = -2 k1 X -4pX²,
    A2_NN =  6 k1 X+20pX²,

    A3_N  = F_sigma v -g v³,
    A3_NN = -2 F_sigma v +4g v³.

The accepted Horndeski braiding function is immediately

    alpha_B=-A3_N/(F H_J)
           =v(2Xg-F_sigma)/(F H_J),

and

    alpha_M=v F_sigma/(F H_J).

The other EFT functions are

    alpha_K = (2X/(F H_J²))[
                  k1+6k2 X-4g_sigma X]
              +12v Xg/(F H_J),

    D=alpha_K+(3/2)alpha_B².

The identity alpha_B=-2 alpha_M is the accepted No-Slip branch
only when Xg=-F_sigma/2. It is not guaranteed by the presence
of geometric pg.

## 444. Full ADM L2 kinetic constraint matrix

Let the Jordan scalar ADM perturbations have one Fourier cosine mode:

    N=1+alpha cos(kx),
    h_ij=a_J² exp[2 zeta cos(kx)]delta_ij,
    N_i=∂_i beta, B=(k²/a_J²) beta_amplitude,
    dot(zeta)(x)=zdot cos(kx).

At zeta=0 for the **kinetic Hessian** of the gravitational action,
integrate over the cosine mode, using <cos²(kx)>=1/2.

The exact quadratic action density in the amplitude variables
(alpha,B,zdot) is

    L2/a_J³ =
         U alpha² + 3Theta alpha zdot
         +Theta alpha B -F B zdot -(3F/2)zdot²,

    Theta=F H_J +A3_N/2
         =F H_J(1-alpha_B/2),

    U = A2_N/2 +A2_NN/4+
        3H_J A3_NN/4 -3F H_J²/2.

Hence the action's symmetric amplitude Hessian is

    H_kin=[
      [2U,     Theta,    3Theta],
      [Theta,  0,        -F],
      [3Theta,-F,        -3F] ].

The shift perturbation is a true constraint with no B² kinetic
entry. Variation with respect to B gives

    alpha = (F/Theta) zdot
          = zdot/[H_J(1-alpha_B/2)].

The lapse constraint then fixes B and gives the
dynamical curvature kinetic coefficient

    Q_s=3F+2U F²/Theta²,

where the *spatially averaged* L2 takes the form
(1/2)Q_s zdot² for a cosine mode.

Substituting A2_N,N N and A3_N,N N gives exactly

    Q_s = 2 F D/(2-alpha_B)²,

    D=alpha_K+(3/2)alpha_B².

This is a fully action-derived kinetic check for the complete
rolling G2/G3/F operator family. Q_s>0 is necessary for no
scalar ghost, provided the constraint block is regular.
Theta=0 or alpha_B=2 is a degenerate-constraint surface and
cannot be replaced by a numerical pseudoinverse.

This is the scalar-only curvature kinetic coefficient, NOT the
physical eigenvalue of the fully coupled geometric pg + Horndeski
ADM fields: the supporting-space trace has its own kinetic,
constraints and time derivative couplings.

## 445. Independent exact ADM numerical variation

Rather than relying only on the known alpha-function formula,
sdmc/pg_rolling_g3_adm.py evaluates the **exact nonlinear ADM
density** for a periodic one-dimensional Fourier profile embedded
in three spatial dimensions:

    h_ij=a² exp(2zeta)delta_ij,
    N=1+alpha cos(kx), N_i=∂_i beta.

The exact conformal 3-curvature for the x-dependent metric is

    R^(3)=a^-2 exp(-2zeta)[-4zeta_xx-2(zeta_x)²].

The extrinsic-curvature matrix has one x eigenvalue and equal
transverse eigenvalues. Differentiating the complete nonlinear
action numerically with respect to alpha, B and zdot gives an
independent 3-by-3 Hessian. It is checked against the analytical
H_kin, not against the EFT expression.

Other unit tests recover:
- the minimally coupled canonical Q_s=k1 X/H_J²;
- the positive k-essence Q_s=X(k1+6k2 X)/H_J²;
- the accepted rolling No-Slip Q_s identity;
- finite-difference A2/A3 lapse derivatives;
- a negative Q_s example even if algebraic
  coefficient identities are satisfied;
- an explicit degenerate braiding constraint.

## 446. Isolated scalar gradient stability checkpoint

For the **isolated** scalar–tensor system without baryonic matter
perturbations or the supporting space, the Horndeski curvature
gradient coefficient is

    F_s=(1/a_J)d/dt_c[a_J F²/Theta]-F,

    Theta=F H_J(1-alpha_B/2).

The scalar propagation speed in that reduced system is
c_s²=F_s/Q_s (for Q_s>0). The expression can also be obtained
from the complete Horndeski generalized-Galileon quadratic
action. It requires the physical Hdot and alpha_Bdot of an
on-shell background.

On the No-Slip branch alpha_B=-2alpha_M, and if the condition
holds in a time interval,

    Theta=F H_J(1+alpha_M),

    F_s=F/(1+alpha_M)[
        -Hdot_J/H_J²
        -alpha_Mdot/(H_J(1+alpha_M))].

The accepted No-Slip kinetic identity gives

    Q_s = F D/[2(1+alpha_M)²].

For the branch 1+alpha_M>0, a necessary
isolated-scalar gradient condition F_s>0 reads

    -Hdot_J/H_J² >
       alpha_Mdot/[H_J(1+alpha_M)].

This is a useful **conditional** stability check involving
the structural scalar's changing gravitational coupling. It
is not the full geometric no-CDM sound speed, nor a positive
condition automatically guaranteed by alpha_B=-2alpha_M.

The code verifies the isolated canonical-scalar cs²=1
when Hdot=-X/F on its background and verifies the No-Slip
gradient identity against an independently differentiated
Theta expression.

## 447. Exact boundary of the result

This work derives the full rolling G2/G3/F scalar–tensor
**ADM representation and kinetic constraint block**.
It is a substantial improvement on the earlier static,
non-braided subclass.

It does NOT yet derive the complete combined ADM
action for F/G2/G3 + L_ij + structural curvature W
with the geometry clock and physical M_u/V_s state.
The following are still needed for the no-CDM prediction:
1. add the *complete* supporting-space ADM lapse and shift
   dependence to the scalar–tensor action, consistently
   distinguishing Jordan and structural clocks;
2. expand the total system to second order including ell,
   traceless Q_ij, zeta, lapse and shift, scalar-matter
   perturbations and the physical inhomogeneous background;
3. solve the coupled constraints, evaluate all scalar and
   tensor eigenvalues and propagation speeds;
4. derive the physical geometric background rho_g and
   its Jordan density, momentum, pressure and shear sources;
5. calculate the physical Newtonian and Weyl potentials
   without inserting a calibrated CDM-like source;
6. run the full CMB, BAO, cosmic-shear and clustering likelihoods.

Files:
    sdmc/pg_rolling_g3_adm.py
    tests/test_pg_rolling_g3_adm.py
    .github/workflows/sdmc_pg_rolling_g3_adm.yml

Source orientation: SDMC Manuscript B Part II Section 65 and
the previously established restricted Horndeski action.
Related perturbative identities: Bellini & Sawicki (2014),
https://arxiv.org/abs/1404.3713.
