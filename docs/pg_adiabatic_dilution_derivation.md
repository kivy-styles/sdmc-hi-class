# SDMC Manuscript B — Part III: Structural-well energy dilution and the pressureless test

## 451. Physical separation of the bare radius, active density and geometry energy

The SDMC structural ruler is S and a_s=S/S0, with
R_s=c t_s and bare structural density rho_v=rho_P/S².
The observable Jordan scale is a_J with H_J=d ln a_J/dt_c.

The structural tracker in an epoch is p_v=d ln S/d ln a_J.
Thus rho_v scales a_J^(-2p_v). In radiation with p_v=2,
bare rho_v~a_J^-4; its gravitationally active
rho_X=chi rho_v can scale a_J^-3 only if chi~a_J.
That background identity does NOT determine structural-well energy
or gravitationally observable cold-like perturbations.

The supporting space V_s and the matter configuration V_u
are not identified as one quantity. M_u is the material
inventory, M_inferred is an observable inference, and
neither is automatically the energy in a scalar trace mode.

The current calculation follows one genuine existing reduced
supporting-space trace (not a new CDM species) to test whether its
*local quadratic energy* is pressureless and conserved. The source
is NOT the physical global FLRW geometric density rho_g.

## 452. Reduced trace oscillator and exact Jordan energy identity

For the real Fourier cosine mode ell(x,t_c)=A(t_c)cos(kx),
with a stationary zero external metric driver u=0,
the existing supporting-space quadratic action averaged
over one comoving cell is

 L_red = a_J³ K/4 [ dot A² - Omega²(t_c) A² ],

 K=B f_L²/3,
 Omega² = m_L²+c_L² k²/a_J²+Lambda,
 Lambda=2C/(9K), C=4W eta_K k⁴/a_s⁴.

All time derivatives are with respect to Jordan proper time t_c.
The trace equation is

 ddot A+(3H_J+dot K/K)dot A+Omega² A=0.

The four reduced energy components are

 Ekin = K/4 dot A²,
 Emass= K/4 m_L² A²,
 Egrad= K/4 (c_L²k²/a_J²) A²,
 Ecurv= C A²/18 =K/4 Lambda A².

 E=Ekin+Emass+Egrad+Ecurv.

Under the restricted homogeneous Jordan scale variation at fixed S,
the isotropic scale stress is

 P_red=Ekin-Emass-Egrad/3-Ecurv.

For arbitrary time-dependent K,c_L²,m_L²,W eta_K,
the general exact reduced energy exchange follows from the
on-shell Euler–Lagrange equation (cf. Section 416):

 J_red = dot E+3H_J(E+P_red)
       =dot K/4 [(m_L²+c_L²k²/a_J²)A²-dot A²]
        +K/4 [dot(m_L²)+dot(c_L²)k²/a_J²]A²
        +dot C/18 A².

For constant m_L²,c_L², parameterize local coefficient histories:

 x=ln a_J,
 H_J=H0 exp(-h x),
 a_s/a_s0=exp(p_v x),
 K=K0 exp(k_run x),
 (W eta_K)=(W eta_K)0 exp(w_run x).

Then

 Lambda=Lambda0 exp[(w_run-4p_v-k_run)x],

 A_xx+(3-h+k_run)A_x+Omega²/H_J² A=0,

 and the EXACT reduced Jordan exchange relation becomes

 J_red/H_J = k_run(Emass+Egrad-Ekin)
             +(w_run-4p_v)Ecurv.                   (452.1)

This last equation is not the actual total covariant J_g of
SDMC; its complementary exchange with F/G2/G3, metric and
material sources is not available in the reduced calculation.

The code independently derives dot E/H_J by a chain rule
using the above Euler–Lagrange equation and verifies the
identity against the direct reduced stress/exchange formula.

## 453. Adiabatic invariant and stress of the existing trace mode

Suppose Omega/H_J >>1 and
abs(dot Omega/Omega²)<<1, with background and couplings
varying slowly over one oscillation. Two independent
quadrature phases of the same structural mode provide
a phase-averaged calculation without adding another
physical particle species.

The adiabatic action invariant for the oscillator is

 I = a_J³ <E>/Omega = constant,

 so

 <E(x)>/<E(0)> ≈ exp(-3x)*Omega(x)/Omega(0).

At leading WKB order, averaged Ekin≈Emass+Egrad+Ecurv.
Consequently the restricted averaged stress obeys

 <P_red>/<E> ≈ G/(3 Omega²),
 G=c_L² k²/a_J².

The energy-exchange term reduces to

 <J_red>/(H_J<E>)
    ≈ [(w_run-4p_v-k_run)*Lambda]/(2 Omega²).

These identities obey

 d ln<E>/d ln a_J =
     -3[1+<P_red>/<E>] +<J_red>/(H_J<E>),

 with the analytical slope

 d ln<E>/d ln a_J =
    -3 +
    [-2G+(w_run-4p_v-k_run)Lambda]/(2Omega²).

Here a vanishing *reduced* pressure does not guarantee
separate conservation: curvature may still exchange
energy with the changing structural ruler.

## 454. Three clear dominant restoring regimes

(A) Mass-dominated, m_L²>>G,Lambda, constant m_L:
 Omega≈m_L, <P>/E≈0, Jred≈0,
 <E>~a_J^-3.

This is the familiar adiabatic behavior of a massive
oscillating field. The supporting-space trace belongs to
the existing SDMC gravitational action; nevertheless,
an independent oscillating trace population is NOT
thereby demonstrated to generate the observed extra gravity.
Its abundance, normalization and relation to physical
metric sources remain undetermined.

(B) Jordan gradient-dominated, G>>m_L²,Lambda:
 Omega~a_J^-1, <P>/E≈1/3, Jred≈0,
 <E>~a_J^-4.

(C) Structural curvature-dominated, Lambda>>m_L²,G:
 Omega~sqrt((W eta_K)/K)*a_s^-2.

For constant positive K and W eta_K:
 <P_red>/E≈0 but
 <J_red>/(H_J<E>)≈-2p_v,

 <E>~a_J^(-3)a_s^-2
     = a_J^(-3-2p_v).

In the matter tracker p_v=3/2:
 <E>~a_J^-6, Jred/(H_J E)≈-3.

In the radiation tracker p_v=2:
 <E>~a_J^-7, Jred/(H_J E)≈-4.

Those are results of local, reduced geometric curvature modes,
NOT a claimed global FLRW density law. They also differ
from the bare rho_v~a_J^(-2p_v).

The WKB assumption becomes less secure as curvature frequency
Omega~a_s^-2 drops toward H_J. The scaling cannot be extrapolated
to arbitrarily late epochs or to the exact Planck boundary.

## 455. Condition required for curvature-sector dust behavior

With running positive K and W eta_K, and curvature dominance,

 Omega~sqrt((W eta_K)/K)*a_s^-2.

Adiabatic dilution <E>~a_J^-3 is possible only if
the frequency remains constant, which requires

 d ln[(W eta_K)/K]/d ln a_J = 4p_v,

 or

 w_run = 4p_v+k_run.

For constant K this would require
W eta_K ~ a_s^4.

This is a **necessary conditional relation**, not a derived
coefficient history of the accepted SDMC structural action.
If a physical UV or late-flat activation makes the curvature
coupling vary differently, the actual source evolution must be
calculated, not tuned to follow an empirical CDM history.

In the compensated conditional case, both the WKB pressure
and the reduced Jred can average to approximately zero,
but the resulting mode would still need to source observable
metric perturbations correctly to replace CDM.

## 456. Numerical experiment and limits

The code uses x=ln a_J and directly solves the exact reduced
trace EL equation for two orthogonal initial phases, each
with equal initial energy. Ensemble averaging removes fast
oscillatory pressure variations to leading WKB order.

Implemented tests check:
 - arbitrary on-shell exact Jordan Ward residual;
 - independent central finite difference of the energy derivative;
 - analytical WKB pressure, exchange and dilution identity;
 - mass-dominated a_J^-3 dilution;
 - gradient-dominated a_J^-4 dilution;
 - curvature-dominated matter a_J^-6 dilution;
 - curvature-dominated radiation a_J^-7 dilution;
 - conditional W eta_K ~ a_s^4 cancellation;
 - consistency when K itself runs;
 - no implication P_red=0 => Jred=0;
 - separation of bare rho_v and reduced trace curvature energy;
 - WKB validity limitations at small Omega/H.

The numerical examples use arbitrary positive quadratic
coefficients and no external dark-matter normalization.

Core files:
  sdmc/pg_adiabatic_dilution.py
  tests/test_pg_adiabatic_dilution.py

## 457. Consequence for the no-CDM objective

The present result does not settle whether the completed
SDMC physical geometric density rho_g scales a_J^-3.
It does establish that a free, curvature-dominated structural
trace with constant couplings, even while behaving with
near-zero averaged *restricted pressure*, does not dilute
as a_J^-3: it exchanges energy through structural expansion.

Mass-dominated oscillation gives dustlike dilution for
the existing trace, but that is not evidence of an
independent gravitational source of the required strength,
clustering or lensing. An allowed compensated curvature
branch likewise requires a dynamically derived coefficient
history and the complete action constraints.

The decisive remaining calculation is to construct an
inhomogeneous supporting-space *background* and derive
its averaged physical Jordan metric stress tensor from
the combined F/G2/G3 + supporting-space ADM action,
without using a CDM-normalized source in the derivation.
