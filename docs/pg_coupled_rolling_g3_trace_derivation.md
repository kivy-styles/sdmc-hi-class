# SDMC Manuscript B, Part III — Rolling G3 plus supporting-space trace on the locked-flat branch

## 448. Joint action and the physical distinction

The accepted scalar–tensor sector is

  S_J=int sqrt(-g)[F(sigma)R/2+
     k1(sigma) Z+k2(sigma)Z²-V(sigma)-g(sigma)Z Box sigma].

The existing supporting-space restricted action is

  S_g=int N sqrt(h)[K/2 (D_n ell)^2
       -K c_L²/(2a_J²)|grad ell|²
       -K m_L² ell²/2
       -W eta_K/2 R_ijkl[gamma]R^ijkl[gamma]],

where K=B f_L²/3 and ell=tr(L), L^i_j=ell delta^i_j/3.

This branch derives the combined quadratic *principal scalar action*
about an exactly locked-flat unexcited supporting-space background
ell_bar=0, grad ell_bar=0, R_ijkl[gamma_bar]=0, with sigma
rolling on the classical scalar-tensor background. It does not
equate the local supporting curvature dilution to global FLRW.

This branch has **no homogeneous geometric rho_g**: the supporting
action starts at second order around this background. Consequently,
the calculated frequencies are structural scalar eigenmodes, not
cosmological transfer functions for a matterlike geometric source.

All observable times are Jordan proper time t_c. The structural ruler
remains a_s=S/S0 and R_s=c t_s.

## 449. Derive the relative-curvature combination exactly

In unitary gauge delta sigma=0 and

  h_ij=a_J² exp(2 zeta)delta_ij,
  gamma_ij=M² h_ik (exp2L)^k_j,
  M=a_s/a_J.

For trace L^i_j=ell delta^i_j/3,

  gamma_ij=a_s² exp[2(zeta+ell/3)]delta_ij.

Hence the scalar conformal structural perturbation is

  q=zeta+ell/3.

For the x-dependent real cosine Fourier mode, the quadratic
supporting-space curvature energy is

  Ecurv=C/2 (zeta+ell/3)²,
  C=4W eta_K k⁴/a_s⁴,

where the factor 1/2 is the spatial average cos²(kx)=1/2.

The positive structural trace Lagrangian kinetic term is
(K/4)dot ell², also due to cosine averaging.

Since all supporting-space operators are at least quadratic
in the *unexcited branch fluctuations*, their N and shift
variations enter at cubic order. The linear lapse/shift
constraints at quadratic order therefore coincide with those
of the rolling G2/G3/F action established in Section 443.
This decoupling is specific to the exactly unexcited background;
a physically inhomogeneous averaged state can mix all constraints.

## 450. Two-field quadratic ADM action

For an isolated scalar–tensor cosmological background satisfying
its equations, the single scalar Horndeski curvature quadratic
action after lapse/shift elimination has

  Qs=2 F D/(2-alpha_B)²,
  Fs=(1/a_J)d/dt_c[a_J F²/Theta]-F,
  Theta=F H_J(1-alpha_B/2).

The standard scalar squared speed is cs0²=Fs/Qs. The expression
for Fs is the scalar-only gradient coefficient; additional
matter and supporting-background fields require their complete
constraint/kinetic matrix.

Combining it with the existing structural trace gives, for
one real cosine Fourier mode,

  L2 =
    (Qs/2)dot zeta² + (K/4)dot ell²
    -(Fs/2)(k²/a_J²) zeta²
    -(K/4)(m_L²+c_L²k²/a_J²)ell²
    -(C/2)(zeta+ell/3)².

Write y=(zeta,ell). Then

  L2=1/2 dot y^T A dot y -1/2 y^T M y,

  A= diag(Qs,K/2),

  M=[
   [Fs k²/a_J²+C,             C/3],
   [C/3,           K omega0²/2+C/9]
  ],

  omega0²=m_L²+c_L²k²/a_J².

The two local frozen-coefficient normal-mode frequencies are
the eigenvalues of

  A^(-1/2) M A^(-1/2).

The complete coupled branch is free of local quadratic kinetic
ghosts only if Qs>0 and K>0. If Fs>=0, omega0²>=0, C>=0, then
the displayed spatial restoring matrix M is positive
semidefinite as a sum of positive quadratic forms. These are
**sufficient conditions for this frozen quadratic system only**,
not the full cosmological no-CDM stability requirements.

The structural higher-derivative curvature term can produce
rapid high-k modes; it must be constrained by the action's
effective-field-theory cutoff and the physical Jordan clock.
No mode in h/Mpc has been calibrated from this test.

## 451. Exact return to the earlier static structural kernel

Vary the joint action with respect to ell while holding zeta
fixed in static equilibrium. The trace solution is

  ell*/zeta = -(C/3)/(K omega0²/2+C/9)
             = -3 Lambda/(omega0²+Lambda),

  Lambda=2C/(9K).

The resulting single-variable zeta stiffness is

  M_eff = Fs k²/a_J²+C-
          (C/3)²/(K omega0²/2+C/9)

        = Fs k²/a_J²+
          C omega0²/(omega0²+Lambda)

        = Fs k²/a_J²+Kgeom(k).

This exactly recovers the earlier action-derived static
structural Schur complement, now from a *combined* rolling
G3 plus structural trace quadratic action.

The structural source supplies a scale-dependent metric
response, but NOT an action-derived nonzero cold-density
background. The strongly constrained zero-mode limit and
finite-frequency kernel require proper treatment of gauge
and propagating structural degrees of freedom.

## 452. Numerical tests of the coupled action

The implementation checks:
 - the exact quadratic action Hessian via finite differences;
 - the scalar-only rolling G3 kinetic identity Qs;
 - the structural trace K/2 positive kinetic block;
 - the relative curvature q=zeta+ell/3;
 - frozen local eigenfrequency positivity on a healthy sample;
 - the independent k=0 and W=0 limiting systems;
 - energy conservation of the frozen autonomous Hamiltonian;
 - an unstable-gradient example even when Qs>0;
 - a scalar ghost example that positive structural kinetic
   cannot repair;
 - exact static elimination against Kgeom(k).

Run:

  python3 -m unittest discover -s tests -p 'test_pg_coupled_rolling_g3_trace.py' -v
  python3 -m sdmc.pg_coupled_rolling_g3_trace

## 453. Conditions still required for no-CDM SDMC

This derivation is the **quadratic principal scalar sector on a
locked-flat, zero-energy structural background**. A nonzero
cosmological geometric source rho_g is expected to correspond to
an inhomogeneous structural state or a dynamical averaged
background, for which:
- spatial averaging must be consistent with full ADM constraints,
- fluctuation and averaging operations generally do not commute,
- scalar, trace, traceless and Jordan constraint mixing may change,
- physical rho_g, P_g, momentum and anisotropic stress must
  be derived from off-shell metric variation,
- a stable P_g~0, J_g~0 regime must arise dynamically,
- baryons, radiation and photons must be initialized without
  a hidden CDM species in a Boltzmann hierarchy,
- Planck CMB, BAO and lensing must test the resulting physical
  transfer functions.

The quadratic trace already exists in the supporting-space
action; introducing it as a propagating geometric metric field
does not by itself establish a new dustlike component or its
cosmological abundance.

Core files:
  sdmc/pg_coupled_rolling_g3_trace.py
  tests/test_pg_coupled_rolling_g3_trace.py
  .github/workflows/sdmc_pg_coupled_rolling_g3_trace.yml

Background source: Manuscript B Part II Section 65 and Part III
structural/Jordan mapping. Scalar-action cross-check:
Bellini & Sawicki (2014), https://arxiv.org/abs/1404.3713.
