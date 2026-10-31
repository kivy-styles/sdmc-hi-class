# SDMC Manuscript B, Part III — From structural wells to an averaged geometric source

## Scope and exact separation of geometries

The structural radius R_s=c t_s is defined using structural time, while
cosmological observables use Jordan proper time t_c and Jordan scale a_J.
The supporting spatial metric is

    gamma_ij = M^2 h_ik (exp(2 L))^k_j ,
    M=(S/S0)/a_J, a_s=S/S0.

V_s is structural space on which the material configuration V_u rests.
There is **no universal V_s - V_u shell assumption**. The special uniform
nested-sphere r_u/R material ratio does not by itself determine geometric
stress energy.

We investigate a restricted **one-dimensional periodic x-dependent trace**
embedded in a 3D spatial metric; this avoids assuming that local curvature
dilution equals homogeneous Jordan FLRW source dilution. A 1D periodic
profile is *not* an isotropic physical primordial spectrum.

## 1. Actual quadratic ADM supporting-space action

In the restricted trace branch of Manuscript B's candidate action,

    S_g = int dt d^3x N_ADM sqrt(h) [
      B f_L^2/2 ( (D_n L)^2 - c_L^2 (D L)^2 - m_L^2 L^2 )
      - W eta_K/2 R_ijkl[gamma] R^ijkl[gamma] ],

with L^i_j = ell delta^i_j/3, q = m + ell/3 - Phi_J,
m=delta ln M, and h_ij=a_J^2(1-2 Phi_J) delta_ij.

For a scalar conformal structural curvature perturbation,

    delta R_s = -4/a_s^2 partial_x^2 q
    R_ijkl R^ijkl = 8/a_s^4 (partial_x^2 q)^2

to quadratic order for this 1D profile.

Assuming **fixed coefficients and no additional explicit ADM-lapse
dependence**, variation of this restricted action gives the local quadratic
Jordan-volume energy proxy

    rho_g^(2)(x) = K/2 [
       dot(ell)^2 + m_L^2 ell^2
       + c_L^2/a_J^2 (partial_x ell)^2 ]
       + C (partial_x^2 q)^2,

    K=B f_L^2/3, C=4 W eta_K/a_s^4.

The code implements **this fixed-coefficient reduced expression**.
It does not establish the full mixed structural + Horndeski ADM
Hamiltonian, momentum or anisotropic-stress sources, and this energy
must not be used as if it were an independently conserved physical fluid.

## 2. Why zero-mean structural wells can give nonzero averaged energy

An inhomogeneous profile ell_s(x) may obey <ell_s>=0 while

    <rho_g^(2)> = (K/2) [
       <dot ell_s^2>+m_L^2<ell_s^2>
       +c_L^2/a_J^2<ell_s,x^2> ]
       +C<q_s,xx^2> > 0.

This is second order in the microscopic structural amplitudes. It
illustrates how a nontrivial physical structural background could be
represented by a nonzero **coarse-grained** energy. Such averaging does
not automatically identify that energy with the normalized cosmological
rho_g or establish its pressure or its gravitational coupling.

The EXACTLY locked flat, unexcited state (ell=dot ell=q=0) has
no first-order source from the quadratic terms:

    delta rho_g^(2) = 0

at zero background fields.

## 3. Derivation of a long-wavelength effective source

At a given epoch, keep K,C,a_J,a_s fixed. The first energy variation is

    delta rho_g =
      K [ dot(ell_bar) delta dot(ell)
          +m_L^2 ell_bar delta ell
          +(c_L^2/a_J^2)(partial_x ell_bar)(partial_x delta ell)]
      +2 C (partial_x^2 q_bar)(partial_x^2 delta q).

For a pre-existing short-scale structural field with slowly modulated
amplitude,

    ell(x) = [1+epsilon cos(k_L x)] A cos(k_s x),
    k_L << k_s,

the quadratic energy after averaging over scales ~1/k_s contains a
long-wavelength contrast proportional to

    delta_g,coarse ~ 2 epsilon cos(k_L x)

up to terms O[(k_L/k_s)^2] and O(epsilon^2).

The coefficient is calculable; it is not an assumed cold-density
perturbation. In the tests, k_s=32, k_L=2, epsilon=0.1, and a Gaussian
coarse-graining width sigma=0.16 are used in a 2*pi periodic domain.
The leading amplitude after Gaussian filtering is

    2 epsilon exp[-(k_L sigma)^2/2].

The numerical code measures this amplitude and checks the leading
separation approximation. Without the modulation, the filtered long
mode vanishes.

This is **mode coupling**: a large-scale variation arises from products
of small-scale inhomogeneous structural fields. It is first order
in the envelope modulation but second order in the microscopic field.

A real 3D cosmological model must define a statistically homogeneous,
approximately isotropic ensemble, its power spectrum and evolution.
No primordial ensemble normalization is supplied by this demonstration.

## 4. Fourier convolution underpinning the long-wavelength energy

Write ell(x)=sum_p ell_p exp(ip x),
q(x)=sum_p q_p exp(ip x).

The Fourier energy coefficient at long wave K is

    rho_g,K =
     K_l/2 sum_p [
       dot(ell)_p dot(ell)_(K-p)
       + (m_L^2 - c_L^2 p (K-p)/a_J^2)
          ell_p ell_(K-p) ]
       + C sum_p p^2 (K-p)^2 q_p q_(K-p).

Every source coefficient is a convolution, not an arbitrary cold fluid
amplitude. Its k dependence can differ strongly from dust, especially
for curvature-dominated configurations with four derivatives.

## 5. On-shell energy evolution: what the reduced action actually says

For a scalar Fourier trace mode with prescribed
u=(delta ln M)-Phi_J, define

    omega0^2=m_L^2+c_L^2 k^2/a_J^2,
    Ck=4W eta_K k^4/a_s^4,
    E=K/2 [dot(ell)^2+omega0^2 ell^2]
      +Ck(ell/3+u)^2.

The action-derived reduced trace equation for constant K is

    ddot(ell)+3 H_J dot(ell)+
      [omega0^2+2 Ck/(9K)]ell = -2 Ck/(3K)u.

Its **on-shell energy identity** (fixed K, c_L, m_L, W eta_K) is

    dot(E) = -3 H_J K dot(ell)^2
             -H_J K (c_L^2 k^2/a_J^2)ell^2
             -4 H_R,c Ck(ell/3+u)^2
             +2 Ck(ell/3+u)dot(u),

where H_R,c=d ln a_s/dt_c = p_v H_J. The last term is work
from changes in the external metric/mapper driver.

Thus for an unforced configuration (u=dot u=0), in an expanding Jordan
and structural epoch with nonnegative quadratic coefficients,

    dot(E) <= 0.

This is an energy-density statement. It does NOT imply the global
Jordan energy E_g,domain=a_J^3 <rho_g> stays fixed. It does NOT determine
P_g or J_g from a completed theory.

If fixed comoving field profiles and constant coefficients are imposed
rather than dynamically solved, the individual mean terms scale as:
    mass ~ a_J^0,
    spatial-gradient ~ a_J^-2,
    curvature ~ a_s^-4.

In the matter tracker a_s~a_J^(3/2), curvature ~a_J^-6.
These are local quadratic-energy proxy scalings, not the asserted
global FLRW gravitational source.

## 5a. Restricted Jordan isotropic pressure is a separate calculation

An energy density is not sufficient to claim pressureless gravitational
behavior. In the same reduced action, consider an **infinitesimal,
homogeneous Jordan scale variation** at fixed structural scale S (thus
fixed a_s=S/S0), fixed coefficients, and fixed fields ell, dot ell,
q and their comoving spatial profiles. Let

    Ekin=K/2 <dot ell^2>,
    Emass=K m_L^2/2 <ell^2>,
    Egrad=K c_L^2/(2 a_J^2) <(partial_x ell)^2>,
    Ecurv=C <(partial_x^2 q)^2>.

The reduced Lagrangian per Jordan physical volume is

    L_red = Ekin - Emass - Egrad - Ecurv,

and the corresponding restricted isotropic scale stress follows from

    P_red = (1/(3 a_J^2)) partial_aJ[a_J^3 L_red]|_(S,fields fixed)
          = Ekin - Emass - Egrad/3 - Ecurv.

It is an **isotropic scale-response diagnostic**. The full physical
pressure must still account for the off-shell coupling of the
structural fields, the full ADM constraints, local metric perturbations
and their averaging. For this explicit diagnostic, a nontrivial
STRICTLY STATIC structural field has Ekin=0 and hence P_red<0 on
the positive quadratic branch. It is not pressureless merely because
its local deformation persists.

A zero instantaneous P_red would require

    Ekin = Emass + Egrad/3 + Ecurv.

This is a necessary instantaneous balance in this restricted stress
calculation, not a sustained attractor or proof of dust dynamics.
The regression suite differentiates the reduced spatial action
numerically with respect to a_J, holding S fixed, and verifies the
pressure formula independently.

## 5b. Integrating out the static trace gives a spatial metric kernel

At a single epoch, for each Fourier mode, write u=m-Phi_J,
omega0^2=m_L^2+c_L^2 k^2/a_J^2,
K=Bf_L^2/3 and Ck=4W eta_K k^4/a_s^4.
The positive static energy is

    E_static=K omega0^2 ell^2/2+Ck (u+ell/3)^2.

The stationary trace ell_star= -3 Lambda/(omega0^2+Lambda) u
with Lambda=2Ck/(9K) gives the exact reduced functional

    E_eff(u)=Keff(k) u^2,
    Keff(k)=Ck omega0^2/(omega0^2+Lambda).

Consequently the fixed-driver response is

    d E_eff/d Phi_J = -2 Keff(k) u,

which has a scale-dependent k^4-derived kernel. This is a spatial
constraint-response term only; its variation must be added to those
of the **full** Jordan metric, structural scalar and matter action.
It is NOT a Poisson equation, a calibrated mu(k,a), a new CDM density,
or a complete cosmological anisotropic stress tensor.

Regression tests independently vary the static energy in ell and the
integrated energy in Phi_J and compare the results.

## 6. What is actually needed for pressureless geometric gravity

The desired Jordan global relation remains

    dot(rho_g)+3 H_J (rho_g+P_g)=J_g.

With J_g=0 and P_g=0, rho_g~a_J^-3.
For rho_g=Gamma_g rho_P/S^2 the necessary relation is

    d ln Gamma_g / d ln a_J = 2 p_v - 3.

These are **tests**, not solved consequences of the quadratic action.
The mode-coupling result supplies a legitimate nonzero averaged energy
from inhomogeneous wells, but the physical normalization, distribution,
pressure, momentum, Weyl potentials and evolution require the full
ADM+Horndeski constraints and primordial structural initial state.

Do NOT confuse:
  (1) rho_g quadratic proxy,
  (2) a statistical homogeneous average of it,
  (3) the physical rho_g entering Jordan FLRW,
  (4) the independent baryonic M_u inventory,
  (5) a separately assumed cold effective fluid.

## 7. Reproduction

NumPy and SciPy required for both the previous linear module and this
quadratic source:

    python3 -m unittest discover -s tests -p 'test_pg_quadratic_source.py' -v
    python3 -m unittest discover -s tests -p 'test_pg_linear_perturbations.py' -v
    python3 -m sdmc.pg_quadratic_source

GitHub Actions workflow:
  .github/workflows/sdmc_pg_quadratic_averaging.yml

References:
- SDMC Manuscript B Part II, Sections 63-67.
- SDMC Manuscript B Part III, Sections 98, 110, 114, 116.
