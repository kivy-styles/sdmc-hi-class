# SDMC Manuscript B, Part III — geometric-source linear perturbation gate

## Scope and provenance

The implementation preserves the difference between conserved material inventory
M_u, supporting space V_s, bare structural density rho_v, active scalar rho_X,
and geometric gravitational source rho_g. The earlier uniform nested-sphere
mass ratio is a special material-partition calculation, not a rho_g perturbation
equation. Structural R_s=c t_s is controlled by structural time; observable
perturbations are evaluated using the Jordan metric and time t_c.

This branch contains (a) a **reduced action-derived structural-trace
response** and (b) a **conditional effective-source conservation/growth
control**. The second must not be claimed as a prediction of the first.

## A. Background and activation

The candidate classical source budget is

  3 F H_J^2 = rho_b + rho_r + rho_X + rho_g,
  rho_v=rho_P/S^2, rho_X=chi*rho_v, rho_g=Gamma_g*rho_v.

If rho_g is separately conserved and pressureless in Jordan time,

  dot(rho_g)+3 H_J rho_g=0  =>  rho_g=rho_g0 a_J^(-3).

The structural identity implies the **necessary reconstructed history**

  d ln Gamma_g/d ln a_J = 2 p_v - 3,
  p_v=d ln S/d ln a_J.

gamma_for_dust implements this algebraic target. It does not evolve Gamma_g
from the structural action.

## B. Exact conditional dust perturbations in the observable Jordan frame

With conformal Newtonian gauge

  ds_J^2=a_J^2[-(1+2 Psi_J)d eta^2+(1-2 Phi_J)dx^i dx_i],

and Fourier velocity divergence theta_g, assuming that the **physical
action-derived** T_g is independently conserved, pressureless and shear-free:

  delta_g' = -theta_g + 3 Phi_J'
  theta_g' = -mathcalH theta_g + k^2 Psi_J,
  mathcalH = a_J H_J/c.

Primes here use conformal distance c*eta (set c=1 in the field equations).
The potentials must follow from the full Jordan Einstein/scalar constraints;
this separate conservation pair does not solve those constraints.

A source with the same background density, initial conditions, vanishing
pressure and shear, and metric coupling as a conventional cold fluid is
**mathematically indistinguishable from CDM at linear order**. That alone
cannot establish a new physical origin.

## C. What the existing quadratic supporting-space action really derives

Manuscript B Part II, Sections 65–67, has the restricted action

  S_g = integral dt_c d^3x N_ADM sqrt(h) [
       B f_L^2/2 * ((D_n L)^2-c_L^2(D L)^2-m_L^2 L^2)
       - W V_K[gamma] ],

  V_K = eta_K/2 * R_ijkl[gamma] R^ijkl[gamma],
  gamma_ij = M^2 h_ik (exp(2L))^k_j.

For scalar L^i_j=(ell/3)delta^i_j, m=delta ln M, define

  q=m+ell/3-Phi_J, a_s=S/S0,
  omega0^2=m_L^2+c_L^2 k^2/a_J^2,
  K_ell=B f_L^2/3,
  C=4 W eta_K k^4/a_s^4,
  Lambda=2 C/(9 K_ell)
        =8 W eta_K k^4/(3 B f_L^2 a_s^4).

For a Jordan-time Fourier mode with externally specified metric/mapper driver,
the reduced scalar-trace Euler-Lagrange equation is

  ddot(ell) + [3 H_J + dot(K_ell)/K_ell] dot(ell)
   + [omega0^2+Lambda] ell = 3 Lambda (Phi_J-m).

trace_dynamic_derivatives implements this dynamical equation, and a test
verifies that its stationary limit matches trace_static_response. The
driver Phi_J-m is NOT determined within this reduced sector.

At fixed coefficients, the reduced quadratic trace energy is

  E = K_ell/2 * [dot(ell)^2 + omega0^2 ell^2] + C q^2

and its first variation is

  delta E = K_ell (dot(ell_bar) delta dot(ell)
            +omega0^2 ell_bar delta ell) + 2 C q_bar delta q.

**If the structural branch is exactly locked and flat**, ell_bar=0,
dot(ell_bar)=0, q_bar=0, so delta E=0 at first order, even with
nonzero linear ell and q fluctuations. Coefficient variations multiply
vanishing background quadratic amplitudes. This is a concrete obstacle
to obtaining a nonzero first-order rho_g by linearizing this unexcited
quadratic action alone.

At fixed background drivers, the quasistatic relative-trace equation gives

  ell = 3 Lambda/(omega0^2+Lambda) (Phi_J-m)
  q   = -omega0^2/(omega0^2+Lambda) (Phi_J-m).

The functions trace_static_response and quadratic_energy_variation check
these equations. The lapse/shift, full F/G2/G3 scalar and structural
constraints have *not* been solved together, and the reduced q energy is
not an independently verified local Jordan stress tensor.

## D. Implementation: conditional subhorizon growth control

For k >> a_J H_J/c, with F constant, smooth X and no neutrino or photon
perturbations, define x=ln a_J, u=theta/(a_J H_J/c),
Omega_b(a), Omega_g(a), and d ln H/dx. Neglect Phi_J' only in this
subhorizon limit:

  d delta_b/dx = -u_b
  d u_b/dx = -(2+d ln H/dx)u_b
               -(3/2) mu [Omega_b delta_b + Omega_g delta_g]

  d delta_g/dx = -u_g
  d u_g/dx = -(2+d ln H/dx)u_g
               -(3/2) mu [Omega_b delta_b + Omega_g delta_g]
               +c_s,g^2 (k/mathcalH)^2 delta_g.

The last pressure term is an approximate **Newtonian subhorizon closure**
and not the exact relativistic nonadiabatic-pressure hierarchy. mu=1
is a GR control, *not* the accepted SDMC metric response.

Cases tested on the same background and simple illustrative post-decoupling
seeds: dust pg; smooth pg with no clustering; and finite effective sound
speed. The initial conditions are not recombination transfer functions,
the code does not contain a Boltzmann photon/neutrino hierarchy, and it
does not calculate CMB or observed cosmic shear.

The reference density assigned to Omega_g equals the archived late300
cold-like fitted density solely to make a quantitative comparator; this
does not show the structural action predicts the source amplitude.

## E. Testing and interpretation

Run with NumPy and SciPy:

    python3 -m pip install numpy scipy
    python3 -m unittest discover -s tests -p 'test_pg_linear_perturbations.py' -v
    python3 -m sdmc.pg_linear_perturbations

The regression tests cover: exactly zero first-order energy on the
unexcited locked branch, nonzero variation with a nontrivial background,
static k^4 curvature stiffness, exact Newtonian-gauge dust conservation,
Gamma_g reconstruction, special r_u/R mass ratio, EdS growth D+=a,
clustering-vs-smooth pg, finite-pressure suppression, and Poisson sign.

## F. What remains before a physical SDMC no-CDM implementation

Independently derive normalized bar(rho_g), delta(rho_g), delta(P_g),
momentum density, anisotropic stress, and metric constraints Phi_J/Psi_J
from the complete supporting-space action on a nontrivial background or
from a consistently averaged structural state. The constraint algebra,
gauge-invariant physical degrees of freedom, sound speeds and initial
conditions must be established. Only then may pg be inserted into a
complete hi_class/CLASS perturbation hierarchy, with separate baryon,
photon, neutrino and X evolution, to test Planck, BAO and lensing.

No conventional CDM species is inserted into this module; however,
assuming a pressureless pg with standard gravity is an **effective cold
fluid closure** and not an independent gravitational mechanism.

References: SDMC Manuscript B Part II Sections 63–67; Manuscript B Part III
Sections 98, 110, 114 and 116; Ma & Bertschinger (1995),
https://arxiv.org/abs/astro-ph/9506072.
