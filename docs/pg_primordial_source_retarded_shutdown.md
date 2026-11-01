# SDMC Manuscript B Part III — Primordial-Well Source, Retarded Geometry, and the Selector-Switch-Off Theorem

## 484. The primordial well is a sourced state, not an assumed free oscillation

The continuous-fabric gravitational-well description in Manuscript B,
Part III Sections 101, 111 and 112 distinguishes:
 (i) the conserved original material inventory M_u and its redistribution,
 (ii) the supporting structural spatial metric gamma_ij,
 (iii) the Jordan potentials Phi_J and Psi_J,
 (iv) the field's source-supported static or slowly varying deformation,
 (v) any homogeneous/free relative-strain mode excited by a changing source.

A well does not need to oscillate to exist. A positive mass term allows
the geometry to relax toward its source-dependent equilibrium; oscillations
are transient homogeneous response about that equilibrium if the system
is underdamped.

At the exact Planck boundary S=1, rho_v,P=rho_P, rho_X,P=0 and
chi_P=0. The classical mode equation begins only after the boundary
on the one-sided classical branch S_A approximately 1. No equation
here is extrapolated to t<t_P, and no primordial mode amplitude
is guessed from rho_P.

## 485. Source exactly present in the accepted Section 106 action

Manuscript B Eq (106.3) supplies the classical relative-trace action
for one mode, with

  K_ell=B f_L²/3,
  omega0²=m_L²+c_L² k²/a_J²,
  C_k=4 W eta_K k^4/a_s^4,
  Lambda_k=2 C_k/(9 K_ell)
          =8W eta_K k^4/(3B f_L² a_s^4),
  u_k=m_k-Phi_J.

The action gives Eq (106.4),

  ddot ell+(3H_J+dot K_ell/K_ell)dot ell
      +(omega0²+Lambda_k)ell=-3 Lambda_k u_k.

The source is proportional to the curvature selector W.
Neither u_k nor Phi_J can be assigned simply from a material
mass fraction; they require the full scalar+metric ADM constraints.

For a fixed instant with slowly varying u_k and coefficients,
the exact algebraic minimizer is

  ell_eq=-3Lambda u/(omega0²+Lambda),

  q_eq=u+ell_eq/3 = omega0² u/(omega0²+Lambda).

The static reduced energy for one REAL cosine mode is

  E(ell)=K_ell omega0² ell²/4+C_k(u+ell/3)²/2.

At its stationary value,

  E_min=[C_k omega0²/(2(omega0²+Lambda))]u².

This is the same positive static Schur kernel derived earlier
in the SDMC supporting-space investigation. Crucially,
positive kernel energy does not imply extra attractive Jordan
gravity without the solved lapse, shift and metric constraints.

## 486. The curvature-source shutdown result

Consider a fixed finite u_k and a nonzero restoring frequency
omega0²>0. If W ->0 at fixed finite K_ell>0,

  C_k ->0, Lambda_k ->0,
  ell_eq ->0,
  E_min ->0,
  K_geom,eff=C_k omega0²/(omega0²+Lambda_k)->0.

For a fully switched-off W and no other direct coupling,
the ACCEPTED reduced action therefore has no stationary
source-supported trace well:

  ddot ell+(3H_J+dot K_ell/K_ell)dot ell
              +omega0² ell=0.

A preexisting excited trace can have nonzero dynamic energy
during subsequent expansion, but its equilibrium is L=0.
This is not yet persistence of a stationary *primordial well*.

A critical geometric distinction: W->0 removes the
curvature-SELECTOR ENERGY and the source of the RELATIVE
strain ell; it does NOT erase the ordinary curved Jordan metric.

Even with exactly L=0, the two metrics are locked as
  gamma_ij=M² h_ij.
In Newtonian gauge a spatial Jordan potential Phi_J and
local mapper fluctuation m_k give
  q_k=m_k-Phi_J+ell_k/3.
For a Fourier mode, Manuscript B Eq (102.5) gives
  delta R_s(k)=4 k² q_k/a_s².

Hence when W=0, ell_eq=0, but u=m-Phi_J !=0,
  delta R_s(k)=4 k² u/a_s² !=0,
  E_selector=0.

A primordial gravitational well carried by the ordinary
Jordan matter/metric sector can therefore persist under
exact relative locking. The above source-shutdown result
rules out a *separate static relative-strain well*
driven by THIS curvature operator after W=0;
it does not rule out ordinary Jordan gravitational
wells or source effects from the complete action.

The k=0 mode cannot be driven by this curvature operator
even when W>0, because C_k proportional to k^4.
This restricts claims about homogeneous "primordial wells"
in the reduced Section 106 operator.

There are three logically different ways to obtain a persistent
primordial well:
  1. W has a nonzero late value, constrained by infrared
     curvature/gradient signatures and the accepted handoff;
  2. a separate explicit physical source coupling survives when
     W=0, derived from the full matter+mapper action rather than
     added as a cold-density proxy;
  3. nonlinear supporting-space dynamics produce a stable
     nonzero configuration not captured by the positive
     quadratic trace operator.

The present accepted reduced source supports none of these
as a proven independent late geometric mass component.
In particular, the locked Jordan gravitational well can
exist without a nonzero independent L source.

## 487. Retarded memory and the adiabatic/nonadiabatic split

The general linear equation can be expressed with a causal
retarded Green function,

  ell(t)=ell_hom(t)+integral_{t_A}^{t} G_ret(t,t')
                                  [-3Lambda(t')u(t')]dt'.

Here ell_hom depends on initial structural data.
The retarded term depends on the full history of the
curvature-driven source.

Introduce the moving equilibrium

  ell(t)=ell_eq(t)+delta ell(t).

The exact deviation equation is

  ddot(delta ell)+Gamma dot(delta ell)+Omega² delta ell
    =-[ddot ell_eq+Gamma dot ell_eq],

  Gamma=3H_J+dot K_ell/K_ell,
  Omega²=omega0²+Lambda.

Slow equilibrium motion produces little freely oscillating
excitation when Omega is high and the initial field tracks
equilibrium. Rapid changes can excite a nonzero delta ell.
The amount retained depends on the whole source history
and the damping, and need not be monotonic in pulse width.

If the source vanishes and the residual equation is damped
with a positive restoring frequency, the free solution
decays in physical amplitude. In a stable constant-coefficient
setting its solutions have exponents

  r_± = -Gamma/2 ± sqrt(Gamma²/4-Omega²).

A stable free mode is not proof of permanent geometric memory.

## 488. Conserved material inventory and a separate metric source

Let the conserved material be partitioned as

  M_u=M_remaining+M_inferred,
  M_remaining=(1-f)M_u,
  M_inferred=f M_u.

For a specified NORMALIZED periodic, wrapped Gaussian
comoving material profile with two widths sigma_c and
sigma_d, the cosine density projection is

  rho_M,k =
    [2 M_u/(a_J³ L_box³)]*
     [(1-f)e^(-k² sigma_c²/2)
       +f e^(-k² sigma_d²/2)].

This calculation conserves material mass exactly while
the inhomogeneous density-mode amplitude changes.
No second material mass is created.

The geometric metric mismatch u_k=m_k-Phi_J is *not*
identical to this density coefficient. Only a physical
coupled Jordan constraint can determine its actual
dependence on mass redistribution and the scalar state.

To stress-test the retarded response, a dimensionless
mismatch u_k is externally PRESCRIBED proportional to the
normalized shape factor in the preceding bracket, with
an arbitrary small amplitude. This is a diagnostic forcing
history, not the physical primordial-to-Jordan coupling.

## 489. Numerical pulse and aftermath

A mathematically convenient *compact-support* W(x) is used
so W=0 exactly before and after a finite interval,
where x=ln a_J. Its envelope is

 W(x)=W_peak[1-((x-x_h)/Delta)²]²
                 for |x-x_h|<Delta,
 W=0 outside.

B(x) rises smoothly from a positive B_floor to 1
to keep the classical kinetic rank regular.

The mode is evolved using Jordan time, with
H_J=H0 exp(-h x), a_s=exp(p_v x), and
full Kdot/K friction. Both B(x) and W(x) histories
are illustrative post-boundary profiles and
are NOT derived from primordial SDMC dynamics.

The actual source in the solved differential equation
is *only* -3 Lambda(x) u(x), exactly the curvature-mismatch
term in Manuscript B. There is no phenomenological
CDM density or an artificially appended long-lived
matter-to-trace coupling in the equation.

The output explicitly distinguishes
  - material conservation and the changes of its
    spatial Fourier profile,
  - quasi-static source-supported structural deformation,
  - transient field mismatch after a changing source,
  - actual W=0 tail with vanishing stationary forcing,
  - positive transient free-mode energy after shutdown.

A second run changes the initial trace velocity while
keeping the same material/source history. Any difference
then belongs to a homogeneous structural mode, rather
than to the source-driven primordial well.

## 490. Field energy versus total interaction energy

A separate important point from Manuscript B Sec 111.2:
when a physical source term is held fixed, a deformed
field can have positive bare quadratic energy while
the complete source-field interaction contributes with
the opposite sign. The effective interaction energy must
come from the original complete source coupling.

For the Section 106 curvature-mismatch term,
the energy minimum E_min is positive as a function of
the externally prescribed metric mismatch u. The
observable sign and magnitude of the Jordan potential
change are NOT determined by this static minimum alone.

Neither the historical material mass M_u nor the
inferred partition M_inferred may be added again as a
separate new mass solely because this geometric
energy is nonzero.

## 491. Consequence for the complete SDMC no-CDM theory

The source-driven primordial-well route and the
freely oscillating geometric-field route are distinct.

The complete physical system must compute:
  - the initial structural metric and material configuration
    in a valid one-sided post-Planck initial state;
  - the conserved material current, not a duplicate mass;
  - the direct interaction of matter with L, if present,
    including its B/W and mapper dependence;
  - the coupled metric mismatch m-Phi_J, from
    F/G2/G3 + supporting-space lapse/shift constraints;
  - the cosmological evolution of actual rho_g,
    P_g and transfer J_g under source-driven expansion;
  - the observed Jordan Phi_J/Psi_J, growth, acoustic
    peaks and lensing after removing external CDM.

Most importantly, if the accepted thermal/flatness
handoff sets W=0 and the full physical parent action
contains no other matter-to-L source, the positive
quadratic trace cannot retain a *stationary source-held*
primordial depression after that handoff. A residual
free oscillation is a different possible physical sector,
whose abundance, clustering and lifetime require
separate physical determination.

Core code:
  sdmc/pg_primordial_well_retarded.py
  tests/test_pg_primordial_well_retarded.py
  .github/workflows/sdmc_pg_primordial_well_retarded.yml

Manuscript reference: SDMC Manuscript B, Part III
Sections 99, 100, 101, 102, 106, 111, 112 and 115.
