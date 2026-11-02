# SDMC Manuscript B — Part III: Conserved Material Redistribution, Correlated Wells, and Geometric Activation

## 478. The material inventory and nonuniform redistribution

The material partition of Manuscript B Sections 83 and 100 is

  M_u=M_remaining+M_inferred,  d M_u/dt_c=0.

M_inferred names material distributed into the expanded region,
not a second gravitational mass or a new particle species.
For nonuniform material positions and physical density profiles,
a finite comoving reference domain may be represented by

  rho_m(x,t)=a_J^-3 [M_remaining W_0(x,t)
                       +sum_i m_i W_i(x-x_i,t)],

where each comoving window integrates to one. Then

  M_inferred=sum_i m_i,   M_u=M_remaining+sum_i m_i.

The comoving mass-current conservation law is

  dot rho_m+3H_J rho_m+(1/a_J) div(rho_m v)=0.

A transfer rate Q_M through the boundary of the retained region
can be defined by relative material flux. Nothing in these equations
provides a unique spatial probability law or primordial correlation
spectrum. In particular the equation M_u=N_P m_P is a *conditional
Planck-cell normalization*, not proof of N_P distinct stable objects
or of an equal number of structural gravitational wells.

## 479. Fourier mass moments and genuine position correlations

For a finite periodic domain of comoving length L and nonzero
integer spatial mode k, define the material density contrast mode

  delta_k=[M_remaining W_core(k) exp(-ik_phys x_core)
                  +sum_i m_i W_i(k) exp(-ik_phys x_i)]/M_u,
  k_phys=2pi k/L.

This Fourier definition is for nonzero k; its zero mode
equals the conserved mean mass and is not a cosmological FLRW
spatial-curvature parameter.

For independent uniform positions, a normalized set of masses
w_i=m_i/M_u (when core mass is absent) gives

  E delta_k=0,
  E|delta_k|²=sum_i w_i²=1/N_eff,
  N_eff=(sum_i m_i)²/(sum_i m_i²).

Equal-mass units have N_eff=N. Unequal masses have N_eff<=N.

For independent positions drawn from the nonuniform one-point
probability law p(x)=[1+epsilon cos(k_phys x)]/L, the mean
Fourier amplitude is epsilon/2, and

  E|delta_k|²=epsilon²/4+(1-epsilon²/4)/N_eff.

This one-point modulation is NOT genuine pairwise clustering:
positions remain independent, but the spatial mean density
itself contains a coherent inhomogeneity.

A second construction has *true two-point correlations*.
Take n_groups independently uniform centers, each containing
n_c equal masses with independent Gaussian displacements
of variance sigma². Marginal individual positions remain
uniform, but members of the same group are correlated.
For N=n_groups n_c,

  E|delta_k|² =
   [1+(n_c-1)exp(-k_phys² sigma²)]/N.

At wavelengths much larger than the group spread, the power
approaches 1/n_groups rather than 1/N. At wavelengths much
shorter than the group spread, it approaches 1/N.

Example: N=4096, n_groups=256, n_c=16, k sigma=.1:
  independent power=1/4096=0.000244140625;
  group-correlated power=0.003869...,
  power enhancement approx 15.85.

Conservation of M_u is exact in both constructions.
Neither fragmentation nor conserved mass by itself fixes
the two-point statistics of the primordial configuration.

## 480. Action-level response of structural curvature to matter inhomogeneity

Manuscript B Section 106 defines the reduced scalar trace
L^i_j=(ell/3)delta^i_j and the relative curvature combination

  u_k=m_k-Phi_J(k),  q_k=u_k+ell_k/3,

where m_k is a local perturbation of the structural–Jordan
mapper and Phi_J is a Jordan metric potential.

For B>0, with a positive restoring branch, the quadratic
static potential energy of a real Fourier mode is

  E(ell_k,u_k)=K_l omega0² |ell_k|²/2
                         +C_k |u_k+ell_k/3|²,

  K_l=B f_L²/3,
  omega0²=m_L²+c_L²k²/a_J²,
  C_k=4W eta_K k⁴/a_s⁴,
  Lambda_k=2C_k/(9K_l)
          =8W eta_K k⁴/(3B f_L² a_s⁴).

For fixed u and positive coefficients the energy is
strictly minimized at

  ell_k=-3[Lambda_k/(omega0²+Lambda_k)]u_k
       =-3(1-T_k)u_k,

  q_k=T_k u_k,
  T_k=omega0²/(omega0²+Lambda_k).

The total on-shell static potential energy is

  E_min(k)=C_k T_k |u_k|²,

with separate contributions

  E_curv(k)=C_k T_k²|u_k|²,
  E_relative(k)=C_k T_k(1-T_k)|u_k|².

Their sum exactly equals E_min. This matters because
the structural-curvature energy alone is NOT the total
energy of the relaxed relative geometry: the relative
field's restoring energy must also be counted.

The off-shell variation of E_min with respect to u
is nonzero during activation, but its translation into
a physical additional Jordan gravitational force requires
the complete Jordan metric, matter, and mapper constraints.

At W=0 and finite B>0,
  C_k=Lambda_k=0,
  ell_k(forced,static)=0,
  E_min(forced,static)=0.
The structural spatial curvature can still inherit
Jordan spatial curvature (q=u). The free homogeneous
solution ell_h(t) need not vanish, but its initial spectrum
and survival must be independently determined.

Thus the temporary selector alone does not produce
a permanent, forced additional geometric source.

## 481. Conditional Poisson substitution and a k-four cancellation

Only for subhorizon, weak-field, nonrelativistic Jordan matter
and ignoring the dynamical mapper perturbation, assume

  -k² Phi_J(k)=A_m delta_m(k),
  A_m=4pi G_eff a_J² bar rho_m.

For u=-Phi_J this gives u=A_m delta_m/k².
At this level the on-shell stationary structural energy is

  E_min(k)=4W eta_K/a_s⁴ *
                    A_m² T_k |delta_m(k)|².

The apparent curvature k⁴ is canceled by the squared
Newtonian potential k^-4; the scale dependence remains
in the transfer T_k and the matter spectrum.

This is a useful cancellation but not a demonstrated
matter-to-structure source in the complete early SDMC
cosmology. At the exact Planck boundary, superhorizon
modes and the relativistic photon-baryon plasma,
the above Newtonian Poisson approximation is not valid.

## 482. Can this mode spectrum supply the required Gamma_g?

Section 108 defines

  rho_g=Gamma_g rho_v, rho_v=rho_P/S²,

and states that a pressureless Jordan effective source
requires

  d ln Gamma_g/d ln a_J=2p_v-3,
  p_v=dln S/dln a_J.

Let Gamma_E=E_min/rho_v be a *diagnostic only*, NOT a
replacement for the physical action-derived rho_g.
For conserved nonrelativistic matter,
bar rho_m~a_J^-3 and fixed G_eff, hence A_m~a_J^-1.

With Fourier contrast power P_delta and the mode transfer T_k,

  Gamma_E(k) proportional to
      W a_s^-2 a_J^-2 T_k P_delta(k),

and therefore

  d ln Gamma_E/dln a_J =
      -2p_v-2
      +dlnW/dlna_J+dlnT_k/dlna_J
      +dlnP_delta/dlna_J.

To meet the manuscript's effective pressureless benchmark,
the CONDITIONAL requirement is

  d ln[W T_k P_delta]/dln a_J = 4p_v-1.

If W and T_k were fixed, the power spectrum would require

  Radiation p_v=2:       P_delta~a_J^7;
  Matter p_v=3/2:        P_delta~a_J^5;
  Mature coasting p_v=1: P_delta~a_J^3.

These are highly restrictive growth benchmarks for the
specified quasi-static curvature energy; they are
NOT requirements on generic independent geometric modes,
nor proof that no alternative SDMC dynamics can be dustlike.

In the illustrative handoff W->0 at late times,
Gamma_E(forced selector)->0 regardless of a finite
density contrast spectrum. Surviving free mode energy
and direct other scalar/matter couplings require a
distinct and complete dynamical investigation.

## 483. Local curvature and global curvature are different quantities

A spatially flat FLRW background has K_FLRW=0 while
inhomogeneous Jordan or structural curvature may remain
nonzero. The Fourier k in the quadratic selector and in
the Poisson equation is NOT the global FLRW K index.
The exact Planck boundary's material distribution and
its long-wavelength spectrum are not fixed simply by
the structural normalization rho_v,P=rho_P.

A scalar mode density covariance, a nonzero relaxed
structural potential energy, or a network of local wells
cannot be directly assigned to the global Friedmann
source rho_g. Such identification needs the full
ADM Hamiltonian and spatial metric stress, both Jordan
potentials, and a complete averaging prescription.

In particular, an induced q_k and curvature energy
need not provide a new attraction on observable scales.
The perturbed Jordan gravitational potential and
its coupling to photons and matter must be calculated.

## 484. Numerical tests and resulting bottleneck

Implemented:
  sdmc/pg_redistribution_spectrum.py
  tests/test_pg_redistribution_spectrum.py

Checks:
 - exact M_u material repartition and normalized Fourier source;
 - N_eff for equal and unequal mass fragments;
 - independent uniform mode variance 1/N_eff;
 - exact modulated one-point density amplitudes;
 - truly pair-correlated clustered-power enhancement;
 - core/dispersed partition with no created mass;
 - bounded structural static transfer 0<=T<=1;
 - finite-difference minimization of E(ell,u);
 - total relaxed energy equals curvature plus relative energy;
 - Poisson k-four cancellation at nonzero k;
 - no forced selector energy after W=0;
 - correct geometric-activation derivative and dust target;
 - exclusion of the k=0 FLRW spatial curvature mode;
 - no inferred CDM normalization or claim of observable
   gravitational success.

This closes the algebraic *conditional* chain
   mass redistribution -> density contrast spectrum
   -> externally supplied Jordan Poisson potential
   -> relaxed structural curvature / local energy.

It does NOT close the dynamical, covariant chain
   initial primordial microphysics
   -> correlated matter spectrum and mapper dynamics
   -> fully self-consistent F/G2/G3/ADM Jordan potentials
   -> independent physical geometric rho_g, P_g, J_g
   -> no-CDM CMB, growth and lensing observations.

The most decisive next derivation is a full matter–metric–
structural scalar source about a nonzero excited state,
with W allowed to switch off, without imposing the
empirical cold-matter normalization.


## 485. A temporary selector can leave a surviving free structural mode

The quasi-static energy-minimization result E_min->0 as W->0
concerns the instantaneous *forced stationary branch*. It
DOES NOT state that every dynamically excited structural mode
must disappear when the external source is switched off.

For B>0, the full reduced scalar equation of Manuscript B
Section 106 is

  ddot ell + (3H_J+dot B/B)dot ell
    +(omega0²+Lambda)ell=-3Lambda u(t),

for constant f_L; more generally add 2dot f_L/f_L.

Write its solution as

  ell(t)=ell_h(t)
          -3∫[ti,t] G_ret(t,t')Lambda(t')u(t')dt',

where G_ret is the retarded Green function of the full
time-dependent linear operator, with the prescribed
initial homogeneous conditions. The source creates an
excitation during the transition. When W and Lambda
vanish afterward, this leaves free oscillatory data
ell(t_end), dot ell(t_end), unless they accidentally cancel.

An initial state ell=dot ell=0 can therefore develop a
nonzero free trace even though the final forced static
response is zero. The kinetic and mass restoring terms
still store energy. This is a possible *transition imprint*,
not a conserved separate material inventory or a claim
about dark-matter abundance.

## 486. Decelerating radiation-tracking numerical experiment

To avoid imposing inflation or a de Sitter background,
the test uses an explicitly decelerating Jordan radiation
tracker with dimensionless time:

  a_J(t)=(1+t/t0)^1/2,
  H_J=1/[2(t+t0)],  q_J=1;
  a_s=a_J^p, p=2,
  a_s=(1+t/t0) (structural coasting).

The assumed smooth activation is

  B(t)=[1+exp(-(t-t_center)/width)]^-1,
  W(t)=4B(1-B),

with t0=5, t_center=5, width=.8,
eta_K=.025, f_L=1, m_L²=1, c_L²=.1,
comoving k=2, and a weak *externally specified* Jordan
metric perturbation u(t)=.01/a_J(t).

The final source is deliberately a weak-field subhorizon
matter-motivated ansatz; a complete radiation-era Jordan
metric, photon-baryon pressure and mapper perturbation
are NOT solved here. The model begins at a finite time
with B>0 and does not extrapolate the classical action
to t<t_P or impose its equation at the singular B=0 limit.

The initial mode ell(0)=dot ell(0)=0 is numerically
integrated using the exact reduced time-dependent ODE.

Results from the normalized example:
 time t=2: ell ≃-0.02087;
 time t=5: ell ≃-0.002665, W=1;
 time t=10: ell ≃-0.0003188, W≃.00769;
 time t=20: ell ≃+0.0005831, W≃2.88e-8;
 time t=60: ell ≃-0.0002900, W≈0.

The amplitude becomes small but is not immediately
identically zero after selector deactivation. The free
mode remains oscillatory and is damped by expansion.

For the stable mass-dominated tail, the free energy is

  E_free=K_l/2 [dot ell²+omega0² ell²],

and the late-time comoving quantity E_free a_J³ is
approximately constant in the adiabatic regime:

  t=20: E_free a_J³ ≃6.89e-7;
  t=30: E_free a_J³ ≃7.00e-7;
  t=40: E_free a_J³ ≃6.81e-7;
  t=50: E_free a_J³ ≃6.70e-7;
  t=60: E_free a_J³ ≃6.77e-7.

This demonstrates the expected conditional reduced
oscillator behavior E_free~a_J^-3 while retaining
radiation-like Jordan deceleration and structural coasting.

Control tests:
  - if the external perturbation u is zero, ell stays zero;
  - if W=0 throughout, ell stays zero from zero initial data;
  - doubling u doubles ell and quadruples free mode energy;
  - after W->0, the forced static mode vanishes while the
    dynamically excited free mode can remain.

The example does NOT establish that physical cosmic
Planck-to-radiation thermalization generates this u(t),
that the surviving mode has the correct initial spectrum,
or that this geometry changes the measurable Jordan
gravitational potentials in the correct way.

## 487. Implication for the primordial-well hypothesis

The gravitational-well analogy now separates three elements:
  1) material redistribution builds inhomogeneous Jordan sources;
  2) structural–Jordan activation can excite trace modes during
     handoff; and
  3) a free stable mass-restoring mode can outlast that activation.

Thus the absence of a *persistent forced curvature-selector source*
does not rule out a *surviving primordial structural excitation*.
This is compatible with the earlier Part III adiabatic result
that mass-restoring geometric oscillations may dilute as a_J^-3.
It does not prove that their averaged complete Jordan pressure
vanishes, or that their gravitational potential growth and
lensing reproduce the no-CDM observational requirements.

The unresolved physical task is to derive the primordial mode
spectrum and amplitudes from the genuine thermal/structural action,
then re-solve the F/G2/G3 plus supporting-space and material
perturbation constraints on the same background without any
reintroduced CDM-normalized source.

Additional implementation and numerical controls:
  sdmc/pg_handoff_survivor.py
  tests/test_pg_handoff_survivor.py
