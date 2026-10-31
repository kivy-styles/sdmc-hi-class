# SDMC Manuscript B — Part III: Excited Supporting-Space ADM Sources and Long-Wavelength Geometry

## 458. Structural wells and the source of the Jordan gravitational constraint

The preceding investigation derived the structural trace's kinetic,
gradient, mass-restoring and curvature energy; its Jordan-time energy
exchange law; and its joint locked-flat quadratic action with the
rolling F/G2/G3 scalar.

The next issue is to derive how an *already excited* supporting-space
well acts as a gravitational source when its ADM lapse and shift are
varied. The calculation must preserve the distinct SDMC quantities:

  R_s=c t_s, a_s=S/S0,
  Jordan a_J, Jordan proper time t_c,
  M=(a_s/a_J),
  material inventory M_u (distinct from M_inferred),
  supporting space V_s (distinct from the distribution V_u).

This treatment makes NO identification of local curvature dilution
with global FLRW, nor an identification of the structural-mode energy
with observed M_inferred. No CDM-like density normalization is inserted.

## 459. The restricted quadratic supporting-space ADM action

For an existing scalar trace ell in the supporting space, the quadratic
action within the accepted positive-energy trace branch is

  S_red=∫dt_c d³x N sqrt(h) [
    K/2 ((dot ell-N^i ell_,i)/N)²
    -K m_L² ell²/2
    -K c_L² h^ij ell_,i ell_,j/2
    -V_K ],

  K=B f_L²/3,
  h_ij=a_J² delta_ij on the source-evaluation slice.

The spatial supporting metric is

  gamma_ij=a_s² exp[2q]delta_ij,
  q=ell/3+zeta,

where zeta is the Jordan spatial curvature perturbation in unitary
gauge. At quadratic order in the structural curvature,

  V_K=D [sum_ij(q_,ij)²+(Delta q)²],
  D=2 W eta_K/a_s⁴.

This expression derives from the 3D spatial
Riemann-squared operator of gamma, and includes both the
Hessian and Laplacian contractions. The spatial coefficients use
Jordan a_J for the ell-gradient term and the structural a_s for
the supporting curvature, while the ADM velocities use t_c.

This action is restricted to the scalar trace. Its exact nonlinear
curvature-squared metric dependence, traceless structural tensor Q_ij,
the ADM scalar–tensor F/G2/G3 couplings, and higher field/metric
vertices are not present in this quadratic truncation.

## 460. First ADM lapse and shift variations

Vary S_red with respect to the lapse N at N=1 and N^i=0,
holding the structural scalar and spatial h fixed. The result
per Jordan physical volume is

  -(1/sqrt(h)) delta S_red/delta N = rho_well^red,

  rho_well^red = (K/2)[dot ell²+m_L² ell²+
                      c_L²/a_J² |grad ell|²]
                + D[sum_ij(q_,ij)²+(Delta q)²].

This is a genuine first variation of the **restricted** ADM
supporting-space action, not a postulated extra mass density.
For a nontrivial quadratic state on the positive branch,
rho_well^red>=0 pointwise.

Similarly, variation with respect to the CONTRAVARIANT ADM shift N^i,
at the flat background, gives

  -(1/sqrt(h)) delta S_red/delta N^i
      =j_i^red=K dot ell ell_,i.

Consequently an evolving spatial well can have a nonzero
momentum source even when its mean trace vanishes.

The physical full Jordan stress T^mu_nu still requires
variation of the COMPLETE action, including structural curvature
beyond quadratic order and the mixed metric/structural clock
operators. These lapse/shift sources must not be mistaken for
a completed covariantly conserved total fluid.

The numerical implementation independently differentiates the
unexpanded lapse/shift dependence of the *truncated* ADM action
for arbitrary periodic spatial lapse/shift modulations, and verifies
that the results equal the local source Fourier projections.

## 461. Physical curvature source versus homogeneous backreaction

Variation of the quadratic curvature functional with respect to
the Jordan spatial-curvature perturbation zeta yields

  (1/sqrt(h))delta S_K/delta zeta
    =-4 D Delta² q
    =-8 W eta_K/a_s⁴ Delta²q

for a periodic 3D domain, with appropriate spatial Fourier signs.
This is the curvature OPERATOR source linear in q.
Its spatial average vanishes identically on a periodic domain:

  <Delta²q>=0.

Yet the averaged lapse energy is quadratic and generally nonzero:

  <rho_well^red> >0

for an excited inhomogeneous well.

This distinction is fundamental:
- a nonzero spatially averaged quadratic energy enters the
  homogeneous lapse tadpole at second order;
- spatial metric curvature variation produces a linear
  high-derivative operator on the nonzero trace;
- metric stress at order ell² requires *cubic* terms of
  the full off-shell action and its nonlinear supporting curvature.

The absence of a homogeneous mean trace, <ell>=0, does NOT
imply a vanishing average of its energy, stress or momentum
variance.

With the background Friedmann equation expressed using SDMC's
existing Jordan conventions and folding Fdot terms into rho_X,

  3 F H_J² = rho_b+rho_r+rho_X+rho_g,

a nonzero geometric source requires re-solving the rolling
F/G2/G3 and supporting-space background constraints. In particular,
the earlier scalar-only solution, which solved

  3 F H_J²=rho_b+rho_r+rho_X,

cannot remain an exact same-background solution after
a nonzero average geometric source has been appended.
The implementation includes a *bookkeeping residual only*,
not an on-shell recalibration or empirical CDM normalization.

## 462. Long-wavelength geometric-energy modulation

Take a nonzero short-wavelength supporting-space configuration

  ell_0(x)=A[1+epsilon cos(k_L x)] f_short(x,y,z),

where f_short has zero spatial mean and k_L << k_short.

Because the source energy is quadratic in ell, it develops
a long-mode Fourier component at k_L despite <ell>=0.
This is a mathematical coarse-grained source generated by the
existing field, not a new matter species.

The exact Gateaux variation of its quadratic density is

  delta rho= K [
      dot ell_0 dot(delta ell)
      +m_L² ell_0 delta ell
      +(c_L²/a_J²)grad ell_0.grad(delta ell)]
     +2D[
       sum_ij q_0,ij delta q_,ij
       +(Delta q_0)(Delta delta q)],

where q_0=ell_0/3 and delta q=delta ell/3 for
fixed Jordan geometry.

The implementation verifies this derivative independently by
centered finite differences of the full quadratic density.
In the exactly unexcited trace state ell_0=dot ell_0=0,
the first density variation vanishes. At nonzero excited
ell_0, delta rho can be nonzero at linear order in
delta ell. Thus the existing structural inhomogeneous
background matters for any physical perturbation closure.

## 463. The spatial stress and isotropy question

The existing trace-gradient action has an explicitly calculable
traceless spatial stress contribution,

  pi^i_j,gradient= K c_L²/a_J² [
        ell_,i ell_,j
       -delta^i_j |grad ell|²/3 ].

This is only the canonical gradient part of the physical stress.
The structural curvature operator's full tensor metric
variation, the traceless structural metric Q_ij and rolling
scalar–tensor mixing remain outstanding.

An equal-axis periodic short-mode configuration has zero
AVERAGED gradient anisotropy in the numerical test, but a
nonzero local gradient stress. This does not prove statistical
isotropy of a cosmological ensemble or an acceptable
CMB/lensing anisotropic stress.

The proper gravitational test requires the complete
T^i_j and its scalar/traceless decomposition, including
its gauge-invariant perturbations.

## 464. A geometric-amplitude consistency check

For a static single cosine trace ell=A cos(kx), curvature-only
and m_L=c_L=0, the averaged source density is

  <rho_K>=2 W eta_K k⁴ A²/(9a_s⁴).

If a geometric background density rho_g were derived
independently from the full SDMC action and cosmological
initial conditions, reproducing it by this restricted
mode would require

  A²=9 rho_g a_s⁴/(2 W eta_K k⁴).

This is an identity, not a CDM-calibrated amplitude choice.
The quadratic truncation further requires |A|<<1,
and a physical EFT cutoff must limit the allowed k.
At long wavelength the k⁴ derivative suppression could
otherwise require large A; high k could exit the EFT regime.

The amplitude, primordial spectrum and coupling
W eta_K have not yet been fixed by the full action.
The construction cannot claim to predict the required
cosmological geometric abundance.

## 465. Numerical tests and significance

The new branch tests:
- exact lapse source versus independent action variation;
- density Fourier mode as the lapse perturbation source;
- exact canonical momentum source versus shift variation;
- structural curvature variation versus bi-Laplacian operator;
- equality of earlier and new 3D structural density forms;
- quadratic density/momentum scaling with the well amplitude;
- linear spatial curvature operator scaling with amplitude;
- first density variation around an excited well versus
  direct energy finite differences;
- no first density variation on the unexcited branch;
- static wells with positive energy and negative restricted pressure;
- zero mean momentum for a standing wave and nonzero momentum
  for a traveling wave;
- gradient anisotropy locally versus equal-axis averaged isotropy;
- direct quadratic long-wavelength density modulation;
- nonzero source incompatibility with an unmodified
  scalar-only Friedmann bookkeeping relation;
- fixed-S restricted Jordan scale-pressure derivative.

No measured mass abundance, inferred M_u, conventional
omega_c or CDM-like gravitational source was inserted.

## 466. Boundary of the present result

This calculation demonstrates that an EXCITED structural well
carries explicit reduced ADM Hamiltonian and momentum sources,
and can generate a long-wavelength density variation from its
underlying short-scale structure.

It does NOT yet prove that the *physical* SDMC global rho_g
exists with an acceptable pressureless, conserved,
cold-clustering cosmological evolution. That requires:
1. the full off-shell nonlinear supporting-space ADM
   action through the cubic metric/trace terms;
2. the complete F/G2/G3 background and perturbation
   constraints with a nonzero geometric state;
3. the physical tensor spatial stress, momentum, and
   curvature coupling, including traceless structural modes;
4. the initial structural power spectrum and positive
   normalization obtained independently of CDM;
5. gauge-invariant Newtonian and Weyl potentials, and
   observable CMB, lensing and structure likelihoods.

Files:
  sdmc/pg_inhomogeneous_adm_source.py
  tests/test_pg_inhomogeneous_adm_source.py
  .github/workflows/sdmc_pg_inhomogeneous_adm_source.yml
