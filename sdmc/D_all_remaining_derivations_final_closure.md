# SDMC D-sector: complete remaining derivations and the microscopic identifiability boundary

**Date:** 8 October 2026  
**Branch:** sdmc-D-forward-canonical-attractor  
**Scope:** The earlier Dfloor + correction D0 tests, finite Horndeski action, radiation tracker, mature canonical SDMC structural-clock action, physical forward evolution, scalar stability, and fair local021 model-selection accounting.

**Executive scientific verdict:** The physically specified *mature canonical action* has now been **evolved forward** (eight independent physical homogeneous initial states) and its coasting attractor derived, including exact global scalar-only attraction, matter/radiation eigenvalues, matter-loaded asymptotic coefficients, the structural-lapse Planck identity, and on-shell cubic decoupling. A separate all-orders action-jet proof establishes why these successful results and the already successful Planck/DESI spectra do **not** determine a unique *finite late300-to-mature operator flow or microscopic kinetic correction*. This is a precise identification of missing physical equations, not a numerical solver failure and not permission to invent UV parameters or a family Bayes factor.

## 1. Observationally accepted finite SDMC kinetic sector

Preserve the author's correction: **D0 is a transition correction to baseline Dfloor**. With n=ln a (which must NOT be confused with the separate manuscript structural lapse N),

    D_fin(n)=Dfloor + D0 * [S_D(n)]^pD,
    S_D(n)=[1+exp(-(n+ln(1+zc,D))/wD)]^-1,
    Dfloor=0.05181542627513409,
    D0=0.34231919445927034,
    pD=1,
    zc,D=3.927876388467848,
    wD=0.33114133956842123.

Today D_fin(0)=0.3913854938989103, finite correction=0.33957006762377623, and the **finite ansatz's extrapolated asymptote** is Dfloor+D0=**0.3941346207344044**. It must NOT be identified with the mature canonical D_infinity=2.

The force-derived transition parameters are **distinct**: AF=0.020520, zc,F=4.034077502899133, wF=0.32925377073762596. The derived F-only workflow **keeps the historical D timing and width**. The current finite action is the on-background reconstructed class

    G2=k1(phi) X+k2(phi) X²-V(phi),
    G3=g(phi) X, G4=F(phi)/2.

The canonical radiation scaling identity is **Dfloor,can=16/lambda_e²=0.047226890271848634 for lambda_e=18.40625**, assuming canonical scalar radiation tracking, negligible early braiding and constant F. No current independent microscopic selector fixes lambda_e, the historical enhancement Dfloor/Dfloor,can=1.0971593932, or the D0 amplitude. The previously cited rK=0.02553009 was inversely inferred from a reduced kinetic model and is not the independently predicted, full action k2 X/k1.

Native full-Planck and DESI true-action replacement replays succeeded; the canonical tracker floor modestly improves the fixed-cosmology profiles. They do not independently derive the finite correction or a model-family complexity discount.

## 2. Canonical mature SDMC action in both scalar coordinates

Manuscript A, pp. 667–669, defines direct structural sigma, Z=-1/2(nabla sigma)^2 and

    G2_IR=(2 F_inf Z - 4 F_inf Z_BI)/sigma²,
    G3_IR=0,
    G4_IR=F_inf/2,
    Z_BI=4pi/(3 tP² S0²).

Use the valid field transformation sigma=C exp(phi), so Z=sigma² X_phi. The same action becomes

    G2_IR=2 F_inf X_phi - V0 exp(-2phi),
    V0=4 F_inf Z_BI/C²,
    G3_IR=0, G4_IR=F_inf/2.

It has k1_IR=2 F_inf, k2_IR=0 and a normalized dimensionless canonical potential slope lambda_IR=sqrt(2), independent of the finite radiation-tracker slope lambda_e. This is a **derived field-coordinate identity**, not a hypothetical new operator.

If F_inf is chosen to equal the newer force-derived profile's asymptote, F_inf=exp(0.020520)=**1.020731982678702**, the same-field canonical k1_IR=**2.041463965357404**. The older canonical manuscript's tested IR branch instead used F_inf≈1.023885427695, a relative difference of 0.0030893957. Both cannot be the *same exact constant-F endpoint* without reconciling the branch normalization or additional physical gravitational running. They are not, however, direct evidence of physical instability; they were calibrated on different branches.

A continuous same-field evolution must connect the currently negative derived-reference k1(today)≈-0.09130788 to k1_IR>0 while suppressing k2 and g. A k1 crossing does **not** imply a ghost if the complete Horndeski scalar kinetic D and sound-speed numerator remain healthy. This crossing cannot be checked with the future continuation of a spline table reconstructed only from past data.

## 3. Genuine forward homogeneous dynamics from the specified mature action

For canonical G2=2 F_inf X - V0 exp(-2phi), define

    n=ln a; u=dphi/dn;
    Om=rho_m/(3 F_inf H²); Or=rho_r/(3 F_inf H²);
    y=V/(3 F_inf H²)=1-u²/3-Om-Or.

Friedmann, Raychaudhuri, scalar and matter/radiation continuity reduce to the **closed forward system**

    epsilon=-d ln H/dn=u² + 3 Om/2 + 2 Or,
    du/dn=-(3-epsilon) u + 3 y,
    dOm/dn=Om*(-3+2 epsilon),
    dOr/dn=Or*(-4+2 epsilon),
    dphi/dn=u,
    d ln H/dn=-epsilon.

Here q=epsilon-1 and, on this constant-F, non-braided canonical branch, D=2u² and c_s²=1. The Friedmann constraint is Om+Or+u²/3+y=1. Scalar kinetic energy is F_inf H²u², so kinetic fraction u²/3.

The scalar-dominated fixed point is

    u=1, Om=0, Or=0, y=2/3,
    epsilon=1, q=0, D=2, c_s²=1,
    H ∝ a^-1, phi=ln a + constant.

This is coasting expansion with w_phi=-1/3. None of these values were numerically inserted as solver targets; they arise from the specified canonical action.

### 3.1 Linear analytic stability and first matter correction

Eliminate y using the constraint. The (u,Om,Or) system has Jacobian at (1,0,0)

        [ -2  -3/2  -1 ]
    J = [  0   -1    0 ]
        [  0    0   -2 ].

Hence eigenvalues **(-2,-1,-2)**, all negative. The radiation mode falls as a^-2, the matter mode as a^-1, and the homogeneous scalar velocity mode as a^-2.

Linear forcing du/dn=-2(u-1)-(3/2)Om-Or, with Om=m1/a+..., gives the particular solution

    u=1-(3/2)m1/a + O(a^-2 log a),
    q=epsilon-1=-(3/2)m1/a+O(a^-2 log a).

This reproduces the exact coefficients in the author's mature structural-clock manuscript for p_structural=u and q. The O(a^-2 log a) resonance is consistent with the coincident radiation and scalar eigenvalue -2. This is **conditional on the already specified IR action**, not a derivation of the finite Horndeski transition.

### 3.2 Global homogeneous stability theorem in scalar-only physical phase space

For Om=Or=0 the field equation simplifies algebraically to

    du/dn=(u-1)(u²-3).

The physical positive-potential scalar-only region has -sqrt(3)<u<sqrt(3). Define L=(u-1)², then

    dL/dn=2(u-1)² (u²-3)<0 for u≠1.

The flow points to u=1 from both sides, so it is **globally attracting within the scalar-only physical interval**. This is stronger than linear stability. The exact separated-variables solution uses the primitive

    I(u)= -1/2 ln|u-1| + 1/4 ln(3-u²)
          + 1/(4sqrt(3)) ln |(u-sqrt(3))/(u+sqrt(3))|,

such that **I(u)-I(u0)=n-n0** for u,u0≠1. Eight independently integrated scalar-only orbits verify the integral, monotonicity and the approach to u=1.

### 3.3 Forward numerical validation: not a prescribed desired curve

Source: [sdmc/run_D_full_forward_IR_and_UV_identifiability.py](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-D-forward-canonical-attractor/sdmc/run_D_full_forward_IR_and_UV_identifiability.py).

The exact first-order equations were integrated via 4th-order Runge-Kutta, Δn=.002, from n=0 to n=12 for **eight distinct physically consistent initial (u,Om,Or)** states (including positive, zero and negative initial scalar velocities). The constraint was independently reconstructed from conserved rho_m, rho_r and the original exponential potential normalization at each step.

Latest CI [37833225296](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37833225296): **8/8 forward evolutions successful**; maximum |u(12)-1|≈2.4962×10^-6; maximum independently evaluated Friedmann residual≈6.8023×10^-12. The script also checked eight scalar-only orbits, seven cubic transition-width cases and five very different Planck clock/stretch normalizations.

**Limit:** The solver starts within the mature canonical action. It does NOT claim that the earlier finite late300 action has been evolved through its unknown operator transition.

## 4. Structural lapse, speed and Planck identity directly from the canonical action

Under scalar domination and coasting, y=2/3 and V=2F_inf H². Since V=4F_inf Z_BI/sigma²,

    H²=2 Z_BI/sigma².
    dot(sigma)=H sigma=sqrt(2 Z_BI).

The manuscript defines structural lapse N_struct=tP S0 dot(sigma) and structural radius R=lP S0 sigma. Therefore

    N_struct,infinity=tP S0 sqrt(2 Z_BI)
                     =sqrt(8pi/3)
                     =2.894405018...,

and

    dot R / c = N_struct,infinity,

because lP/tP=c. The dependence on F_inf, S0 and tP cancels **exactly** for Z_BI=4pi/(3tP²S0²). This is a true action-to-structural-kinematics identity within the adopted canonical Balanced Identity normalization. Five clock/stretch choices were tested in the forward audit; all reproduce the same limit.

This exact mature **structural lapse N_struct** is distinct from both n=ln a and the current structural clock in the finite era. It must not be substituted for a finite-epoch input just because their notation overlaps.

## 5. Cubic No-Slip decay: the right operator-level condition

The past-fitted linear-G3 relation **on its selected cosmological trajectory** is g(phi)=-F_phi/H² and G3=gX, for X=H²/2. Assume a new future trajectory independently approaches coasting H²∝a^-2, phi=ln a, and continues to satisfy No-Slip. Then

    F=exp(AF S_F), F_phi/F=alpha_M,
    g∝a^[2-1/wF] as a ->∞.

Hence the **coefficient g** vanishes if wF<1/2. Derived width wF≈0.32925377 gives coefficient exponent -1.03717098747, as verified by a source-guarded audit [37831641596](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37831641596).

**New stronger identity:** With □phi=-(ddot phi+3Hdot phi)=-2H² on exact coasting,

    G3=gX=-F_phi/2,
    (-G3 □phi)/(F H²)=-F_phi/F=-alpha_M.

For logistic F, alpha_M∝a^(-1/wF), which vanishes for **every wF>0**. Therefore the on-shell *cubic action-density contribution relative to the gravitational kinetic scale* decays even when the coefficient g would grow for wF>1/2. Seven widths [.25,.32925,.45,.5,.6,1.,1.5] passed the exact on-shell relation in CI [37832971681](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37832971681). **Width <1/2 is sufficient for coefficient decay but NOT necessary for the cubic contribution to become small on the coasting trajectory.**

Neither statement is an off-shell theorem that the existing spline-reconstructed finite action must find this future history. The finite action's cubic coefficient was constructed from H and F on a **preselected trajectory**, and the C evaluator extrapolates its final spline segment past the calibrated domain.

## 6. Why no additional unique finite kinetic derivation follows: two exact no-go mechanisms

### 6.1 Positive smooth endpoint matching is nonunique

Earlier verified exact C-infinity family:

    D_bridge(n;b,c) = D_fin(n) + [2-D_fin(n)] U(n;b,c),

where U=0 for n<=b, U=1 for n>=c and for b<n<c,

    t=(n-b)/(c-b),
    U=exp(-1/t)/[exp(-1/t)+exp(-1/(1-t))].

For every b>0,c>b this retains the entire finite past *exactly*, becomes D=2 in the mature future, has positive D and smooth matching to all orders. Distinct b,c select different transitions. Four representative examples passed [37830802991](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37830802991).

These are mathematical counterexamples proving **endpoint conditions and positive D alone do not select a history**, not candidate solutions of the coupled Horndeski EOM. Additional dynamics and initial/boundary physics are required.

### 6.2 Complete background and linear CMB also cannot determine the nonlinear off-shell action

An even stronger exact kinetic-action-jet theorem is tested in [37833138647](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37833138647):

For any smooth selected monotonic background trajectory X=Xbg(phi), deform

    delta G2_m = nu_m(phi) [X-Xbg(phi)]^m,

with integer m>=2 and a smooth nonzero nu_m. Then the new term and its first derivatives vanish on the *entire background trajectory*, so it leaves the homogeneous equations unchanged.

- At m=2, deltaG2_XX=2nu_2 !=0: scalar quadratic kinetic support, D and linear perturbations change despite identical homogeneous expansion.
- At m=3, **all scalar action derivatives up to order 2 vanish on the background**, so the *full quadratic action* and resulting linear scalar perturbations are identical, while deltaG2_XXX=6nu_3!=0: the cubic/nonlinear interactions differ.
- For any m>3, the first new interaction arises at perturbative order m.

The analytic identities were checked for 7 field values and m=2,3,4,5: **28 combinations passed**. They prove that linear Planck, linear DESI and matching H(z) are insufficient to uniquely fix the nonlinear scalar completion. Physical UV selection, symmetries and appropriate nonlinear data are necessary; generic deformations are not being proposed as viable dark-energy theories.

This also means the good finite-action observational fits and the canonical future attractor **do not logically determine the complete finite SDMC microscopic coefficients**. A genuine *forward-derived operator flow* must be supplied independently of the already selected benchmark trajectory.

## 7. Correct statistical comparison — no double counting and no fictitious Bayes factor

Already-completed candidate-specific true-action full Planck nuisances (run 37807119146):

| D choice | full Planck chi2 | difference vs accepted D |
|---|---:|---:|
| accepted Dfloor and D0 | 2784.825495386080 | 0 |
| D0=2pi/lambda_e | 2784.821140904252 | -0.004354 |
| Dfloor=16/lambda_e² | 2784.559302430430 | -0.266193 |
| combined canonical floor plus D0=2pi/lambda_e | 2784.556171262576 | -0.269324 |

Genuine separately reconstructed DESI DR1 full-shape (run 37807287174), with fixed LCDM control chi2≈340.5324:

| D choice | SDMC DESI chi2 | difference vs fixed LCDM |
|---|---:|---:|
| accepted Dfloor and D0 | 329.541566860942 | -10.990837 |
| D0=2pi/lambda_e | 329.541183024805 | -10.991220 |
| Dfloor=16/lambda_e² | 329.524673102794 | -11.007730 |
| combined canonical floor and 2pi/lambda_e | 329.524692885322 | -11.007711 |

The within-SDMC **fixed-parameter** delta chi2 for combined structural floor and D0 prescription vs accepted kinetic baseline is -0.269324 (Planck) + -0.016874 (DESI) = **-0.286198**. This very small shift is not a comparison to *independently optimized local021* and not a model-family evidence calculation.

For separate previously derived-bundle fixed-spectrum nuisance-only tests, full Planck delta lnZ=-2.12549±0.69167 vs local021 (run 37764900113), and HSC NLA-z delta lnZ=+2.03904±0.61484 (run 37764410555). These are different from the canonical-floor D replacement. HSC TATT conditional best-found delta chi2=-2.66689 (run 37764359939); both optimizers hit their iteration ceilings. Do **not** sum distinct dataset/model/profile proxies and claim full posterior evidence.

For any real family-level comparison, a predeclared UV action and priors must define the SDMC family and parameters; the local021 and SDMC families must be jointly sampled under matched Planck/DESI/SN/HSC likelihood conventions (including shared cosmological coordinates, survey cross-correlations where relevant, and dataset-specific nuisances). Nested evidence is

    Z_M = integral pi_M(theta) L_joint(data|theta,M) dtheta,

not exp[-best-fit chi2/2]. If extra dimensions were historically tuned on data, they cannot disappear simply because a designer covariant action is later free-evolved at the selected spectrum.

**Status:** no newly estimated full marginal family Bayes factor is claimed because the current UV coefficient functions, normalized physical priors and final joint model family are not determined. It would be scientifically invalid to fabricate the missing likelihood integral.

## 8. Definitive physical status ledger after doing all derivations supported by current equations

| Required derivation / calculation | What has genuinely been completed |
|---|---|
| Radiation canonical Dfloor | **Derived conditionally:** 16/lambda_e², requires independent radiation-track slope origin and canonical early regime. |
| Historical noncanonical floor | **Not independent:** previous rK solves for the known accepted floor, not for unknown UV coefficients. |
| Historical finite D0 correction magnitude | **Not uniquely derived:** action m=2 degeneracy and broad causal lower bound. |
| D0 transition centre, width and pD | **Not dynamically selected:** historical D window differs from derived F window; leading linear pD is minimality ansatz. |
| Early finite action -> mature IR operator flow | **Not specified:** source k1,k2,V,g,F are inverse-calibrated along selected finite history; endpoint matching is nonunique. |
| Mature canonical action coefficient conversion | **Derived exactly:** k1_IR=2F_inf, k2_IR=g_IR=0, V_IR∝e^-2phi. |
| Canonical future field and expansion | **Forward integrated:** 8/8 homogeneous initial conditions to ln a=12, Friedmann residual ~6.8e-12. |
| Canonical mature attractor stability | **Analytic local and scalar-only global proofs:** eigenvalues (-2,-1,-2) and Lyapunov L=(u-1)². |
| Matter-loaded mature asymptotics | **Derived:** u=1-1.5*m1/a and q=-1.5*m1/a, consistent with manuscript. |
| Structural-lapse and R propagation | **Derived from canonical action + Balanced Identity:** N_infinity=Rdot/c=sqrt(8pi/3). |
| Future cubic No-Slip operator limit | **Conditional coasting derivation:** cubic action density /F H²=-alpha_M->0 for every wF>0; g coefficient itself decays for wF<1/2. |
| Full action uniqueness from linear observables | **No-go theorem:** m>=3 off-shell deformations preserve all linear tests but change nonlinear interactions. |
| Candidate-specific Planck / DESI | **Already completed:** full Planck 6/6; DESI 4/4 after retry, with original candidate-action input validation. |
| Final AIC/BIC and family Bayes vs optimized local021 | **Unresolved:** requires independent physical action/prior specification and full consistent joint sampling. |

**No further algebra can uniquely compute a finite D0, the UV coefficient flow, or a family Bayes factor from the present definitions and selected-target reconstruction alone.** These are *missing independent physical input/data specifications*, not unfinished numerical integrations. The existing mature canonical theory is genuinely forward-solvable and internally stable in its studied regime; the finite-to-mature matching is the scientifically open SDMC problem. The original accepted late300 action and cosmological benchmark remain unchanged.

## Provenance and reproducibility

- [Forward canonical IR equations and eight-seed numerical trajectories](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-D-forward-canonical-attractor/sdmc/run_D_full_forward_IR_and_UV_identifiability.py)
- [Final forward-action CI run 37833225296](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37833225296): JSON verdict, 8-seed trajectories CSV, canonical lapse and exact scalar orbits.
- [Higher-order UV jet non-identifiability proof](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-D-forward-canonical-attractor/sdmc/audit_D_offshell_kinetic_higher_order_degeneracy.py)
- [Jet-proof CI run 37833138647](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37833138647): 28 algebraic cases.
- [Full earlier kinetic endpoint/likelihood audit run 37830802991](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37830802991).
- [Source-guarded coefficient-level finite-to-canonical matching run 37831641596](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37831641596).
- [True-action full Planck run 37807119146](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807119146), [DESI run 37807287174](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807287174), [independent high-l cross-check 37817678495](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37817678495).
