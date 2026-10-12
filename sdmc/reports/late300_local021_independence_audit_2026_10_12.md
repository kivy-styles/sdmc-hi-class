# SDMC Late300 versus optimized local021: structural-parameter independence and model-family complexity audit

**Date:** 12 October 2026  
**Scope:** Reproducible, source-constrained audit of the *published accepted Late300 target function family* and the actual ΛCDM local021 six-coordinate optimization. This is not a new maximum-likelihood search, a new early-universe microscopic derivation, or an independent Boltzmann likelihood. **Canonical late300 and local021 are unchanged.**

## Executive result

- **Both** local021 and Late300 were selected through observation-driven cosmological scans and then dynamically evolved. A frozen evolution uses **zero fresh retuning steps for either model**. That is **not** the same as either model having zero free parameters during its selection.
- Late300 has eleven listed SDMC-specific quantities in `late300_candidate.json`; `Omega_x0` is **exactly algebraically fixed** by flatness and shared background densities, with numerical residual `1.11e-16`. The ten remaining displayed SDMC-specific coordinates represent locally **ten distinguishable directions** in the **chosen phenomenological target ansatz**. Separate numerically normalized sensitivity ranks are **4/4 for H**, **3/3 for F**, **3/3 for D given its window**, and **10/10 jointly for sampled H/F/D**. This is **not** a measurement of their theoretical independence, the number of viable action coefficients, or a likelihood Fisher rank.
- The accepted `linear-G3` covariant action was reconstructed **after generating a Late300 fitted target**. It then genuinely freely evolved and was stable to tested initial-condition perturbations. This demonstrates a robust *covariant realization of a previously selected target*, not an *upstream autonomous determination of the target values*.
- The canonical radiation-tracker derivation gives `D_floor,can=16/lambda_e^2=0.04722689027185` **conditional on** a canonical early tracker and `lambda_e=18.40625`; accepted `D_floor=0.05181542627513` differs by `+9.71594%`. The plausible candidate `D0=2pi/lambda_e=0.34136151074660` differs from accepted `0.34231919445927` by `-0.27976%`. Neither identity is an exact microphysical reduction for the accepted finite action.
- Off-shell Horndeski action constructions prove that a background and even its complete quadratic linear-perturbation action do not uniquely determine nonlinear operator coefficients. In particular, `delta G2 = mu(phi)*(X-Xbg(phi))^2` preserves the selected background but changes quadratic kinetic response; using the **third** power preserves both the background and quadratic perturbations but changes higher-order interactions.
- **Consequently:** `delta_k=0` for additional *post-fit replay retuning* is established. `delta_k=0` for the **full observationally selected SDMC family** remains an open physical and statistical question. The verified ten-dimensional phenomenological target rank must **not** be naively equated with a ten-parameter fundamental action penalty; equally, it cannot be discounted to zero without an upstream action-level constraint or predeclared prior family.

## 1. Equal-footing pipeline audit

| Phase | ΛCDM local021 | SDMC Late300 |
|---|---|---|
| Family specified | Flat ΛCDM perturbation evolution with 6 ordinary cosmological parameters | SDMC early-tracker + two-lobe background, No-Slip F window and finite kinetic D window |
| Observational selection | Explicit 6-dimensional optimizer: `H0, omega_b, omega_cdm, n_s, ln10As, tau_reio` | Staged searches for structural/late parameters and shared ordinary parameters; `lambda_e`, `D0`, `power` also explored in earlier lineages |
| Freeze optimized target | Yes | Yes |
| Independent evolution with frozen implementation | Yes | Yes; an accepted linear-G3 covariant action was reconstructed from target background, then evolved |
| Initial-condition robustness | ΛCDM well-defined evolution | Late300 reconstructed fixed action shown to attract tested perturbations within a finite basin |
| Fundamental origin of optimal parameter values | Cosmological parameters fitted, not fixed by mathematical ΛCDM equations | Shared cosmological parameters fitted; independent source of exact finite SDMC target shape remains unproven |
| Retuning performed after freezing | 0 | 0 |

The accepted action's workflow explicitly invokes the following sequence:

1. Generates `output/linear_cov_target.ini` containing the selected values `A_F`, `z_c`, `width`, `D0`, `power`, `D_floor`, `Omega_x0`, `lambda_e`, `z_t`, `A_late`, `B_late` and fixed shape constants.
2. Runs hi_class to generate `linear_cov_target_00_background.dat`.
3. Calls `sdmc/run_late300_linear_covariant.py`, which **reads** the fitted target background and reconstructs `g(phi)`, `k1(phi)`, `k2(phi)`, and `V(phi)`.
4. Compiles a fixed spline-table covariant action and evolves the scalar and background without prescribing the original H(z) as a runtime override.
5. Replays the CMB/power observables and varies initial conditions while leaving the reconstructed action fixed.

**Verdict:** It is valid to call step 4 a **free covariant evolution**. It is not valid to call step 3 a first-principles determination of the already selected Late300 numerical target.

## 2. Which quantities are algebraically fixed?

The published Late300 values include `H0=69.71482083084993`, `omega_b=0.022083219194622913`, `omega_cdm=0.12299536722293603`, and fixed physical-radiation density `omega_r=4.17998772e-5` in this implementation. Put `h=H0/100` and define

```
Omega_m0=(omega_b+omega_cdm)/h^2;
Omega_r0=omega_r/h^2;
Omega_X0=1-Omega_m0-Omega_r0.
```

Results: `Omega_m0=0.2985060132835426`, `Omega_r0=0.00008600521280791505`, `Omega_X0=0.7014079815036495`, matching the published `Omega_x0` within `1.11e-16` in floating-point arithmetic. This is an **exact algebraic closure** conditional on spatial flatness and the assumed inventory; it does not explain why the primordial cosmology should be flat.

Other exact *conditional* identities are not parameter-origin laws:

- Early canonical tracker: `f_r=4/lambda_e^2`, `f_m=3/lambda_e^2` **given** the slope.
- In canonical radiation tracking: `D_floor,canonical=16/lambda_e^2`, but not equal to the selected full Horndeski floor.
- In the chosen low-redshift two-lobe background, `A_late ≈ q0-q0,flat`; `q0` itself must still be predicted or fitted.
- The two-lobe model's second-order cosmography allows conditional inversion of `B_late` once `j0` and other shape inputs are supplied; it does not predict the jerk.
- Given `D_floor, z_c,D, width_D, power_D` and an independently predicted kinetic midpoint slope `D'(N_c)`, `D0=[2^(p_D+1)*width_D/p_D] D'(N_c)`. The current code does **not** predict this slope without prescribing the kinetic trajectory.
- On the selected No-Slip trajectory, `X G3_X= -F_phi/2` (reduced units); with `G3=g(phi)X`, `g=-F_phi/(2 X_bg(phi))`. If `X_bg` was itself generated from the fitted target, this constrains the reconstructed `g` but does not independently predict `F(z)` or `H(z)`.

A relation among the coefficients is not a demonstrated loss of a free dimension until **all inputs on its right-hand side are themselves fixed by independent fundamental laws or shared fitted parameters**.

## 3. New numerical local-identifiability test

Freeze the **six shared cosmological parameters**, accepted radiation convention, and three specified shape constants `DeltaN=0.5`, `tau_A=0.25`, and `tau_B=1.5`. Use precisely the source-coded functions

\[
E^2(z)=\frac{[\Omega_m(1+z)^3+\Omega_r(1+z)^4+\Omega_X]
[1+A_{\rm late}ze^{-z/0.25}-B_{\rm late}z^2e^{-z/1.5}]^2}
{1-W_t(z)[3+\Omega_r(1+z)^4/(\Omega_m(1+z)^3+\Omega_r(1+z)^4)]/\lambda_e^2},
\]

\[
W_t=\tfrac12[1-\tanh((\ln a-\ln a_t)/0.5)],\qquad
U=\{1+\exp[-(\ln a-\ln a_c)/w_F]\}^{-1},
\]

\[
F=\exp(A_F U),\qquad D=D_{\rm floor}+D_0U^{\nu}.
\]

For each parameter \(\theta_j\), compute the finite-difference column

\[
J_{ij}=\frac{\partial \ln Y_i(z)}{\partial \ln \theta_j},
\qquad Y_i\in\{E,F,D\}.
\]

Each Jacobian column is normalized to unit Euclidean length before singular-value decomposition, avoiding spurious dependence on parameter scale. The test uses 1001 redshift nodes spaced uniformly in `ln(1+z)` from z=0 to 100,000. This is a *functional* derivative matrix, NOT the Planck/DESI/SN likelihood Jacobian.

| Model-function block | Changed inputs | Numerical rank (tol=1e-7 relative) | Interpretation |
|---|---|---|---|
| ln E | `lambda_e,z_t,A_late,B_late` | **4/4** | Four separate background profile directions |
| ln F | `A_F,z_c,width` | **3/3** | Three separate gravitational-profile directions |
| ln D (holding F window fixed) | `D_floor,D0,power` | **3/3** | Three separate kinetic-profile directions |
| Joint ln E, ln F, ln D | All ten | **10/10** | No exact algebraic redundancy within this target ansatz |

Joint column-normalized singular values are approximately

```
[1.74910, 1.09923, 1.04857, 1.00019, 0.98834,
 0.94650, 0.91452, 0.89237, 0.35182, 0.05787]
```

The minimum divided by the maximum is `0.03310`. None vanishes near the accepted point. The rank also remains 10/10 when sampling `0≤z≤1`, but is poorly conditioned there: smallest/largest singular value is `6.20e-6`, compared to `~0.021` for `0≤z≤30`. A low-z-only *empirical* likelihood would therefore have difficulty distinguishing all directions. Actual posterior dimension and degeneracy require the complete observational covariance and parameter priors.

**Meaning:** The code's ten target coefficients can genuinely be varied separately. This DOES NOT refute a deeper theory tying them together. To show that deeper reduction, an upstream action/physical boundary condition must predict relationships equivalent to eliminating ten, six, or another number of degrees of freedom, *without fitting the eliminated values afterward*.

## 4. Two action-level counterexamples to mistaken inference

For the accepted restricted Horndeski class, consider the known background kinetic curve `X_bg(phi)`, where `phi` is a monotonic coordinate along the selected trajectory. Introduce an arbitrary smooth coefficient `mu(phi)`:

\[
\Delta G_2^{(2)}=\mu(\phi)[X-X_{\rm bg}(\phi)]^2.
\]

Along that trajectory, the deformation and its first partial derivatives vanish, so the homogeneous field equations remain unchanged. But `delta G2_XX=2mu`, changing the scalar quadratic kinetic coefficient; on the reconstruction's field convention, `delta D=4mu X_bg/F`. Hence **H(z) alone cannot select the finite kinetic correction**.

Now use

\[
\Delta G_2^{(3)}=\mu(\phi)[X-X_{\rm bg}(\phi)]^3.
\]

Its full action jets through second order vanish along the selected trajectory. It preserves both the background and quadratic perturbation action, yet modifies cubic/nonlinear interactions. Therefore even reproducing full linear spectra after free evolution does not fix the nonlinear completion without an independent microscopic selection law.

The earlier exact GitHub action-jet theorem verifies these counterexamples on multiple finite trajectories. These mathematical deformations need not all be stable or observationally viable; they demonstrate **non-uniqueness without an upstream law**, not that the accepted action is physically wrong.

## 5. Canonical-to-Horndeski D correction: meaningful, but not yet unique

SDMC's accepted finite sector is `D=D_floor+D0 U^power`, with baseline `D_floor=0.05181542627513409`, correction `D0=0.34231919445927034`, `power=1`. In No-Slip restricted Horndeski, `D=alpha_K+6 alpha_M^2`. Because the braiding-squared contribution near the F transition is much smaller than the full kinetic correction, `D0` is not simply the braiding-squared term.

The canonical radiation tracker predicts `16/lambda_e^2=0.047226890271848634` *if* braiding vanishes and F is effectively fixed, 9.71594% below accepted floor. The candidate `2pi/lambda_e=0.34136151074659893` differs from accepted D0 by 0.27976%. Neither numerical closeness is an exact selection equation.

The mature canonical IR action can derive `D_infinity=2` under its coasting closure. The finite Late300 logistic profile, however, approaches `D_floor+D0=0.3941346207344044`, and there are infinitely many smooth positive bridges between a finite profile and the future value of 2. Two endpoints cannot uniquely predict the intermediate height, width or handoff. The prior D-sector audits explicitly prove this and validate several nearby D replacements with distinct covariant actions and observables.

**Necessary physical next step:** Supply the full `k1(sigma), k2(sigma), g(sigma), V(sigma), F(sigma)` and action-boundary conditions **independently of the optimized Late300 profile**, evolve them and check whether `D0`, Dfloor, the force profile and handoff appear as outputs. Merely refitting the action tables from revised target values does not provide this proof.

## 6. Consistent model selection versus local021

The local021 source optimized six ordinary coordinates, not zero:

```
H0, omega_b, omega_cdm, n_s, ln(1e10 As), tau_reio
```

Late300's accepted family has the **same six ordinary inputs**, plus its additional target functions. Exactly how many SDMC parameters remain free **under a completed microscopic theory** has not been established. Quantities held fixed in a final replay may have been selected upstream by an observational search.

For identical likelihoods and genuinely independent fitted parameters,

\[
\Delta{\rm AIC}=\Delta\chi^2_{\min}+2\Delta k,
\qquad\Delta{\rm BIC}=\Delta\chi^2_{\min}+\Delta k\ln n_{\rm eff}.
\]

Use the reported **conditional fair** Planck+DESI+Pantheon+ `delta_chi2 = -0.20248` for a simple AIC *sensitivity example*. This is not a fully marginalized evidence calculation and a rigorous AIC would require truly optimized minima on exactly the same likelihood:

| Assumed `delta_k` | Illustrative `delta_AIC` |
|---:|---:|
| 0 | -0.2025 |
| 1 | +1.7975 |
| 6 | +11.7975 |
| 10 | +19.7975 |

AIC therefore depends critically on whether the structural parameters are **actually derived**. BIC additionally requires a justified independent-data count and regular-model assumptions, which can fail for a highly constrained/nonlinear EFT model. Full family Bayesian evidence requires independent priors on the genuinely free action coefficients, a normalized likelihood and prior integration across the *entire* cosmological parameter family; a fixed-spectrum nuisance-only Bayes factor is not enough.

**Neither** assigning the ten target coordinates mechanically as `delta_k=10` **nor** counting zero simply because the final action is frozen is a demonstrated fair model-family penalty.

## 7. Falsifiable closure criterion and current stopping point

Zero **additional structural-family** parameters could be justified in principle if SDMC's independently specified (pre-fit) complete action and boundary/normalization conditions implied unique *testable* H, F, D trajectories once its six ordinary cosmological parameters were set. Specific checkable prerequisites:

1. Define the action coefficient functions and finite-to-mature canonical/Horndeski matching without referencing preselected Late300 `A_F,z_c,width,D0,Dfloor,z_t,A_late,B_late,lambda_e` or the already selected target `H,F,D` files.
2. Show early structural-clock normalization, tracker fraction and the handoff epoch follow from action/boundary conditions, rather than inherited vacuum-density normalizations or likelihood-selected transition redshifts.
3. Solve Friedmann, Raychaudhuri, scalar-field and linear perturbation equations genuinely forward. Predict and report, *before observational refitting*, the amplitudes/timings of H/F/D and the same Late300 observables.
4. Specify the source of any remaining dimensionless potential amplitude, ratio, kinetic slope or initial-state choices. If those are independently chosen, count or marginalize them consistently; a fixed action may still have independent constants.
5. Construct a preregistered joint local021 versus candidate SDMC inference (Planck, DESI, SN, lensing, BBN as appropriate) with consistent shared six-parameter fitting and a physical structural prior family. Derive **model-family**, not frozen replay, AIC/BIC/evidence.

The earlier October 10 autonomous double-exponential countermodel independently evolved and produced **emergent** handoff/acceleration epochs, but failed to reproduce Late300 and lacked its exact No-Slip braiding relation. This demonstrates that forward dynamics can generate transitions but do not generically pick out the accepted Late300 history.

### Final verdict

- **Verified:** a fixed, dynamically robust, reconstructed Late300 covariant action; zero extra tuning at its replay stage; exact flatness closure of OmegaX; theoretical conditional relations; 10 locally identifiable directions in the existing target functions; six-parameter local021 optimization.
- **Not verified:** a full early-to-late microscopic SDMC action uniquely generating the accepted *numerical* H/F/D transition profiles before fitting, a zero structural-family delta-k, or a genuine fully marginalized local021-versus-Late300 family Bayes factor.
- **Next scientific blocker:** an independently specified finite-era canonical/Horndeski action **and** a microscopic handoff law, rather than additional target-function fit scans. Reconstructing the action from an already optimized target is not sufficient to establish fundamental parameter independence.

## Primary reproducibility references

1. Accepted parameter file: `sdmc/late300_candidate.json`, branch `sdmc-activation-handoff-equations-audit-20261010`.
2. Historical parameter lineage: `sdmc/run_late300_parameter_provenance.py` and `sdmc/run_late300_deep_parameter_provenance.py`, same branch.
3. Full structural derivation: `sdmc/run_late300_full_structural_derivation.py`, successful workflow https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37698426992.
4. Accepted action reconstruction: `sdmc/run_late300_linear_covariant.py`, validated by workflow https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37581925240.
5. Fixed-action initial-state robustness: workflows https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37642049222 and https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37642609079.
6. Action-jet nonidentifiability: `sdmc/audit_D_offshell_kinetic_higher_order_degeneracy.py`, https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37833138647.
7. Kinetic end-state and full-source limitations: `sdmc/D_complete_kinetic_origin_and_local021_model_selection_verdict.md`, same branch.
8. True-action D replacements: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807119146, https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807287174.
9. Independent autonomous countermodel audit: `sdmc/reports/activation_autonomous_countermodel_results_2026-10-10.md`, https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38042838049.
10. local021 six-parameter screen: `sdmc/run_lcdm_joint_ordinary6d_screen.py`, branch `lcdm-local021-fullplik`.
11. New tests: `sdmc/run_late300_independence_audit.py`, `sdmc/test_late300_independence_audit.py`, alongside the workflow file, on the present audit branch.