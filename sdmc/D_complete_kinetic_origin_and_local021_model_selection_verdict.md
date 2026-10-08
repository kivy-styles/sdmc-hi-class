# SDMC D-sector one-pass physical-origin and local021 verdict

**8 October 2026.** Scope: kinetic baseline and its D0 correction, transition shape, radiation and mature endpoints, action-level identifiability, independently rebuilt Planck/DESI likelihoods and statistical standing against local021.

**Result:** Calculations and source-supported tests closed. A unique microscopic kinetic transition and full jointly marginalized local021 family evidence remain **unresolved by the supplied equations**. No claimed derivation or Bayes factor is fabricated.

## Source and conventions

Primary manuscripts: SDMC_A_COMPLETED(7).pdf pp. 640–645, 666–669; SDMC_Structural_Radius_Jacobian_Investigation_Appending_Continuation_VI(2).pdf secs. 122–129. Repository data: sdmc/late300_candidate.json, sdmc/run_late300_linear_covariant.py, GitHub runs documented below.

The author's correction is preserved: **D0 is an additive correction to Dfloor, not a second unrelated kinetic origin.** Use x=ln a for the computational e-fold variable, **not** the distinct structural lapse N in the structural-clock manuscript.

Dfinite(x)=Dfloor+D0[S_D(x)]^pD, S_D=[1+exp(-(x+ln(1+zc,D))/wD)]^-1.

The force-derived Planck-mass window zc,F=4.0340775029, wF=0.3292537707 differs from the historical kinetic correction window zc,D=3.9278763885, wD=0.3311413396. The tested F-only patch retains the latter. Deriving the former does not derive the latter.

## 1. Canonical radiation floor

Assuming a canonical scalar radiation tracker, negligible early braiding, constant Planck mass F and wphi=1/3, rho_phi=3 F H^2 fr and phi_dot^2=rho_phi+p_phi=(4/3)rho_phi. Consequently D_rad,can=phi_dot^2/(F H^2)=4fr. With fr=4/lambda_e^2, this gives **D_rad,can=16/lambda_e^2**.

The historically specified lambda_e=18.40625 gives D_rad,can=0.047226890271848634. The accepted finite Dfloor=0.05181542627513409 is larger by 9.7159393% (ratio 1.097159393237). The old reduced-model rK=0.02553009 was **inverted from the selected floor**, not independently predicted from the full Horndeski coefficients. The finite reconstructed action generally contains a nonzero G3, nonminimal G4, and sign-changing k1; a positive-k1,k2 reduced prethermal identity cannot be transplanted uncritically.

The canonical floor is thus **conditionally derived from the radiation tracker**. This calculation neither selects lambda_e nor proves that the historical accepted full noncanonical Dfloor has that value. The equation fr=4/lambda_e^2 returns different valid fr for different input slopes, absent an independent radiation structural-fraction prediction.

## 2. Finite D0 correction

Historical late300: Dfloor=0.05181542627513409; D0=0.34231919445927034; pD=1.

S_D(0)=0.9919691128; correction today=0.3395700676; D(0)=0.3913854939; **finite ansatz's late asymptote Dfloor+D0=0.3941346207344044**.

At the kinetic transition midpoint xc,D=-ln(1+zc,D), the slope is dD/dx=D0 pD/[2^(pD+1) wD], hence D0=[2^(pD+1)wD/pD] Dprime(xc,D). This extracts D0 *from a preselected profile*; an independent equation predicting its slope is still missing.

The corrected hi_class subluminal scan imposes approximately D0>=0.019188325854 under a fixed F background, historical D window and sound-speed numerator. This is a **lower bound**, not a selector of the accepted 0.342319 value. Positive D admits many pD>0; linear pD=1 is the simplest nonzero Taylor term, but not a microscopic theorem.

In the reconstructed linear-G3 Horndeski action G2=k1(phi)X+k2(phi)X^2-V(phi), G3=g(phi)X and G4=F(phi)/2, consider delta G2=mu(phi)[X-Xbg(phi)]^2. On the background, deltaG2=deltaG2_X=deltaG2_phi=0, while deltaG2_XX=2mu and deltaD=4mu Xbg/F. Therefore a family of covariant actions reproduces the **same selected background** while changing kinetic response. This is an action-level non-identifiability proof, not proof every alternative is perturbatively stable. GitHub verification: run 37803894563.

## 3. Mature D=2 and finite-to-mature matching

The manuscript's mature canonical infrared branch has D_infinity=2, cs2_infinity=1, G3_infinity=0, G4_infinity=F_infinity/2, under the author's adopted mature closure conditions. Its long-horizon background and scalar-perturbation audits were passed for tested positive C3 branches. But the manuscript states that the microscopic mechanism and the unique finite C3/handoff route into this endpoint remain unresolved.

**Do not set D0=2-Dfloor.** The finite late300 logistic law extrapolates to **0.39413**, not 2. These are different stages of SDMC evolution. Obtaining the complete physical D(x) requires a derived covariant coefficient flow/transfer across that later stage.

Strong non-uniqueness calculation: choose a smooth step U(x;b,c), identically 0 for x<=b, identically 1 for x>=c. On b<x<c, write t=(x-b)/(c-b) and U=exp(-1/t)/[exp(-1/t)+exp(-1/(1-t))]. This is infinitely differentiable, with every boundary derivative vanishing. Define Dbridge(x;b,c)=Dfinite(x)+(2-Dfinite(x))U(x;b,c).

For **each** b>0,c>b:
- the exact accepted finite history for every x<=0 remains unchanged;
- Dbridge becomes exactly 2 for x>=c, reproducing the mature kinetic endpoint;
- Dbridge is positive throughout and joins smoothly to all orders;
- changing b,c changes the intervening trajectory.

Four explicit pairs were numerically checked: (0.5,3.5), (1.5,6), (1.1,8.5), (3,9). These are **mathematical kinetic examples, not native freely evolving SDMC solutions**. They establish rigorously that finite and mature endpoint constraints plus positivity and smoothness do **not** select the dynamical handoff. Full-background equations, G3/G4 evolution, sound-speed numerator and matter/thermal transfer might rule some out, but no unique independent connecting action was supplied.

## 4. Genuine candidate-specific Planck and DESI results

Six genuinely candidate-specific linear-G3 actions passed covariant replay and full Planck nuisance profiling (run 37807119146):

| D case | Full Planck chi2 | Delta from accepted D |
|---|---:|---:|
| accepted baseline/correction | 2784.825495386080 | 0 |
| D0=2pi/lambda_e | 2784.821140904252 | -0.004354 |
| D0=Xi inverse/sqrt(F0) | 2784.819502467613 | -0.005993 |
| Force D0 | 2784.824684047699 | -0.000811 |
| canonical tracker floor | 2784.559302430430 | -0.266193 |
| combined canonical floor and 2pi/lambda_e | 2784.556171262576 | -0.269324 |

Six independent native-to-covariant replays, background/spectrum stability and high-l Planck Plik-lite results succeeded (run 37817678495). They confirm these are not stale spectra generated by the previously hardcoded D0/Dfloor workflow.

The corrected official DESI DR1 full-shape run 37807287174 finished successfully for all four actions after a data-download retry:

| D case | Covariant SDMC DESI chi2 | Delta vs fixed LCDM control |
|---|---:|---:|
| accepted | 329.541566860942 | -10.990837 |
| D0=2pi/lambda_e | 329.541183024805 | -10.991220 |
| canonical floor | 329.524673102794 | -11.007730 |
| combined canonical floor+2pi/lambda_e | 329.524692885322 | -11.007711 |

The **within-SDMC fixed-cosmology** combined-change values are Delta chi2_Planck=-0.26932412 and Delta chi2_DESI=-0.01687398, total **-0.28619810** (assuming the experiment-specific likelihood factors as implemented). This is small, positive evidence of observational viability of a structural replacement, **not** its physical derivation and **not** an optimized SDMC versus optimized local021 comparison. DESI's fixed LCDM pipeline control is not a newly reoptimized local021.

## 5. What the existing local021 evidence actually measures

Prior **derived-bundle** conditional nuisance-only results, which are **not the same candidate** as canonical-floor+2pi D0:

- Full-Planck fixed-spectrum nuisance nested log-evidence difference (derived minus local021) **-2.1254874 +- 0.6916684**; run 37764900113; nuisance parameters marginalized, cosmological/structural parameters fixed.
- HSC NLA-z fixed-cosmology nuisance nested log-evidence difference **+2.0390447 +- 0.6148358**; run 37764410555; HSC nuisance coordinates marginalized, cosmological/structural parameters fixed.
- Separate HSC full-TATT best-found chi2 difference **-2.66689485** against local021; run 37764359939; both TATT minimizers stopped at their iteration limits, so global-minimum confirmation is missing.

If the two fixed-spectrum evidence likelihoods have no shared nuisance parameters and their priors factor, their formal product gives a combined conditional log BF of roughly **-0.08644**, with quoted propagated integration error roughly **0.9254**. That does not sample the joint cosmological/structural priors and **does not apply to the canonical-floor D replacement**.

There is no valid numerical shortcut from a fixed-spectrum full-Planck and fixed-control DESI chi2 pair to a **joint fully marginalized Bayes factor against independently optimized local021**. One would need a predeclared parameter family/UV prior, common cosmological parameters sampled **once**, compatible Planck+DESI+SN+HSC likelihoods and covariance/nuisance definitions, independent best fits, and converged nested sampling or equivalent.

For any **specific shared likelihood**, Delta AIC=Delta chi2+2 Delta k and Delta BIC=Delta chi2+ln(n_eff)Delta k. In the distinct HSC TATT example, Delta chi2=-2.666895 and n_eff=60: Delta k=0 gives AIC/BIC=-2.6669; Delta k=1 gives AIC=-0.6669/BIC=+1.4274; Delta k=2 gives AIC=+1.3331/BIC=+5.5218. These are sensitivity scenarios, **not** the measured independence rank. For the within-SDMC D replacement, if it somehow required one additional free fitted parameter, its Delta AIC relative accepted under the fixed Planck+DESI screen would be -0.2862+2=+1.7138. Neither credit nor debit follows without an action-level origin and independent prior choices.

## 6. Final status: full one-pass investigation

| Target | Result |
|---|---|
| Verify D0 corrects Dfloor | **Verified** |
| Canonical tracker Dfloor=16/lambda_e^2 | **Derived conditionally**; full action and lambda_e origin remain open |
| Noncanonical accepted floor enhancement | **Not independently derived** |
| D0 numeric correction | **Not independently derived**; conditional causal lower bound and action degeneracy |
| Kinetic exponent pD and D-switch timing/width | **Not uniquely derived** |
| Future mature D=2 and cs2=1 | **Manuscript closure result** for tested canonical IR branches |
| Finite-to-mature kinetic handoff from endpoints alone | **Not unique**; four exact C-infinity counterexamples |
| True-action Planck checks | **6/6 complete** |
| True-action DESI checks | **4/4 complete** |
| Full joint local021 AIC/BIC and family Bayes factor for D replacement | **Not established**; fixed-parameter comparisons and nuisance-only evidence cannot answer it |

**Verdict.** The quantitative kinetic checks are closed and reproducible. A unique SDMC dynamical evolution from radiation canonical tracker through the finite correction to the mature canonical endpoint **cannot be derived from current inverse-reconstructed action coefficients and endpoint conditions alone**. Claiming that all free structural coordinates have been erased, or that full Bayesian model selection has closed, would not follow from the source material. The next *necessary new physical input* is an independently determined finite scalar action, structural/thermal transfer and UV/initial conditions that jointly select D(x) and its matching, followed by a genuine family-level joint likelihood.

## Reproducibility

Executable audit: sdmc/audit_D_physical_endpoint_and_joint_selection.py. GitHub Action **37830802991** passed eight test groups and archived JSON, conditional Planck/DESI CSV and kinetic handoff curves. Branch: sdmc-late300-D-covariant-true-inputs. Original accepted late300 benchmark unchanged.

Proof runs: 37803894563 (action degeneracy), 37804665278 (corrected causal bound), 37807119146 (full Planck), 37807287174 (DESI), 37817678495 (independent replay), 37764900113/37764410555/37764359939 (local021 conditional evidence and TATT).


## Addendum: source-guarded canonical covariant action-matching theorem (8 October 2026)

**New completed test:** [GitHub run 37831641596](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37831641596), eight checks successful; source [sdmc/audit_D_canonical_IR_action_matching.py](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-late300-D-covariant-true-inputs/sdmc/audit_D_canonical_IR_action_matching.py). These are **analytic necessary conditions on a possible physical continuation**, not independent forward covariant future evolution.

### Mature canonical structural field in the finite field convention

The manuscript's mature canonical action for sigma and Z_sigma=-1/2(nabla sigma)^2 is

    G2_IR=(2 F∞ Z_sigma-4 F∞ Z_BI)/sigma²,
    G3_IR=0; G4_IR=F∞/2; D_IR=2; c_s2_IR=1.

On the mature coasting branch p=d ln sigma/d ln a->1, define sigma=C*exp(phi), phi=ln a, C>0. Then Z_sigma=sigma² X_phi, with X_phi=-1/2(nabla phi)^2. This implies

    G2_IR(phi,X_phi)=2 F∞ X_phi-(4 F∞ Z_BI/C²)*exp(-2 phi).

The exact same-field canonical IR coefficients required for the mature action are

    k1_IR=2 F∞,  k2_IR=0, g_IR=0,
    V_IR(phi)=(4 F∞ Z_BI/C²)*exp(-2 phi).

For phi=ln a coasting, X_phi=H²/2 and alphaK=2X_phi*G2_X/(H² F∞)=2; alphaB=0, so D_IR=2. These are new explicit coefficient boundary conditions, not a guessed kinetic interpolation.

For force-derived AF=0.020520, F∞=exp(AF)=1.020731982678702, so k1_IR=2.041463965357404. Independently validated derived-reference coefficients currently have k1(today)=-0.0913078779588896 and k2*X(today)=+0.0825143261780817. If the same field convention and continuous operator coefficients apply, k1 must cross zero before the canonical IR limit, while quadratic kinetic effects must decouple. The sign change **does not by itself imply a ghost**; the full Horndeski scalar kinetic and sound-speed functions are decisive.

### Necessary on-shell No-Slip G3 decay on coasting

The selected linear-G3 reconstruction gives *on its target trajectory* g(phi)=-F_phi/H² for

    F=exp(AF S_F(phi)), S_F=[1+exp(-(phi-phi_c)/wF)]^-1.

Assume a future coasting history H² proportional to a^-2 and phi=ln a, on which this same No-Slip trajectory relation persists. Then F_phi is proportional to a^(-1/wF), giving the *required* asymptotic on-shell cubic coefficient

    g(phi) proportional to a^[2-1/wF].

The coefficient g decays if **wF<1/2**. For derived wF=.32925377073762596, 2-1/wF=**-1.03717098747**. Direct logistic-derivative sampling from phi=10 to 12 reproduces the exponent with ~1e-15 difference, and alternative widths .45, .50 and .60 correctly classify decay/borderline/growth.

**Critical scope:** This relation was reconstructed from a selected background, **not independently specified off-shell as a future-action law**. One cannot infer that the original fixed coefficient table will follow the future coasting H. The reconstruction source sdmc/run_late300_linear_covariant.py builds CubicSpline coefficient tables from the fitted target grid; its C interpolator extends beyond that grid by reusing the **last cubic segment**. This is a numerical extrapolation, not physical action closure. No automatic convergence to the mature state was proven.

### Inter-branch gravitational normalization

The derived finite force transition has F∞=1.020731982678702, but the older structural-clock manuscript's matured branch reports F∞ approximately 1.023885427695. Their ratio differs from one by **0.003089395718**, approximately 0.309%. An action with one fixed gravitational normalization must reconcile these calibrated branches or derive a new future F evolution; simply matching D cannot satisfy both exact limits. Different calibrated model branches are being compared, so this is a theoretical matching issue, **not** a direct observational exclusion.

### What is still needed before claiming closed first principles

The radio-era canonical tracker gives the conditional baseline and the mature canonical action supplies necessary operator limits. A unique entire kinetic flow still needs independently supplied functional coefficients k1(phi),k2(phi),V(phi),g(phi),F(phi), microscopic energy transfer and boundary conditions, and a genuine forward (not inverse-matched) homogeneous/perturbative evolution satisfying D>0, 0<=c_s2<=1, and consistent No-Slip. The source-defined out-of-domain cubic continuation is insufficient. No full family Bayes factor against local021 follows from this analytic matching audit.
