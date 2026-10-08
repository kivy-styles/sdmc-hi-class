# SDMC late300: genuine Dfloor and D0 correction validation

**Date:** 8 October 2026. This report documents the corrected kinetic-parameter workflow results and remaining theoretical and likelihood limitations.

## Original error and repair

D0 **is the correction to Dfloor**, not an unrelated scalar kinetic sector:

D(N) = Dfloor + D0 S_D(N)^pD.

The derived F-force transition (zc,F=4.034077502899133, wF=0.32925377073762596) differs from the historical kinetic correction transition (zc,D=3.927876388467848, wD=0.33114133956842123). They must not be silently conflated.

Older full-Planck/DESI wrappers modified the native target D0/Dfloor but called a separate legacy covariant-reconstruction script with *hardcoded* kinetic coefficients. Those outputs cannot validate the replacement action. New true-input workflows rebuild the action coefficients from each candidate-specific target and check the actual free evolution.

## Genuine six-candidate full Planck (run 37807119146)

Nuisance-parameter-profile chi2, fixed chosen cosmological coordinates (not optimized model-family Bayesian evidence):

| Candidate | D0 | Dfloor | Full Planck chi2 | Delta chi2 vs accepted |
|---|---:|---:|---:|---:|
| accepted D | 0.34231919445927034 | 0.05181542627513409 | 2784.825495386080 | 0 |
| D0=2pi/lambda | 0.34136151074659893 | 0.05181542627513409 | 2784.821140904252 | -0.004354 |
| D0=Xi inverse/sqrtF | 0.342002216921491 | 0.05181542627513409 | 2784.819502467613 | -0.005993 |
| Force candidate D0 | 0.3419932122356156 | 0.05181542627513409 | 2784.824684047699 | -0.000811 |
| Dfloor=16/lambda^2 | 0.34231919445927034 | 0.047226890271848634 | 2784.559302430430 | -0.266193 |
| Combined Dfloor canonical and D0=2pi/lambda | 0.34136151074659893 | 0.047226890271848634 | 2784.556171262576 | -0.269324 |

Each job fetched and validated its own distinct reconstructed-action spectrum artifact from run 37806771106 before profiling full Planck. All six jobs succeeded.

## Independent cross-check (run 37817678495)

This branch's original proof-of-input script is sdmc/replay_D_covariant_true_inputs.py and its GitHub Actions runner is .github/workflows/sdmc_D_true_input_covariant_planck.yml.

It genuinely regenerates the target using each D0/Dfloor pair, reconstructs G2=k1X+k2X^2-V, G3=gX and G4=F/2, checks recovered D0 and Dfloor, free-evolves the action and compares the CMB/matter perturbations. All **six jobs passed**. The free covariant H(z) agreed with the target at ~5e-10 to 7e-10 for z<=100, with distinct reconstructed kinetic coefficients and sound-speed stability for all six.

Real high-l TTTEEE Planck Plik-lite (calibration prior included, *not* full Planck):

| Candidate | Plik-lite chi2 on free covariant action |
|---|---:|
| Derived reference, D0=0.341993212 | 591.838436036280 |
| Historical D0=0.342319194 | 591.832548744154 |
| D0=2pi/lambda | 591.834184410076 |
| D0=Xi inverse/sqrtF | 591.839763295824 |
| Dfloor=16/lambda^2 | 591.582542065711 |
| Combined canonical floor and 2pi/lambda correction | 591.583043830562 |

This independent result supports the **direction** of the true full-Planck finding. It must not be numerically combined with full-Planck values because the likelihood coverage differs.

## Genuine DESI DR1 full shape

Run 37807287174 used candidate-specific rebuilt covariant spectra and the official DESI DR1 full-shape data, compared with the workflow LCDM control:

| Candidate | SDMC DESI chi2 | LCDM control chi2 | Delta |
|---|---:|---:|---:|
| Canonical floor | 329.524673102794 | 340.532403392675 | -11.007730 |
| Combined canonical floor + 2pi/lambda | 329.524692885322 | 340.532403392527 | -11.007711 |

The initial accepted-D and D0=2pi/lambda jobs failed **before likelihood evaluation** on external DESI HDF5 downloads (timeouts on DESI host, 404 on mirror). A failed-jobs-only rerun was started. These failures are unrelated to physical stability. Verify the rerun job results before declaring all DESI cases closed.

Note: the control is fixed in this workflow; it is not demonstrated to be an independently reoptimized local021.

## Correct physical and statistical interpretation

1. D0 is a time-dependent correction added to a baseline Dfloor. If truly independent dynamical equations predict both endpoint values then D0=Dlate-Dfloor, but these current action reconstructions remain inverse mappings of a selected target.
2. The canonical radiation-tracker Dfloor=16/lambda_e^2 improves Planck slightly and is compatible with DESI in the tested branch.
3. D0=2pi/lambda_e and D0=Xi inverse/sqrtF remain *proposed origin hypotheses*: near-equality and likelihood fit do not themselves derive them from a microscopic theory.
4. No model-family AIC/BIC parameter reduction or full Bayesian evidence credit can follow solely from this numerical replay; the selected action and structural priors must be independently specified and all likelihood parameters fairly sampled.
5. The remaining goal is the independent UV/kinetic coefficient origin, correctly handling the separate D and F transition windows; then optimize and compare each actual model family to local021.

### Sources

- Full Planck: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807119146
- Candidate true-action spectrum input: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37806771106
- Independent six-candidate replication: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37817678495
- DESI full shape: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807287174
