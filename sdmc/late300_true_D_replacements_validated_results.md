# SDMC late300 true-D-replacement validation — corrected covariant action and likelihoods

**Status: 2026-10-08, live work.** Source of truth is the linked GitHub run artifacts, not historical claims about full Plik/DESI replay.

## Why the old replacement likelihood checks required repair

The historical full Planck and DESI one-coordinate branch workflows replaced D0 or Dfloor in an **upstream parameterized target** and then called \`sdmc/run_r032_partial_zeq_ztbest_covariant.py\`. That reconstruction file includes hardcoded AF, zc, width, D0, power, Dfloor, H0-related scale etc. Thus the downstream physical action/spectrum was **not guaranteed to carry the replacement**. Indeed the nominal D0 and Dfloor full-Planck runs had exactly the same full-Planck profile chi2=2789.7048259470976, an obvious concern. The historical results from those branches must not be promoted into independently checked D-coordinate likelihood evidence.

## Corrected approach

- Use the successfully re-evolved accepted linear-G3 covariant action reconstruction \`sdmc/run_late300_linear_covariant.py\`, where **AF, zc, width, D0, power and Dfloor are parsed from the actual \`output/linear_cov_target.ini\`**.
- Preserve separate Planck-mass and kinetic windows exactly as in \`sdmc/patch_F_only_force_center.py\`: F uses zc=4.034077502899133, width=0.32925377073762596 and AF=0.02052; the D0 correction retains its historical zc,D=3.927876388467848, width,D=0.33114133956842123, pD=1.
- Rebuild the complete action **separately for every candidate**; \`sdmc/verify_late300_D_replacement.py\` proves that the target native D(N), reconstructed action D0/Dfloor, and free-evolution trajectory each match the candidate, and verifies healthy scalar propagation. Distinct action C_l SHA256 values verify the cases did not silently replay the same spectra.
- Results are **conditional model comparisons**, not independent derivations of \`D0\` or \`Dfloor\`. Accepted late300 is the **kinetic-sector reference on the derived-F background**, not the historical full late300 profile.

## Candidates

| Label | Dfloor baseline | D0 correction | Formula / provenance |
|---|---:|---:|---|
| accepted_D | 0.05181542627513409 | 0.34231919445927034 | Accepted data-selected baseline and correction |
| D0_2pi_lambda | 0.05181542627513409 | 0.34136151074659893 | Candidate \(D0=2\pi/\lambda_e\), not independently proven |
| D0_Xi_sqrtF | 0.05181542627513409 | 0.342002216921491 | Candidate \(\Xi_v^{-1}/\sqrt{F_0}\), not independently proven |
| D0_force_candidate | 0.05181542627513409 | 0.3419932122356156 | Previously proposed force-normalized candidate |
| floor_16_lambda2 | 0.047226890271848634 | 0.34231919445927034 | Canonical radiation tracker floor \(16/\lambda_e^2\) |
| combined_2pi_floor | 0.047226890271848634 | 0.34136151074659893 | Both canonical floor and \(2\pi/\lambda_e\) correction |

All six passed the covariant action target, spline/free-evolution and high-l Planck tests. Free H(z) errors over z<=100 were below 8e-10. Run **37806771106**, six separate archived artifacts.

## Planck 2018 high-l TTTEEE Plik-lite with calibration prior

| Case | chi2 | Delta chi2 vs accepted |
|---|---:|---:|
| accepted_D | 591.838788 | 0 |
| D0_2pi_lambda | 591.834702 | -0.004086 |
| D0_Xi_sqrtF | 591.832565 | -0.006223 |
| D0_force_candidate | 591.838436 | -0.000352 |
| floor_16_lambda2 | 591.575133 | -0.263655 |
| combined_2pi_floor | 591.571978 | -0.266810 |

High-l evidence alone cannot rank global SDMC vs local021.

## Full Planck 2018 (TTTEEE, low-l TT, low-l EE, lensing, nuisance priors)

Run **37807119146** completed successfully; every nuisance Minuit fit reports valid. Same fixed cosmological parameters and 21 nuisance freedom across candidates.

| Case | Full Planck profile chi2 | Delta chi2 vs accepted |
|---|---:|---:|
| accepted_D | 2784.82549538608 | 0 |
| D0_2pi_lambda | 2784.821140904252 | -0.004354481828 |
| D0_Xi_sqrtF | 2784.8195024676133 | -0.005992918467 |
| D0_force_candidate | 2784.8246840476986 | -0.000811338381 |
| floor_16_lambda2 | 2784.5593024304303 | -0.266192955650 |
| combined_2pi_floor | 2784.5561712625763 | -0.269324123504 |

The **canonical floor** accounts for virtually all improvement; alternative D0 prescriptions have minuscule impact. Differences are too small to establish statistical superiority or a first-principles formula. Not a posterior-marginalized Bayes factor. This is a *within-SDMC conditional* comparison; independent local021 must be reoptimized and the correct parameter accounting applied.

## DESI DR1 full-shape P0/P2/P4

Official-covariance DESI workflow **37807287174** recompiles the same candidate-specific linear-G3 action for each candidate, with updated consistent primordial tilt. Its fixed ΛCDM comparison is **not** a newly optimized local021 control. Results confirmed as of report creation:

| Case | DESI profile chi2 | Delta vs workflow fixed LCDM control |
|---|---:|---:|
| floor_16_lambda2 | 329.52467310279377 | -11.007730289881465 |
| combined_2pi_floor | 329.5246928853224 | -11.007710507204138 |

**The reference accepted_D and D0_2pi_lambda jobs remained active in their official DESI dependency stage at this checkpoint.** Do not infer a full DESI delta relative to accepted_D until those results are obtained. A difference of 1.98e-5 in chi2 between the two canonical-floor variants is operationally negligible at present.

## Physical interpretation / parameter provenance

D0 is **a transition-dependent correction to** Dfloor, not a separate, unrelated kinetic baseline. D(N)=Dfloor+D0*S_D(N)^pD. The possible canonical-floor relation, and the candidate D0 formulae, are still conditional on externally selected or reconstructed structural quantities. The earlier reconstructed action was designer-matched to selected D(N), and the action-preserving perturbation construction \`deltaG2=mu(phi)*(X-Xbg(phi))^2\` proved this background does not uniquely select the kinetic correction. AIC/BIC/Bayesian parameter penalties cannot be reduced without a genuinely independent SDMC action predicting the correction and its kinetic window.

## Links

- True-D reconstruction and Planck Plik-lite: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37806771106
- True-D full Planck: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807119146
- True-D official DESI full-shape: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37807287174
- Independent correction & action non-uniqueness: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37803894563
- Refined causality bound \`D0_min≈0.01918832585\`: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37804665278

## Outstanding before a model-selection verdict

Complete remaining DESI profiles, require a fair joint Planck+DESI+SN+HSC inference with both models free to optimize and consistent likelihoods, recompute AIC/BIC and **fully marginalized** family evidence accounting for historically selected D0, Dfloor and D window. Separately derive the kinetic coefficient functions and early/late endpoints from action and boundary conditions, rather than retrofitting them to the chosen late300 history.
