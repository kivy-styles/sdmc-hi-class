# CMB inference of H0 and equal-budget ΛCDM local021 versus SDMC Late300 candidate evolution

**Date:** 2026-10-12  
**Verified completed GitHub run:** [38274914841](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38274914841) (**SUCCESS**)  
**Branch:** \`sdmc-lcdm-matched-free-evolution-search-20261012\`  
**Scope:** A matching exploratory search; **not** a full Plik nuisance optimization, DESI full-shape fit, cosmic parameter posterior or Bayes factor. Historical local021 and Late300 remain unchanged.

## 1. How H0 follows from the CMB

CMB maps determine the acoustic angular scale
\[
\theta_\ast=\frac{r_s(z_\ast)}{D_M(z_\ast)},
\quad r_s(z_\ast)=\int_{z_\ast}^\infty c_s(z)H(z)^{-1}\,dz,
\quad D_M(z_\ast)=c\int_0^{z_\ast} H(z)^{-1}\,dz
\]
in spatially flat models, alongside peak relative heights, baryon loading, equality, diffusion damping and lensing. The model's combined fit is what determines \`H0\`; acoustic angle alone cannot separate all degeneracies.

Planck 2018 finds \`100 theta_star = 1.0411 ± 0.0003\` and \`H0 ≈ 67.4 ± 0.5 km/s/Mpc\` *under flat six-parameter ΛCDM*. A useful approximate ΛCDM-only relation is \`Omega_m h^3 ≈ 0.09633\`, so \`H0 ≈ 100*0.09633/(omega_b+omega_c)\` once densities are independently constrained. Neither relation can be transplanted directly as an SDMC \`H0\` determination.

Our actual paired Boltzmann outputs:

| Observable | optimized local021 | accepted Late300 (unchanged) |
|---|---:|---:|
| H0 (km/s/Mpc) | 68.568587447 | 69.714820831 |
| z_star | 1088.78 | 1089.2775 |
| z_drag | 1059.8408 | 1059.5880 |
| r_s(z_star), Mpc | 144.949016 | 143.253971 |
| r_drag, Mpc | 147.530125 | 145.865180 |
| D_M(z_star), Mpc | 13913.062894 | 13753.251524 |
| 100 theta_star | 1.041819596 | 1.041600752 |
| ell_A | 301.548624 | 301.611980 |

These are derived from each model's CLASS background and thermodynamics and used alongside full predicted temperature/polarization/lensing *in the screening likelihood*. Do not confuse r_drag with r_s(z_star).

## 2. Controlled paired candidate experiment

- Exactly **25 evaluations for local021-family LCDM** and **25 evaluations for Late300-family SDMC** using their own correct background and perturbation physics in CLASS; same Planck-lite lowT lowE lensing+profiled calibration, same DESI DR1 Gaussian BAO likelihood, exact-covariance Pantheon+, and separately evaluated Union3 and DESY5.
- Candidate budgets: one unmodified baseline, sixteen Sobol variations of the six shared ordinary parameters with the same dimensionless offset design, plus eight axis probes per family (LCDM ordinary axes, SDMC \`A_F, z_c, width, D_floor\`). Fixed search budget is NOT equal prior volume or equivalent cosmological dimension.
- Primary screening objective = Planck-lite + DESI DR1 BAO + Pantheon+, **raw χ²**. Union3/DESY5 use separate alternative objective totals; none is double-counted alongside Pantheon+. SH0ES proxy reported separately.
- SDMC candidates must have \`D>0\` and \`0<c_s²<=1\` throughout the accepted background grid. Retain No-Slip residual diagnostics. Four candidate points failed in the scalar gradient-stability branch.

**This is a candidate-ranking screen**. Plik-lite profile is not the official full high-ℓ nuisance treatment; DESI BAO is not the DR1 full-shape P0/P2/P4 likelihood. AIC/BIC/Bayesian evidence and local distance-ladder likelihood have **not** been newly computed for the searched points.

## 3. Exact matched-screen output

| χ² component | local021 ΛCDM baseline/best | SDMC Late300 baseline | improved SDMC \`zc_minus\` |
|---|---:|---:|---:|
| Planck-lite+low-ℓ+lensing | 1018.698453918 | 1020.954475161 | **1020.840604316** |
| DESI DR1 BAO | 13.393566872 | 13.156730316 | **13.156730316** |
| Pantheon+ | 1406.214941392 | 1406.128540211 | **1406.128540211** |
| **Primary sum** | **2438.306962182** | **2440.239745688** | **2440.125874842** |
| Union3 (alternative SN) | 28.783451446 | 28.793004596 | 28.793004596 |
| DESY5 (alternative SN) | 1650.636343271 | 1650.569310464 | 1650.569310464 |

The new SDMC winner differs from accepted Late300 only by the transition centre \`z_c: 3.927876388467848 → 3.747876388467848\`. Thus its screening improvement over Late300 is **−0.113870846** raw χ² (Planck contribution only).

Against local021:
- Planck + BAO + Pantheon+: **Δχ²(SDMC−LCDM) = +1.8189126604**.
- Planck + BAO + Union3: **+1.9148669915**.
- Planck + BAO + DESY5: **+1.8382810347**.

On all three alternative screens the evaluated local021 baseline remained the winner over the sampled ΛCDM variants; the same \`zc_minus\` candidate was the highest-scoring stable SDMC point. **No overall screened SDMC victory** was obtained. These numbers do not negate the prior official full-Plik/DESI full-shape comparison; they use different likelihood approximations.

Evaluated points: **50** total. Passed ΛCDM: **25/25**; stable SDMC: **21/25**. Four SDMC gradient-instability failures occurred for \`A_F−0.0013\`, \`z_c+0.18\`, and \`width±0.027\` at their other baseline coordinates; their reported minimum \`c_s²\` was respectively −0.0128108, −0.00332303, −0.00461002, −0.00216486. These are source-labeled local perturbations, not general parameter exclusion bounds.

The uncalibrated Gaussian \`SH0ES_H0=73.49±1.048\` proxy gives:
- LCDM local021: χ²≈22.0525;
- SDMC candidate: χ²≈12.9764.
This explicitly indicates lower local-H0 tension for SDMC, but it is NOT a re-fit to calibrated SH0ES likelihood.

## 4. Independent field-equation validation of screened SDMC point

After the screening, the winner was used only once to reconstruct and **freeze** its restricted covariant Horndeski functions:
\[
G_2=k_1(\phi)X+k_2(\phi)X^2-V(\phi),
\quad G_3=g(\phi)X,\quad G_4=F(\phi)/2.
\]
The original parameterized \`expansion_model\` and fitted kinetic-profile input were then removed from the free-run INI. Field equations and perturbations were evolved dynamically under the reconstructed frozen action.

**Verified covariant replay status: PASS.**
- Maximum relative H trajectory error for z≤100: **3.5271962×10⁻⁹**.
- Maximum relative F error: **1.4796136×10⁻¹⁰**.
- Maximum relative D error: **7.3318368×10⁻⁸**.
- Maximum peak-normalized lensed Cl differences: TT **3.08045×10⁻⁶**, TE **5.36941×10⁻⁶**, EE **5.06870×10⁻⁶**, PP **8.59803×10⁻⁷**.
- Scalar stability: min D **0.0518154262**, min c_s² **0.0102068567**, max c_s² **0.6286484283**; max |alpha_B+2alpha_M| **3.17944×10⁻¹⁰**.
- Covariant replay Planck-lite χ² **1020.84295475**, BAO **13.15673342**, Pantheon+ **1406.12854173**; sum **2440.12822991**.
- Covariant versus target screening Δχ² = **+0.00235507**; covariant versus local021 Δχ² = **+1.82126773**.

Free evolution is **genuine**, but the action is still **inverse reconstructed from the screened target**, not obtained independently from UV/Planck structural boundary conditions. The test therefore establishes a dynamically admissible candidate rather than zero upstream model-family freedom.

## 5. Verdict and next decision gate

1. **No new ΛCDM model** within this sampled ordinary 6D envelope outperformed optimized local021 under the joint screen; this gives additional local support for its status as a strong baseline, not a global-optimum theorem.
2. SDMC has a **new, slightly improved, scalar-stable, dynamically replay-validated** candidate \`zc_minus_3\`; −0.114 χ² is very small. Do **not** relabel or replace accepted Late300; further full-precision tests are required.
3. The matched screening favors local021 by ∼1.82 χ². This is not the same as the previously reported full-Plik/full-shape Δχ²≈−0.2025 for the accepted Late300, which must be retained as a separate pipeline/result.
4. **Next necessary verification before promotion:** candidate-specific official full high-ℓ Planck nuisance profiling, official DESI DR1 full-shape, calibrated SN/SH0ES, weak lensing and growth. Then test whether an action-derived structural theory predicts its changed \`z_c\` instead of fitting it.
5. No parameter-complexity penalty was applied to either family's **raw observational χ²**. AIC/BIC or evidence requires independent family-level priors; freezing either action does not erase upstream fitting.

## Reproducibility

[Workflow and all steps](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38274914841). Artifact: **paired-sdmc-local021-cmb-h0-candidate-search** containing all 50 candidate records and free-action validation JSON. Source:
- \`sdmc/run_paired_sdmc_lcdm_search.py\`
- \`sdmc/run_paired_covariant_promotion.py\`
- \`sdmc/test_paired_sdmc_lcdm_search.py\`
- \`.github/workflows/sdmc_lcdm_matched_free_search.yml\`

No new matter creation and no changes to accepted Late300 or local021 candidate records.
