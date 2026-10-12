# SDMC versus ΛCDM: standalone supernova tests FIRST; joint CMB + DESI BAO tests SECOND

**Date:** 12 October 2026  
**Source run:** [50 fully evolved candidates, matched CLASS+BOLTZMANN screen](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38274914841)  
**Two-stage postanalysis:** [38277644360 — SUCCESS](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38277644360)  
**Actual official DESI DR1 P0/P2/P4 full-shape:** [38277970929 — distinct extended workflow](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38277970929), status **SUCCESS for both matrix cases**, verified against separate job logs and official DESI FS artifact CSVs; do NOT copy BAO chi2 into full-shape column.

## Scientific correction

The prior 50-candidate run used **DESI DR1 Gaussian BAO distances only**. Neither the full-shape power-spectrum data \(P_0,P_2,P_4\), nor nuisance-parameter marginalisation/profiling thereof, were run. This is now explicitly separated from the official DESI FS run.

Standalone uncalibrated SN datasets have overlap and different nuisance handling and must be **separate**:
* Pantheon+ covariance likelihood (uncalibrated SN, profiled intercept).
* Union3 covariance likelihood (uncalibrated SN, profiled intercept).
* DES-Y5 covariance likelihood (uncalibrated SN, profiled intercept).

These are independent **computational tests**, NOT statistically independent overlapping astronomical samples suitable for summing their chi2 as one combined dataset.

## STAGE 1: Supernova-only candidate search (NO Planck, NO BAO)

All 50 source likelihood values were recomputed **only as rankings**, without rerunning or re-fitting the 50 spectra. The only metric for selection at this stage is the named SN's \(\chi^2\). Across the explored 25-member population in each model family, the standalone SN best is **ordinary_sobol02** for both families.

| Standalone SN likelihood | ΛCDM best χ² | SDMC best χ² | Δχ² SDMC−ΛCDM | Best ID |
|---|---:|---:|---:|---|
| Pantheon+ | **1404.327979** | **1404.215307** | **−0.112672** | ordinary_sobol02 both |
| Union3 | **27.030829** | **26.965104** | **−0.065725** | ordinary_sobol02 both |
| DES-Y5 | **1646.574803** | **1646.453556** | **−0.121246** | ordinary_sobol02 both |

Best ΛCDM SN-only H0 input = **67.5425601** km/s/Mpc; best SDMC SN-only H0 input = **68.6887935** km/s/Mpc. **Important:** This is the preassigned trial parameter, *not* an SN determination of H0, because the intercept was profiled. Report age/rdrag/H0 only after an external calibration or absolute-distance information.

With these SN-only winners, Planck and BAO scores are much worse than their baselines. For example, the standalone Pantheon+ winner's Planck-lite χ² is **1087.824 (ΛCDM)** / **1097.857 (SDMC)**, versus model baselines **1018.698** / **1020.954**. SN-only optimization is not a valid combined cosmology score by itself.

## STAGE 2: Planck + DESI DR1 BAO + ONE SN likelihood, separately

This second stage only ran after stage 1 finished and archived. The exact same 50 independently evolved model candidates are ranked anew for each exclusive joint score.

| Joint objective | ΛCDM best χ² | SDMC best χ² | Δχ² SDMC−ΛCDM | Best IDs |
|---|---:|---:|---:|---|
| Planck-lite + BAO + Pantheon+ | **2438.306962** | **2440.125875** | **+1.818913** | local021, SDMC zc_minus |
| Planck-lite + BAO + Union3 | **1060.875472** | **1062.790339** | **+1.914867** | local021, SDMC zc_minus |
| Planck-lite + BAO + DES-Y5 | **2682.728364** | **2684.566645** | **+1.838281** | local021, SDMC zc_minus |

The SDMC joint winner is a **frozen-covariant replay verified** \(z_c=3.747876388...\) candidate, slightly improved over accepted Late300 \(z_c=3.927876...\) under these **screening** approximations. The original accepted Late300 remains unchanged.

No additional complexity or "extra parameter" penalty is applied to either raw χ². No novel full-Plik, DESI full-shape or calibrated SN likelihood is implied by the stage-1/stage-2 data.

## Actual DESI DR1 full-shape extension (completely separate)

A separate two-case covariant/Horndeski validation workflow compares:
1. local021 ΛCDM full-shape transfer predictions,
2. **accepted Late300** frozen linear-G3 action,
3. **candidate \(z_c=3.747876\)** frozen linear-G3 action.

It uses the existing official DESI DR1 six-tracer covariance/window files, multipole \(P_0/P_2/P_4\), and the repository's established EPT model plus multistart nuisance profiling. Per-case DESI results must be labeled **official DR1 full-shape data analysed using a custom nuisance/PT implementation**, *not* an official DESI collaboration posterior.

It is a long, separate workflow and its own success/failure, diagnostics, final χ² and possible comparison to prior full-Plik/FS work must be checked from **its new run artifacts**, never inferred from BAO-only tables. Do not replace local021 or Late300 based on this screening result.

## Source and artifacts

Two-stage tests: \`sdmc/run_separate_sn_then_joint_ranking.py\`; \`sdmc/test_separate_sn_then_joint_ranking.py\`; \`.github/workflows/sdmc_sn_first_then_joint.yml\`.

Official DESI full shape extension: \`.github/workflows/sdmc_matched_three_candidates_fullshape_desi.yml\`, adapted from the existing SDMC six-bin official DESI FS profile test and the verified accepted Late300 linear-G3 reconstructed action.

**Verdict:** SN-only samples each weakly prefer the best SDMC screen by raw chi2 of order 0.07–0.12; joint Planck-lite + DR1 BAO + one SN favors local021 by ≈1.82–1.91. These are different questions, and actual DESI full-shape must be evaluated separately.


## FINAL FOLLOW-UP: official DESI DR1 galaxy full-shape results (actual jobs completed)

[Full-shape GitHub Actions workflow #38277970929](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38277970929) successfully ran both matrices, **accepted_Late300** and **candidate_zc_minus**, using the restricted accepted linear-G3 covariant action reconstructed once per candidate and then evolved through the CLASS transfer solver with expansion-model input removed. The same locally optimized flat-ΛCDM **local021** was matched in both cases. DESI six-tracer official DR1 input HDF5/window/covariance files, P0/P2/P4, the existing EPT/REPT transfer treatment and multistart nuisance profiling were used.

### Standalone DESI DR1 full-shape ONLY

| Fixed model | Custom profiled official-DR1-data FS χ² | Δχ² versus local021 |
|---|---:|---:|
| ΛCDM local021 | 331.937262326103 | 0 |
| **Accepted SDMC Late300** | **329.451426389657** | **−2.485835936446** |
| **Modified SDMC zc-minus** | **329.456368232153** | **−2.480894093950** |

The accepted Late300 fits DESI FS better than zc-minus by χ²=**0.0049418425**, a negligible difference. Both FS calculations have six tracer blocks with optimizer-best success; some independent optimization seeds converge to slightly different local minima, so finer reproducibility/nuisance profiling is desirable. These results come from **official DR1 survey data in the repository's custom PT+profile model**, not a DESI collaboration full posterior.

Full-shape workflow artifacts:
- \`desi-fullshape-accepted_Late300-versus-local021\` (per-tracer CSV, profile summary, full CLASS transfers).
- \`desi-fullshape-candidate_zc_minus-versus-local021\`.

### Standalone and joint must not be conflated

**SN-only** results above were independently ranked before Planck or BAO were added. The **BAO-only three-joint results** above remain as the separate BAO-likelihood screen.

Because DESI DR1 full-shape contains related clustering and BAO information, **do not naively add DESI DR1 BAO χ² on top of its full-shape χ²**. A separate exploratory joint comparison instead combines Planck-lite + DESI DR1 *full shape* + one of the three exclusive SN datasets:

| One-SN Planck-lite + DESI FS joint | local021 χ² | accepted Late300 χ² | modified SDMC χ² | Δ χ² modified SDMC − local021 |
|---|---:|---:|---:|---:|
| Pantheon+ | 2756.850657636 | 2756.534441762 | **2756.425512759** | **−0.425144877** |
| Union3 | 1379.419167690 | 1379.198906147 | **1379.089977143** | **−0.329190546** |
| DES-Y5 | 3001.272059515 | 3000.975212015 | **3000.866283012** | **−0.405776503** |

These are arithmetic sums of separately profiled fixed-candidate likelihood terms. They are **not** new full joint re-optimizations or Bayesian evidences; SN overlap means never sum all three joint columns. They differ from BAO-only screens, which favored local021. The sign change demonstrates why galaxy clustering full shape is an important discriminator and why the CMB alone does not choose a universally preferred H0.

Do not replace accepted Late300 with the zc-minus candidate on χ² differences ~0.1 or claim a robust Hubble-tension resolution. Full official nuisance-prior likelihoods, nonlinear validity and family-level parameter counting remain essential before scientific promotion.
