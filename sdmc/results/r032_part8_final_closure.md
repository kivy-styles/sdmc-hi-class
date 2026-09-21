# SDMC r032 — Part VIII final validation / closure ledger

## Current superseding Planck + DESI closure

The previous fair-control stopping statement below is retained as a historical snapshot, but it is superseded for the Planck + DESI question by the exact-covariant **eref015** result.

### Exact-covariant SDMC eref015

- H0 = 70.79146512830793
- omega_b = 0.022083219194622913
- omega_cdm = 0.12299536722293603
- n_s = 0.9625227132590487
- tau_reio = 0.055202901571989066
- A_s = 2.1218418557926507e-9
- ln(10^10 A_s) = 3.054869604391977
- A_F = 0.03597474167915061
- z_c = 3.580509699881077
- width = 0.3986928062886
- D0 = 0.34231919445927034
- D_floor = 0.04906211170135066
- lambda_e = 18.40625
- z_t = 16.141839363323523

Exact-covariant stability:
- min D = 0.04906211167807
- min c_s^2 = 0.007939174161602
- max c_s^2 = 0.6674784881772
- max |alpha_B + 2 alpha_M| = 1.0005885009434223e-14

### Exact likelihood comparison

eref015:
- full Plik chi2 total = 2786.1184593309654
- full Plik delta vs common LCDM reference = -17.875731750041723
- DESI DR1 raw full-shape chi2 total = 327.0937584010189
- DESI delta vs common LCDM reference = -13.438644991802846
- P+D delta vs common reference = **-31.31437674184457**

LCDM local021:
- full Plik chi2 total = 2781.4370287015954
- full Plik delta vs common LCDM reference = -22.557091181603937
- DESI DR1 raw full-shape chi2 total = 331.93726232610265
- DESI delta vs common LCDM reference = -8.59514105435187
- P+D delta vs common reference = **-31.152232235955807**

Direct fair comparison:
- Planck: eref015 is +4.681430629370035 chi2 above local021
- DESI: eref015 is -4.843503925083752 chi2 below local021
- Planck + DESI: eref015 is **-0.162073295713716 chi2** relative to local021

Using each run's stored delta-to-reference values gives -0.162144505888762; the ~7e-5 difference is only the tiny independent Planck reference-profile optimizer offset. The direct absolute likelihood comparison is the clean fair number.

**Conclusion: the formerly positive Planck+DESI gap is closed and surpassed.**

Primary exact runs:
- eref015 full Plik: run 35667433363, job 106556203694
- eref015 DESI: run 35667451332, job 106556258785
- LCDM local021 full Plik: run 35585790426, job 106288695887
- LCDM local021 DESI: run 35587244801, job 106293308905

### Reconstruction audit rule

A promoted exact run is admissible only when its covariant reconstruction summary echoes the intended A_F, z_c, width, D_floor and H0 of the promoted candidate.

Runs excluded from the closure ledger because they reconstructed an older action include:
- 35666650275: stale edge028 full-Plik reconstruction
- 35666652503: stale edge028 DESI reconstruction
- 35666724761: initial localref008 full-Plik run using the prior lowridge action
- 35665667528: initial lowridge021 DESI run using D_floor=0.041875
- 35666871548 and 35666851338: stale edge028-Qbest DESI reconstructions

The exact O4D030 promotion is still queued and may improve the post-closure leader; it is not required to establish the closure above. SN likelihoods have not yet been re-promoted for eref015, so this superseding statement is specifically the fair Planck + DESI result.

---

## Superseded historical snapshot

This file records the terminal result of the r032 fair-control investigation.
All deltas below are relative to the common original fixed LCDM reference unless explicitly stated otherwise.

## Final exact-covariant SDMC leader

Model: local002 ordinary cosmology + expansion refinement

- H0 = 70.5653567390982
- omega_b = 0.022011983189284802
- omega_cdm = 0.12404331885203719
- n_s = 0.9625227132590487
- tau_reio = 0.055202901571989066
- ln(10^10 A_s) = 3.0598696043919773
- A_F = 0.06017362505197525
- z_c = 2.827464461401105
- width = 0.5514493708219379
- D0 = 0.34231919445927034
- D_floor = 0.045
- lambda_e = 17.7
- z_t = 17.1
- late-background A = 0.01105624999, tau_A = 0.25
- late-background B = 0.01951933685, tau_B = 1.5

Exact-covariant stability:
- min D = 0.04499999923422
- min c_s^2 = 0.08758419507446
- max c_s^2 = 0.7349576876011
- max |alpha_B + 2 alpha_M| = 1.0005885009434223e-14

## Exact likelihood ledger

SDMC expansion-best:
- Planck full Plik delta chi2 = -9.002764487901459
- DESI DR1 raw full-shape delta chi2 = -12.517050371798064
- Planck + DESI delta chi2 = -21.519814859699522
- Pantheon+ delta chi2 = +4.375243047717959
- Union3 delta chi2 = +4.018381844944088
- DES-Y5 delta chi2 = +9.245343506336212
- P+D+Pantheon+ = -17.144571811981564
- P+D+Union3 = -17.501433014755435
- P+D+DES-Y5 = -12.27447135336331

Fair LCDM local021:
- Planck full Plik delta chi2 = -22.557091181603937
- DESI DR1 raw full-shape delta chi2 = -8.59514105435187
- Planck + DESI delta chi2 = -31.152232235955807
- Pantheon+ delta chi2 = +2.3297537905164063
- Union3 delta chi2 = +2.236117954074871
- DES-Y5 delta chi2 = +5.168786585330963
- P+D+Pantheon+ = -28.8224784454394
- P+D+Union3 = -28.916114281880937
- P+D+DES-Y5 = -25.983445650624844

Fair SDMC-minus-LCDM gaps:
- Planck + DESI = +9.632417376256285
- P+D+Pantheon+ = +11.677906633457837
- P+D+Union3 = +11.414681267125502
- P+D+DES-Y5 = +13.708974297261534

## Where the fair Planck gap comes from

Exact full-Plik SDMC expansion-best minus LCDM local021:
- high-l TTTEEE = +9.136206038926503
- lowT = +2.96073065188034
- lowE = +0.2553836016185187
- lensing reconstruction = +0.5111577737826387
- nuisance priors = +0.6909198253025965
- total Planck difference = +13.554397891510234

DESI favors SDMC relative to local021 by:
- -3.921909317446194

Hence the net fair P+D gap is +9.632417376256285.

## Search-closure evidence

1. Ordinary-cosmology transplant toward LCDM local021 fails immediately.
A 6.25% move from local002 toward the local021 six-parameter coordinates worsened the screening Planck score by about +7.4 chi2; larger moves rapidly became catastrophic.

2. Re-optimizing the expansion coordinates did help, but only modestly.
The promoted point (lambda_e,z_t)=(17.7,17.1) improved full Planck from about -8.19 to -9.00 and kept the DESI gain near -12.52. This recovered about one chi2 point in the fair P+D comparison.

3. The widened ordinary-parameter closure hypercube did not reveal an open edge.
The expansion-best anchor remained the best Planck and best Planck+SN point in the widened 35-point scan that explicitly opened lower omega_b, higher A_s and higher n_s directions. Thus the earlier box-edge concern is closed locally.

4. D0 is a weak likelihood direction.
The prior wide D0 audit changed the screening Planck likelihood only at the few-hundredths chi2 level across D0=0.24...0.36.

5. Earlier late-background A/tau_A scans changed the joint likelihood at only order-one level in their tested basin and did not expose a large hidden rescue direction.

## CMB structural diagnostic

Expansion-best SDMC versus LCDM local021:
- ell_A(SDMC) = 301.6798271349195
- ell_A(LCDM local021) = 301.5486238021766
- delta ell_A = +0.13120333274287077
- r_s,star(SDMC) = 142.9566350658675 Mpc
- r_s,star(LCDM local021) = 144.94901560720697 Mpc
- r_d(SDMC) = 145.61046443300566 Mpc
- r_d(LCDM local021) = 147.53012524715427 Mpc

TT and EE acoustic peaks align to roughly 0--1 multipole through the measured acoustic range. Therefore the residual is not a gross acoustic-phase / peak-position failure.

The remaining spectral mismatch is coherent and sub-percent:
- typical TT/EE differences are about 0.3--1% through much of high-l,
- the SDMC lensing-potential spectrum is roughly 5--10% higher than local021 from ell~500 upward,
- lensed and unlensed spectra show similar primary-shape differences, so lensing alone is not the whole residual.

Full Plik is sufficiently precise that these small correlated differences accumulate into the remaining ~9.14 high-l chi2 deficit.

## Part VIII stopping statement

Within the exact No-Slip r032 model family tested here, the ordinary six cosmological parameters, the principal local structural parameters (A_F,z_c,width), and the expansion coordinates (lambda_e,z_t) have been locally closed around the final leader under matched Planck/DESI/SN validation.

The remaining fair preference of LCDM local021 is therefore not explained by a missed ordinary-cosmology optimum or by the previously frozen expansion coordinates. The evidence localizes the residual primarily to high-l CMB spectral/transfer structure. Closing that gap would require a new perturbation/lensing-transfer degree of freedom (or another explicit model extension) and should be treated as a new model question rather than continued tuning of the closed exact No-Slip r032 parameterization.

## Primary validation runs

- exact-covariant expansion-best full Plik: 35593253133
- exact-covariant expansion-best DESI: 35593282108
- exact-covariant expansion-best SN: 35593665183
- widened ordinary-edge closure hypercube: 35592519217
- CMB structural diagnostic: 35593964504
- LCDM local021 full Plik: 35585790426
- LCDM local021 raw DESI: 35587244801
