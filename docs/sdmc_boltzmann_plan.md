# SDMC Boltzmann implementation plan

This branch is reserved for direct Einstein-Boltzmann tests of the SDMC observational candidates.  The upstream `master` branch is kept untouched as the hi_class reference baseline.

## Models to run

The comparison suite is intentionally split into four runs:

1. **LCDM control** — unmodified hi_class/CLASS reference.
2. **NKp-v2 B0** — SDMC background + tracker + endpoint-clean No-Slip gravity, with standard photon-baryon acoustics.
3. **Kp core** — frozen Kp background + tracker + frozen No-Slip gravity, but with the terminal reciprocal-acoustic operator disabled.
4. **Kp full** — the same frozen Kp core with the terminal reciprocal-acoustic operator enabled in the photon-baryon hierarchy.

This split is required so that the CMB cost of the Kp terminal operator can be isolated from the common SDMC early/late sectors.

## Frozen background registry

### Common values

- `H0 = 70.8514 km/s/Mpc`
- `omega_b = 0.02239952`
- `Omega_r = 8.3268e-5`
- early tracker slope `lambda_e = 20`
- early structural sound speed `c_phi,e^2 = 0.003`
- default tracker handoff width `Delta N_t = 0.5`

### Kp branch

- `Omega_m0 = 0.2836698`
- `r_s(z*) = 143.779 Mpc`
- `r_d = 145.2577 Mpc`
- late correction

  `deltaH(z) = 0.874 * [0.024283 z exp(-z/0.25) - 0.0132767 z^2 exp(-z/1.5)]`

- tracker handoff diagnostic: `z_t ~ 21.26` for `Delta N_t = 0.5`
- frozen No-Slip targets from the manuscript: `M_*^2(0)/M_Pl^2 ~ 1.02559`, `alpha_M,max ~ 0.01447`, transition centered near `z ~ 4.38`

The exact frozen numerical Kp gravity interpolation must be reconstructed/validated against these reported targets before the final likelihood run.  An analytic proxy may be used for code smoke tests only and must never be labeled as the final frozen Kp result.

### NKp-v2 branch

- `Omega_m0 = 0.2925181427`
- `r_s(z*) = 142.688139 Mpc`
- `r_d = 145.395978 Mpc`
- late correction

  `deltaH(z) = 0.01105624999 z exp(-z/0.25) - 0.01951933685 z^2 exp(-z/1.5)`

- tracker handoff diagnostic: `z_t ~ 21.10` for `Delta N_t = 0.5`
- endpoint-clean No-Slip closure is branch-specific and must not inherit the Kp profile unchanged.

The current B0 construction uses an early GR endpoint (`F=1`, `alpha_M=0`), exact No-Slip (`alpha_B=-2 alpha_M`), `alpha_T=0`, and returns to `alpha_M(0)=0` while keeping the scalar stability numerator positive.

## No-Slip implementation

Define

`F(a) = M_*^2(a)/M_Pl^2`

and

`alpha_M = d ln F / d ln a`.

For exact No-Slip gravity use

`alpha_B = -2 alpha_M`,

`alpha_T = 0`.

The scalar stability numerator used by the B0 shooting construction is

`N_s = -2(1+alpha_M)(h+alpha_M) - 2 alpha_M' - (3 Omega_m + 4 Omega_r)/F`,

where prime denotes `d/d ln a` and `h = d ln H / d ln a`.

The Boltzmann implementation must use the code's native stability variables and checks; the equation above is the reconstruction target/diagnostic, not a request to bypass hi_class stability logic.

## Early tracker

For the early exponential tracker with slope `lambda_e`, use the standard tracking fraction as a background diagnostic,

`Omega_phi ~ 3(1+w_b)/lambda_e^2`,

with a smooth handoff near `z_t ~ 21` to the late-time SDMC reconstruction.  The full implementation must evolve the perturbations consistently with `c_phi,e^2 = 0.003`; it is not sufficient to modify only `H(z)`.

## Kp terminal acoustic operator

The retained Kp terminal equation is represented schematically by

`Theta0'' + [R'/(1+R) - C'/C] Theta0' + C^2 k^2/[3(1+R)] Theta0 = S_K`.

The operator is active only in the terminal photon-baryon acoustic window and overlaps the visibility/recombination epoch.  It must therefore be implemented in the perturbation/tight-coupling evolution, not as an a-posteriori shift of `r_s`, `r_d`, or the C_l peak positions.

The terminal profile must satisfy the frozen integral constraint associated with the Kp drag/recombination tail.  Until the recombination-resolved profile has been recalibrated in code, the terminal operator should remain behind an explicit runtime flag.

## Kp baryon-CDM relative sector

The terminal phase map produces non-standard relative transfer functions.  For the eventual DESI raw-multipole test retain and output

- `T_delta_bc(k) = T_b(k) - T_c(k)`
- `T_theta_bc(k) = T_theta_b(k) - T_theta_c(k)`

and do not replace the relative velocity by a naive baryon-velocity amplitude rescaling.  The final galaxy full-shape likelihood must include tracer-dependent relative-density/relative-velocity bias operators.

## Source-code map

The first implementation audit identified these upstream modules as the main touch points:

- `include/background.h` — add explicit SDMC model identifiers/state if needed.
- `gravity_smg/gravity_models_smg.c` — register and parse a branch-specific SDMC No-Slip parameterization.
- `gravity_smg/gravity_functions_smg.c` — return the branch-specific EFT/alpha functions.
- `gravity_smg/background_smg.c` / `source/background.c` — branch background and Planck-mass evolution.
- `source/perturbations.c` and `gravity_smg/perturbations_smg.c` — tracker perturbations and Kp photon-baryon modification.
- `source/thermodynamics.c` — visibility/tight-coupling timing diagnostics; recombination physics itself should remain standard unless a derivation explicitly changes it.
- `source/transfer.c` / `source/output.c` — expose baryon/CDM relative transfer diagnostics.
- `source/harmonic.c` / `source/lensing.c` — no ad-hoc SDMC corrections should be added here; these modules should consume the physical perturbation sources produced upstream.

## Validation ladder

Every new stage must pass the previous one before the model is interpreted physically:

1. Unmodified LCDM build and spectrum reproduction.
2. SDMC runtime flag off reproduces the same control spectra.
3. NKp-v2 background-only distances reproduce the frozen `r_s`, `r_d`, and `D_M(z*)` targets.
4. Add tracker; verify early fraction and handoff closure.
5. Add No-Slip functions; require native hi_class ghost/gradient/tensor stability checks to pass.
6. Produce unlensed and lensed TT/TE/EE plus lensing potential for NKp-v2.
7. Repeat Kp core with the terminal operator off.
8. Enable terminal Kp operator and inspect visibility-weighted TT/TE/EE residuals.
9. Only after the spectra are stable, evaluate Planck likelihoods.
10. For DESI raw full shape, include the Kp relative baryon-CDM transfer operators before quoting a final Kp likelihood.

## Scientific-status rule

Compressed CMB, distance-prior, phase-shift, ShapeFit, broad-band, and growth screens are useful diagnostics, but none is to be relabeled as a full TT/TE/EE or raw DESI multipole likelihood.  Final claims require the corresponding solver and likelihood calculation.
