# SDMC early massive galaxies and SMBH: microscopic-coefficient closure audit (10 October 2026)

**Status: no unique microscopic coefficient pair derived; no new candidate promoted.**  
**Baseline:** frozen `sdmc/late300_candidate.json`; accepted Jordan/No-Slip IR action unchanged.  
**Reproducible workflow:** [GitHub Actions 38049191916](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38049191916), success; artifact `sdmc-early-galaxy-smbh-uv-consistency`.  
**Code:** `sdmc/audit_early_galaxy_smbh_microphysics.py`.  
**Manuscript basis:** Manuscript B, spatial-map §§50–58, nonlinear Jacobian §§25–36, accepted-action/JWST §§58–59 and verdict §72. Earlier conditional JWST workflow: [37581925341](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37581925341). Independent observational baseline: Comini, Vagnozzi & Loeb, arXiv:2604.13866v2, §4.1, and Xiao et al., Nature 2024.

## A. What the accepted action actually fixes

The **accepted linear Horndeski/Jordan action** determines the fitted late300 background and linear modes. The **candidate microscopic spatial-map sector** is extra ultraviolet physics; Manuscript B §54 explicitly states that the infrared Horndeski action does not determine it uniquely.

Let `h_ij` be the Jordan spatial metric and `gamma_ij=M^2 exp(2u) h_ij` the isotropic spatial structural metric. Define `u=delta ln M_loc=ell/3` (ell is the spatial trace). The candidate quadratic local potential reduces **exactly within this sector** to
`V_map^iso=(mu_V^4/2)u^2`.
This potential has a unique minimum `u=0`: without an additional source, it cannot create local matter-era mapper displacement.

The proposed nonlinear operator is
`Delta V_nl = mu_V^4 g3 W_m(p) h_loc kappa^3 u`,
with `W_m(p)=4(p-1)(2-p)`, `kappa=(3)R[h]/(6H^2)`. A conditional quasistatic minimum yields
`u=-g3 W_m h_loc kappa^3`, `Q_M=-u/delta` and `Q_M~(125/27)g3 W_m delta^2` for small delta. That is why the operator is negligible at **linear order**.

**Non-identifiability of g3:** this `u*kappa^3` term is higher than quadratic order around flat FLRW, so its value cannot be inferred from matching homogeneous H(z), linear TT/TE/EE, BAO, linear power or No-Slip. The pre-existing spatial action contains no coefficient-fixing identity tying `g3` to the late `A_F` in `F(N)`. Varying `g3` leaves the above accepted linear data unchanged at the displayed perturbative order. Therefore `g3=A_F` is a *numerical hypothesis*, not a derivation.

## B. New geometric reduction of the curvature selector: alphaK_spatial

Manuscript B §§54–56 adds `(W eta_K/2) ^(3)R_ijkl[gamma] ^(3)R^ijkl[gamma]` to the positive spatial potential, with `eta_K>0`.

For a single isotropic Fourier trace perturbation on flat `h_ij=a^2 delta_ij` and locally homogeneous `M`,

`^(3)R_ij[gamma]= (k_i k_j+delta_ij k^2) u` (indices lowered, linear order).

Consequently, `Ricci_ij Ricci^ij=6 (k/a)^4 u^2/M^4`, `R^2=16 (k/a)^4u^2/M^4`; the 3D identity `Riem^2=4 Ricci^2-R^2` gives

`Riem[gamma]^2 = 8 (k/a)^4 u^2 / M^4` (quadratic trace action).

The trace kinetic convention of manuscript §56.2 is `(B f_ell^2/2) dot(ell)^2 = (9B f_ell^2/2)dot(u)^2`. With the same outer locking `B` multiplying the selector term (§54.7), matching the fourth-derivative contribution to `(9B f_ell^2/2) alphaK_spatial (k/a)^4 u^2` gives

`**alphaK_spatial = 8 W eta_K / (9 f_ell^2 M^4).**`

This derivation fixes the **relative geometric numerical coefficient 8** for the pure isotropic conformal trace and this normalization. It **does not** fix `eta_K/f_ell^2`, the time-dependent `W`, or numerical `alphaK_spatial`; mixed/anisotropic and constraint sectors require a full off-shell Hamiltonian derivation. Do **not** confuse this coefficient with Horndeski/Bellini-Sawicki `alpha_K`, a different variable.

**Critical timing problem:** Manuscript B §55.2 explicitly requires the prethermal spatial curvature selector `W -> 0` in the infrared after selecting the flat branch. Therefore the inherited selector has `alphaK_spatial -> 0` by the galaxy-formation era. Using it as a nonzero dwarf-scale filter at z~5–20 requires either a **new matter-era reactivation of W**, or a separate allowed spatial-gradient invariant with an independent coefficient. Neither law is present as a parameter-free prediction. The earlier claimed availability of a k^4 galaxy filter is therefore *architectural*, not achieved for the frozen original selector history.

## C. Independent quasistatic halo screening test

Introduce the dimensionful trace mass `m`, propagation speed `c_ell`, and the spatial quartic coefficient `alphaK_spatial`. The subhorizon static response of the trace to a specified local nonlinear source has the conditional normalized form

`F(k,z) = m^2/[m^2+c_ell^2(k/a)^2+alphaK_spatial (k/a)^4]`.

Use units c=1 for inverse lengths. At z=5.5, applying the frozen late300 H(z) prescription (`H0=69.71482083`, `Omega_m0=0.2985`) yields
`H(z)=628.34065 km/s/Mpc`, `H(z)/c=0.00209591882 Mpc^-1`, `a=1/6.5`.

The manuscript's `m/H=300` is an **illustrative heavy decoupling hierarchy**, not a measured or action-derived mass; `m/H=5–10` was earlier considered as a *softened halo mode*. The following direct screening audit makes the missing gradient normalization important.

| assumed m/H | comoving halo k (Mpc^-1) | F(k) if c_ell=1, alphaK=0 | c_ell needed for F >= 0.9, alphaK=0 |
|---:|---:|---:|---:|
| 10 | 0.553 (massive galaxy) | 0.000033998 | <0.001944 |
| 10 | 1.191 (galaxy halo) | 0.000007330 | <0.000902 |
| 300 | 0.553 | 0.029691 | <0.05831 |
| 300 | 1.191 | 0.006554 | <0.02707 |

The inequality is `c_ell <= m a/(3k)` if F >= 0.9 when alphaK=0. Positive quartic gradients make F *smaller*. Therefore a gradient term with `c_ell~1` would suppress the proposed halo response by 30–100000-fold in the examples; the previous local algebraic minimum effectively assumed gradients could be neglected. To keep the galaxy-scale mechanism alive, a **very slow trace propagation mode, a much higher trace mass, or another derived nonlocal source structure** would be needed. None is fixed by late300.

The *conditional* diagnostic high-k cutoff `k_c=3 Mpc^-1` from Manuscript B can be translated into a coefficient if `c_ell=0`, `W>0` and one sets `F(k_c)=1/2` at a **specified** epoch:

`alphaK_spatial = m^2 a^4/k_c^4`.

At z=5.5 this gives `3.03816e-9 Mpc^2` for m/H=10 or `2.73434e-6 Mpc^2` for m/H=300. These **are not two derived predictions**: both satisfy the *same arbitrarily selected* k_c when m/H differs. For nonzero c_ell, `alphaK=(m^2-c_ell^2 (k_c/a)^2)/(k_c/a)^4` can even become **negative**, violating the desired stabilizing sign. At m/H=10, positive alphaK with this k_c requires `c_ell<0.001075`; for m/H=300 it requires `c_ell<0.03225`.

The inverse problem remains underdetermined even given a preferred cutoff: unknowns include `eta_K,f_ell,M,W,m,c_ell` and the microphysical source. Stability inequalities `f_ell^2>0, c_ell^2>0, m^2>0, eta_K>0` restrict *signs*, not values.

## D. The only quantitative g3 values currently available are *target fits*

The completed [late300 JWST likelihood-equivalent audit](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37581925341) gives:
- **Conditional diagnostic:** `g3 = A_F = 0.0204814649`, `Qeff=0.07028`, `Aeff=1.21226`, inferred FRESCO `epsilon95~0.315`.
- **Target for roughly epsilon95=0.20:** `g3~0.0348160`, `Qeff~0.11947`.
- **Target for roughly epsilon68=0.20:** `g3~0.0442275`, `Qeff~0.15177`.

These numbers assume the effective local source remains unscreened at relevant halo scales. If the conditional `F(k,z)` is small, holding the same effective response requires `g3_required ~ g3_target/F`. At k=1.191 and c_ell=1, alphaK=0 this implies `g3_required~4.75e3` for the softened m/H=10 case or `~5.31` even for m/H=300. Such replacements are **not predictions or validated perturbative couplings**.

Furthermore, the manuscript's proposed modified top-hat source `delta_S=delta_J-3u` is still a **phenomenological local-volume gravitational closure**. It has not been obtained by varying a fully coupled nonlinear Jordan+spatial-triad+matter action. The EdS top-hat curvature identity must itself be re-derived when structural backreaction is included.

The real [Comini–Vagnozzi–Loeb 2026](https://arxiv.org/html/2604.13866v2) flat-LambdaCDM FRESCO 95% lower efficiency **0.63** is correctly reproduced in the manuscript's reference. However, the manuscript's *transferred SDMC limits* are calibrated interpolations, not independently rerun FRESCO posteriors; the original authors emphasize observational uncertainties and the importance of astrophysical efficiency assumptions.

## E. Independent SMBH age and ideal growth envelope

Under the *late300 background* and adopting the **earlier manuscript's conditional** collapse shift z_seed 20 -> 23.26, the proper ages are `t20=175.6714 Myr`, `t23.26=141.2866 Myr`, `t7=753.2615 Myr`. The illustrative gain is `34.3848 Myr`, not additional structural-clock accretion time.

For ideal Eddington exponential accretion `M=Mseed exp(duty*Delta t/45 Myr)`:
- 100% duty from z23.26 to z7 gives growth `8.0568e5`, needing an initial seed of at least `1.241e3 Msun` to reach `1e9 Msun`.
- 50% duty gives growth `897.6`, needing a seed of at least `1.114e6 Msun`.
- At the same 50% duty, z20 seed formation requires `1.632e6 Msun`.

These are **idealized mathematical requirements**, not a solved SMBH formation channel. No independent SDMC black-hole seed mass function, accretion rate/duty distribution, radiative-efficiency evolution, feedback or host-merger population has been derived, and no high-z SMBH likelihood was run.

## F. Scientific verdict and next non-circular requirement

1. **Derived:** the conformal-trace geometry of the spatial curvature selector yields a symbolic `alphaK_spatial` relation (with convention-dependent normalization), and the static screening and age bounds above.
2. **Not derived:** numerical `g3`, numerical `alphaK_spatial`, a nonzero halo-era selector `W_halo`, its compatible `c_ell`, an exact joint nonlinear Einstein/trace/halo force, and SMBH seed/accretion microphysics.
3. **Conditional retained target:** g3 of a few hundredths may improve rare-halo abundance *if* the nearly algebraic local mapper response survives the gradient and full constraints. This is not presently demonstrated.
4. **Do not promote** any new candidate or modify late300 based on this audit. A fully defined microscopic UV action (couplings, scales, activation, kinetic and curvature terms) must first supply numerical parameters *without* inspecting JWST; then solve coupled nonlinear perturbations, compute mass variance and halo mass function, build galaxy selection/SED/volume likelihood, and integrate SMBH population with seeds and duty cycle. Independently recover full low-z/CMB likelihoods.

**Reproducibility:** this branch adds the Python audit and isolated GitHub workflow and leaves the baseline candidate untouched. Automated success verifies the described mathematical checks, not the existence of an observational solution.
