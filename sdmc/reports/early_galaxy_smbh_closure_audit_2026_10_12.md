# SDMC — Early Massive Galaxies and SMBHs
## Microscopic closure audit, 3-D Lagrangian diagnostics, JWST reproducibility, and seed-growth calculations
**12 October 2026 · Research continuation · No changes to canonical Late300**

### Executive result

The previously reported nonlinear mapper and JWST results survive as **conditional model diagnostics**, not independent microscopic predictions. This follow-up provides (i) an action-level **non-identifiability proof** for the local curvature coefficient `g3`, (ii) a numerically verified three-dimensional anisotropic Jacobian test (40³ Lagrangian cells, plus independent grid/seed checks), (iii) a **new quantitative screening-versus-SMBH-seed-host tension**, and (iv) an SMBH seed-growth matrix using Late300's own expansion rate. The exact 3-D covariant N-body simulation and JWST catalog posterior **cannot be legitimately marked completed**: no autonomous, uniquely normalized nonlinear structural-map/particle-matter field equations nor the complete public advanced observational selection likelihood are available in the accepted implementation.

No new matter is introduced: Lagrangian volume concentration is the only assumed mass source, as requested. The historical Late300 parameters, likelihood baseline, and observational success ledger remain unchanged.

**Evidence register**

| Investigation | Status after this run |
| --- | --- |
| Algebraic microscopic coefficient audit | **Completed:** `g3` not determined by linear/background data; UV `alpha_K` numerical value not fixed by sign/stability. |
| Fully covariant 3-D gravitational SDMC simulation | **Not closed:** needs off-shell action and physical nonlinear Poisson/constraint dynamics. |
| 3-D anisotropic Lagrangian dust diagnostic | **Computed**, 64,000 cells and six robustness checks; not an N-body prediction. |
| Mass/redshift screening compatibility | **Computed:** identifies new selectivity and SMBH tension at `1e8–1e9 M_sun`. |
| Exact FRESCO/CEERS catalog posterior | **Not performed**, because fully specified selection/covariance and SDMC nonlinear halo mass function are unavailable. |
| Published JWST likelihood-equivalent transfer | Previously performed; independently audited, **not relabeled exact likelihood**. |
| Black-hole seed mass and growth-time bounds | **Computed** under explicit Eddington-accretion assumptions, not a seed-formation solution. |

### 1. Provenance and conventions

Canonical input: [`late300_candidate.json`](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/late300_candidate.json). Background H0 = 69.71482083084993 km/s/Mpc, omega_b = 0.022083219194622913, omega_cdm = 0.12299536722293603, lambda_e = 18.40625, z_t = 16.19317622240633, A_late = 0.024822032167576252, B_late = 0.01264704344701022. We evaluate the published phenomenological `H(a)`; we do not change the fitted action.

Primary source: *SDMC Manuscript B, Structural Radius/Jacobian Investigation Appending Continuation*, Sections 26–36 (local mapper, isotropic trace potential, cubic source, nonlinear spherical collapse and JWST likelihood-equivalent transfer); and *Manuscript B Part I*, Section 56 (UV spatial k⁴ dispersion). Distinguish the **dimensionless local-map coefficient** `g3` from the **Horndeski cubic kinetic-braiding operator** `G3`; identical notation does **not** imply identical coupling.

### 2. The exact microscopic limitation: g3 is not determined

The adopted map potential on a local isotropic relative-metric trace is

```
V_eff/mu_V^4 = (1/2)u^2 + g3 W_m(p) h_loc kappa^3 u .
```

At its quasistatic minimum: `u = -g3 W_m h_loc kappa^3`. In perturbation theory `u` and `kappa` are first order; thus the source `u*kappa^3` begins **at fourth order in fields**. At `u=kappa=0`, the quadratic Hessian satisfies

```
V_uu = mu_V^4;
V_u,kappa = 3 mu_V^4 g3 W_m h_loc kappa^2 = 0;
V_kappa,kappa (from this term) = 6 mu_V^4 g3 W_m h_loc u*kappa = 0.
```

Accordingly, the homogeneous background and all linear observables have **zero local sensitivity** to `g3` under this operator. The positive quadratic potential stabilizes the local trace without choosing `g3`. `g3 = A_F` remains an explicit **hypothesis**, not a derived identity. The on-shell numerical fit cannot resolve this higher-order coefficient. There is also no demonstrated smooth covariant completion of the piecewise `h_loc > 0` switch or the phenomenological `W_m(p)` window; these need to arise from a local invariant action before a complete relativistic solver can be constructed.

The UV dispersion is distinct from Late300's Bellini–Sawicki kinetic alpha_K:

```
omega^2 = m_ell^2 + c_ell^2 k^2/a^2 + alpha_K^UV k^4/a^4.
```

`alpha_K^UV > 0` is a stability sign, **not** a coefficient determination. The microscopic curvature-selector coefficient `eta_K`, the trace-mode kinetic normalization, and the source normalization must be independently specified to compute it. The finite Late300 scalar kinetic function `D = alpha_K^EFT + (3/2)alpha_B^2` cannot be inserted as `alpha_K^UV` without deriving the relation between different action sectors.

This is an identifiability boundary, not a solver failure: a free coefficient multiplying a term invisible to lower-order observables cannot be inferred from them. A true prediction needs independent UV/symmetry input or nonlinear data likelihood with physical nuisance and selection priors.

### 3. New three-dimensional Jacobian calculation

For Zel'dovich/Lagrangian particle displacement `x=q-D grad(phi)` define symmetric tidal tensor `T_ij=partial_i partial_j phi`, eigenvalues `lambda_i`, and

```
J = det(I-D T) = product_i(1-D lambda_i)
  = 1-D I1 +D^2 I2 -D^3 I3.
```

Prior to shell crossing, `1+delta=1/J`, with no particle generation. At fixed trace `lambda_1+lambda_2+lambda_3`, the arithmetic–geometric mean inequality gives

```
J <= (1-D (lambda_1+lambda_2+lambda_3)/3)^3,
```

provided all three principal factors are positive. **Anisotropy can therefore yield stronger local density compression than the isotropic deformation with the same linear trace**. That is a geometrical statement, *not* proof of earlier virialized-halo or black-hole formation.

At `D=0.5` with the same trace `I1=1.5`:

| Tidal eigenvalues | J | (1+delta) | First-axis caustic D |
|---|---:|---:|---:|
| (0.5, 0.5, 0.5) | 0.421875 | 2.3704 | 2.0 |
| (1.2, 0.3, 0.0) | 0.340000 | 2.9412 | 0.8333 |

**Independent 40³ numerical diagnostic**, seeded random smooth Lagrangian Gaussian density field with no Late300 transfer function (N=40, seed=20261012):

- `64,000` Lagrangian cells; maximum principal deformation `D lambda_max=0.85`, so no crossed axes.
- Determinant/invariant maximum residual `1.33e-15`.
- The cellwise `J*(1/J)=1` check has max error `1.11e-16` (**algebraic**, not a simulation of dynamical continuity).
- Mean J = `1.0002377`, within 0.024% of periodic volume conservation; finite-grid aliasing remains.
- Mean `J_anisotropic/J_equal-trace = 0.94515`, with **no AM–GM violation**.
- 8.18% of cells have `delta > 1`; these percentages are **seed and synthetic-amplitude dependent**, not a prediction of halo statistics.
- Imposing the previous *spherical* curvature-source equation as an **anisotropic diagnostic only** produces median `u=-0.114` among active cells and 95th percentile `|u|=0.416`. This relatively large response flags **backreaction risk**, not a validated physical enhancement: the anisotropic Hamiltonian constraint has shear terms and the metric/scalar constraints were not solved.

Six grid/seed checks (24³, 32³, 40³, 48³, and alternative random seeds) retain `J` identity residuals below `1.4e-15`, the anisotropic inequality, and similar synthetic mean Jacobian ratios `0.934–0.945`. Files `results/robustness.json` and `test_sdmc_early_halo_closure.py` reproduce these checks.

**What remains to be solved before calling this fully coupled:** 3-D gravitational constraints, the spatial-map trace/traceless field equations, matter trajectories, baryonic gas physics, and consistent nonlinear propagators through shell crossing. An algebraic 3-D Jacobian—even one evaluated on a random field—is not a self-consistent 3-D N-body simulation.

### 4. New high-k screening and halo-host consistency test

The earlier report adopted an **illustrative comoving** filter:

```
T(k) = [1 + (k/k_c)^4]^-1,  k_c=3 Mpc^-1.
```

Map halo mass to a Lagrangian comoving radius using the Late300 present matter density:

```
r_L(M)=[3M/(4 pi rho_m0)]^(1/3), k_proxy=1/r_L,
rho_m0 = 2.77536627e11 * (omega_b+omega_cdm)
        = ~4.03e10 Msun/Mpc^3.
```

| Halo mass (M_sun) | k≈1/r_L (Mpc^-1) | T(k) fixed comoving k_c=3 | T(k,z=20) if cutoff instead fixed physically |
|---:|---:|---:|---:|
| 1e8 | 11.903 | 0.0040 | 0.00004 |
| 1e9 | 5.525 | 0.0800 | 0.00080 |
| 1e10 | 2.565 | 0.6519 | 0.0169 |
| 1e11 | 1.190 | 0.9758 | 0.2702 |
| 1e12 | 0.553 | 0.9989 | 0.8886 |
| 3e12 | 0.383 | 0.9997 | 0.9719 |

For a simple requirement `T(1e11 Msun)>=0.9`, `T(1e9 Msun)<=0.1`, the comoving cutoff must satisfy

```
2.06173 Mpc^-1 <= k_c <= 3.18990 Mpc^-1.
```

This provides a **bounded screening window**, conditional on those acceptance thresholds and the filter exponent, *not* a derivation of the microscopic `alpha_K^UV`.

To compare with the microscopic physical dispersion, ignore the k² term and define the physical crossover by `alpha_K^UV k_phys,c^4 = m_ell^2`. Then

```
k_c,comoving(a) = a (m_ell^2/alpha_K^UV)^(1/4).
```

If physical `m_ell²/alpha_K^UV` were time independent and `k_c=3` at `z=5.5`, then `k_c(z=20)=0.92857 Mpc^-1`: the `1e11 Msun` response falls to 0.270, **violating the 90% massive-halo criterion at z=20**. Keeping >=90% at z=20 instead requires the physical `m_ell²/alpha_K^UV` ratio to be >=24.3 times its z=5.5 normalization, under these exact threshold assumptions; holding *k_c=3 comoving at both epochs* would require a still stronger ~109-fold ratio evolution. Such time evolution has **not** been derived.

Moreover, the choice which suppresses dwarf halos also suppresses `1e8–1e9 Msun` candidate SMBH-seed hosts. The previous unscreened `z=20 -> z=23.26` collapse advance therefore **does not automatically apply to small seed-host masses**. Under the purely illustrative `A_M(k)=1+0.2123 T(k)` plus matter-era mapping `(1+z_new)=A_M(1+20)`:

| Mass | Advance in cosmic time from z=20, fixed comoving cutoff | Same with fixed physical cutoff |
|---:|---:|---:|
| 1e8 M_sun | 0.226 Myr | 0.002 Myr |
| 1e9 M_sun | 4.407 Myr | 0.045 Myr |
| 1e10 M_sun | 31.22 Myr | 0.95 Myr |
| 1e11 M_sun | 43.46 Myr | 14.19 Myr |
| 1e12 M_sun | 44.26 Myr | 40.34 Myr |

**Qualification:** These times are kinematic redshift shifts from a filtered algebraic `A_M`, not recalculated fully coupled collapse thresholds. The historical spherical result has a *different* response calibration (`A_M≈1.1555` for the old 0.01574 source), while `0.2123` represents the later Late300 likelihood-transfer normalization. They must not be silently equated. The diagnostic nevertheless exposes the main issue: **the galaxy-scale preferred filter and early low-mass SMBH seed-host boost do not automatically coexist**.

For a given epoch `z=5.5`, the exact normalization corresponding to the toy cutoff is

```
(alpha_K^UV/m_ell²) = a^4/k_c^4 = 6.9161e-6 Mpc^4,
```

but **neither** `m_ell` nor this ratio is predicted by the source action; the relation merely recasts the chosen cutoff.

See CI artifacts `results/screening_mass_redshift.png` for the screening figure.

### 5. JWST/FRESCO: what is and is not reproducible

The actual Late300 GitHub file `sdmc/run_g022a034_jwst_likelihood_equivalent.py` begins with an explicit warning that it is **not an exact catalogue-level JWST MCMC**. It samples the known published baseline efficiency-to-halo-response curve and adjusts it for the new Late300 linear spectrum.

The canonical JSON reports:

| Historical trial | g3 | A_eff | 95%-equivalent epsilon |
|---|---:|---:|---:|
| g3=A_F (hypothesis) | 0.020481465 | 1.2074346 | 0.3190523 |
| 20% efficiency **95% threshold target** | 0.0353433 | 1.36 target after transfer | ~0.20 target |

An older manuscript segment rounds the central likelihood-equivalent threshold to ~0.315 and reports `g3≈0.03482` for a previous exact-spectrum normalization. The **canonical candidate JSON** above is the correct version-specific input for the present table; mixing the two without their linear-spectrum correction introduces an avoidable discrepancy.

Independent public evidence: Comini, Vagnozzi & Loeb (2026), arXiv:2604.13866 (full Bayesian cosmological/stellar-efficiency analysis); Xiao et al., Nature 635, 311–315 (2024), DOI 10.1038/s41586-024-08094-5, 36 spectroscopically confirmed massive dusty galaxies at 5<z<9, including three extreme stellar masses. The original FRESCO paper states the advanced analysed data are available from the corresponding author on request; the raw survey images live at MAST JWST program 1895. The accepted Late300 repository does not contain a selection-complete, posterior-ready catalog covariance or a physical nonlinear SDMC halo mass function.

Thus **no new exact JWST posterior, p-value or Bayes factor is claimed here**. Reusing published epsilon thresholds to estimate a new epsilon threshold is a likelihood *transfer*, not a substitute for a joint catalog selection likelihood.

A useful independent mass-budget calculation uses `f_b = omega_b/(omega_b+omega_cdm) = 0.15221557` and `M_star = epsilon f_b M_halo`. For the three previously used FRESCO stellar-mass thresholds (1.096, 1.514, 2.344)e11 M_sun, the required host masses are:

| Stellar mass | Minimum host mass at epsilon=0.315 | Minimum host mass at epsilon=0.20 |
|---:|---:|---:|
| 1.096e11 | 2.286e12 | 3.600e12 |
| 1.514e11 | 3.158e12 | 4.973e12 |
| 2.344e11 | 4.889e12 | 7.700e12 |

These are baryon-budget identities **conditional on the efficiency**, not observed halo-mass measurements or catalog likelihood fits. Rare 1e12–1e13 M_sun halos at z~5.5 still require a genuine nonlinear abundance prediction.

### 6. SMBH mass-growth calculation with Late300 cosmic ages

Using the accepted fitted background integrated from the radiation era:

| Epoch | Age (Gyr) |
|---|---:|
| z=0 | 13.596894 (canonical result 13.596889) |
| z=6 | 0.921493 |
| z=7 | 0.753245 |
| z=10 | 0.465885 |
| z=20 | 0.175668 |
| z=23.26 | 0.141284 |

The historical **unscreened spherical** advance z=20 to z=23.26 provides an additional **34.384 Myr** of ordinary cosmic accretion time.

With classical Eddington-limited exponential growth,

```
M_BH(t) = M_seed exp[(t-t_seed)/t_ef]
 t_ef = t_E (epsilon_rad/(1-epsilon_rad))/(lambda_Edd * duty),
 t_E = c sigma_Thomson/(4 pi G m_proton) = 0.45049 Gyr.
```

At epsilon_rad=0.1, lambda_Edd=1, duty=1, t_ef=50.0545 Myr. Earlier seeding provides a maximum **1.99-fold mass advantage** at fixed later epoch if the two seeds have identical accretion histories. This is helpful but is not a replacement for a viable initial seed/gas supply.

For a `1e9 M_sun` target at observation redshift z=7:

| Accretion assumption | Minimum seed if formed z=20 | Minimum seed if formed z=23.26 |
|---|---:|---:|
| epsilon=0.10, lambda=1, duty=1 | 9.74e3 M_sun | 4.90e3 M_sun |
| epsilon=0.10, lambda=1, duty=0.5 | 3.12e6 M_sun | 2.21e6 M_sun |
| epsilon=0.15, lambda=1, duty=1 | 6.99e5 M_sun | 4.54e5 M_sun |

For a 1e9 M_sun target already at z=10, with epsilon=0.1 and continuous Eddington accretion, the minimum seed would be `3.03e6 M_sun` at z=20 or `1.53e6 M_sun` at z=23.26. Thus, **earlier halo collapse makes seed growth easier by about a factor of two under the unscreened reference but does not itself specify the physical seed mass or accretion history**.

The filtering test in §4 sharply qualifies this: at 1e8–1e9 M_sun the **screened** early-collapse advance may become negligible. The 34 Myr advantage should only be applied to an actual halo mass if the full nonlinear solver shows that halo receives it. No new matter production is required or assumed.

See CI artifacts `results/smbh_seed_constraints.png` for SMBH figure.

### 7. Required physical completion and falsification sequence

What would actually close the four requested tasks:

1. **Microscopic coefficient:** independently specify a smooth local covariant spatial-map/trace action including the source coefficient `g3` and positive spatial-curvature selector `eta_K`, together with the trace kinetic normalization. A parameter may be empirically constrained, but it cannot be advertised as derived just because the successful value `g3=A_F` was inserted.
2. **True 3-D dynamics:** evolve the **same** selected action through matter+metric+scalar+relative-spatial field constraints, include anisotropic shear and feedback, maintain the Hamiltonian constraints, ensure mass conservation until shell crossing, and proceed with resolved nonlinear gravity and baryonic subgrid physics. Test convergence at multiple resolutions and volume sizes. The present Zel'dovich verification is a geometry unit test, not this solution.
3. **Joint halo/galaxy observation:** derive the mass- and redshift-dependent collapse response and `dn/dM(z)` from that forward solution; incorporate stellar-mass scatter, selection completeness, volume, stellar-population systematics and correlated data to fit FRESCO+CEERS; require explicit likelihood inputs or access to the requisite advanced catalog. Do not convert the published epsilon transfer into an exact posterior by relabeling it.
4. **SMBH:** for the *screened* host mass distribution, calculate seed production per halo, fragmentation/cooling, infall, metal enrichment, angular-momentum transport, BH accretion rate and duty cycles, mass density and quasar luminosity functions. The current analytical Eddington bounds only constrain necessary seed/growth combinations.
5. **Global constraints:** any solution must preserve Late300's linear CMB/DESI/SN histories, BBN, weak lensing, halo counts at lower redshift, local tests and the galaxy/SMBH abundance simultaneously. A strong response in dwarf or linear sectors is a possible falsification, not a success.

**Verdict:** The geometric calculations, screening conditional window, and SMBH growth bounds are reproducible and completed here. A unique action-derived `g3`, a physical high-k cutoff, a fully coupled 3-D nonlinear cosmological run, and an exact JWST catalog posterior are **not closed with the present SDMC theory and accessible datasets**. The most important new physical result is the **mass/redshift screening compatibility test**: a single fixed cutoff tuned to enhance massive galaxies can negate the same collapse enhancement in small early SMBH seed-host halos. The next genuine first-principles model must account for this instead of assuming the unscreened z=20→23.26 advance applies at every mass.

### 8. Reproducibility

Files:

- `run_sdmc_early_halo_closure.py` — 3-D field diagnostic, screening, exact Late300 age quadrature and BH growth.
- `test_sdmc_early_halo_closure.py` — five unit tests; run `python -m unittest -v test_sdmc_early_halo_closure`.
- `plot_sdmc_early_halo_closure.py` — graphs of screening mass/redshift dependence and BH seed thresholds.
- `results/results.json` — input provenance and derived diagnostics.
- `results/robustness.json` — independent grid and random-seed checks.
- `results/screening_mass_table.csv`, `results/smbh_growth_table.csv`, `results/stellar_host_mass_table.csv` — editable outputs.
- `results/screening_mass_redshift.png`, `results/smbh_seed_constraints.png` — generated figures.

To reproduce: `python run_sdmc_early_halo_closure.py --output results --grid 40 && python -m unittest -v test_sdmc_early_halo_closure && python plot_sdmc_early_halo_closure.py` in a local Python environment with NumPy, SciPy and Matplotlib.

### Public source anchors

- SDMC Late300 candidate: https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/late300_candidate.json
- Historical conditional JWST script: https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/run_g022a034_jwst_likelihood_equivalent.py
- SDMC activation/covariant audit: https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/reports/activation_handoff_field_equations_2026-10-10.md
- Comini, Vagnozzi & Loeb (2026), arXiv: https://arxiv.org/abs/2604.13866
- Xiao et al., *Nature* 635, 311–315 (2024): https://doi.org/10.1038/s41586-024-08094-5
- FRESCO archive DOI: https://doi.org/10.17909/gdyc-7g80