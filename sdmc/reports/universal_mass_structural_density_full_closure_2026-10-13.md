# SDMC universal \(M_G/m_P=\kappa S\) hypothesis: full era, conservation, Jacobian, gravity and acoustic closure

**Investigation date:** 13 October 2026.  
**Run:** [Verified full scientific regression test #38327737402](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38327737402) (**SUCCESS**).  
**Status:** Necessary-condition mathematical and numerical closure; historical SDMC chronological \(S/S_0=t/t_0\) reference; **NOT** an independent microscopic \(S(t)\) derivation, a new primordial BBN integration, a new Planck likelihood or a new Hubble parameter determination.  
**Preserve accepted Late300 and ΛCDM local021. Investigation G remains on hold.**

## 1. Assumption, normalization and the exact inverse-square result

The *observed-present-epoch* Part II gravitational-source diagnostic is
\[
M_{G,0}=\frac{4\pi}{3}(\ell_P S_0)^3\rho_{m,0}\simeq9.74\times10^{53}\ {\rm kg}.
\]
With \(S_0=2.71678534\times10^{61}\) and \(m_P=2.176434\times10^{-8}\ {\rm kg}\),
\[
\kappa_0\equiv\frac{M_{G,0}}{m_PS_0}=1.647244742218.
\]
A normalization evaluated instead from precise accepted Late300 ordinary physical matter \(\omega_b+\omega_c=0.14507858642\) gives \(\kappa_{m0}=1.634411081486\). The difference is about 0.785% in **density**, reflecting historical rounding/normalization. Neither coefficient is independently predicted at all epochs.

**Global hypothesized continuation** \(M_G(S)=\kappa_0 m_PS\), together with exact \(R_{\rm struct}=\ell_PS\), gives
\[
\rho_G(S)\equiv\frac{M_G}{(4\pi/3)R_{\rm struct}^3}
=\frac{3\kappa_0}{4\pi}\frac{\rho_P}{S^2},
\quad \rho_P=\frac{m_P}{\ell_P^3},
\]
\[
\boxed{\rho_G/\rho_P=0.393234\ S^{-2}.}
\]
If one uses a spherical Planck reference \(\rho_{P,\rm sph}=3m_P/(4\pi\ell_P^3)\), the equivalent coefficient is \(\rho_G/\rho_{P,\rm sph}=\kappa_0/S^2\). The familiar rounded \(\kappa=1.64\) gives \((3\kappa)/(4\pi)=0.39152\); do not mix the *rounded* formula and *more precisely normalized* numerical tables.

This relationship is a consequence of three stipulated quantities: present *observed* matter mass, Planck units, and the structural radius. In particular
\[
\frac{M_{G,0}}{m_PS_0}
=\frac{4\pi}{3}\frac{\rho_{m,0}}{\rho_P/S_0^2}.
\]
Thus the apparent new constant **is algebraically fixed by the accepted present matter density**. Unless a new structural action predicts \(\kappa\) independently, this relation cannot by itself predict \(\omega_b,\omega_c\), \(R_b\) or \(H_0\).

At the Planck boundary \(S=1\), universal \(\kappa=1.64724\) would give \(M_G=1.64724m_P\), not the strict reference \(m_P\). A strict \(M_G(1)=m_P\) condition requires an epoch-dependent \(\kappa(S)\) or a revised initial-boundary definition.

## 2. Two distinct mass histories MUST NOT be conflated

Under the stipulated linear mass relation the structural-domain mass **changes**:
\[
\frac{M_G(S)}{M_{G0}}=\frac{S}{S_0}.
\]
For the historical chronological \(S/S_0=t/t_0\) anchored to fully integrated accepted Late300 \(H(z)\):

| Epoch | \(M_G(S)/M_{G0}\) | Implied effective structural mass |
|---|---:|---:|
| z=1090 | \(2.70451\times10^{-5}\) | \(2.63419\times10^{49}\ {\rm kg}\) |
| z=20 | \(0.0129197\) | \(1.25838\times10^{52}\ {\rm kg}\) |
| today | 1 | \(9.74\times10^{53}\ {\rm kg}\) |

By contrast, the conserved ordinary matter mass inside a **specified comoving patch** remains \(M_{\rm patch}(z)=M_{\rm patch,0}\) at all epochs. Therefore **one cannot simultaneously impose both \(M_G\propto S\) and fixed \(M_G=M_{G0}\) at all epochs**. The present Part II mass assignment does not settle which quantity describes all matter in the evolving structural region. If \(M_G\) is effective gravitational source mass, it need not equal conserved rest mass, but its stress-energy interpretation must be derived.

## 3. Mathematical equation-of-state closure and what 'era-dependent exponent' really means

Define \(p=d\ln S/d\ln a\). Then \(\rho_G\propto S^{-2}\) yields
\[
n_G(z)\equiv\frac{d\ln\rho_G}{d\ln(1+z)}=2p(z),\qquad
\rho_G(z)=\rho_{G0}\exp\left[2\int_0^z\frac{p(z')}{1+z'}dz'\right].
\]
The ideal match to an isolated conserved \(w\) fluid requires \(p=\frac32(1+w)\). This gives radiation (\(p=2,n=4\)), matter (\(p=3/2,n=3\)), and when \(p\approx1\), \(w_G\approx-1/3\), **not the \(w=-1\) vacuum equation of state**. This is a *conditional effective-fluid identification*, not SDMC's vacuum-pressure derivation.

For a conserved **matter-plus-radiation mixture**,
\[
n_{mr}=3+x_r,\quad x_r=\frac{\rho_r}{\rho_m+\rho_r},\quad
p_{mr,\mathrm{required}}=\frac{3+x_r}{2}.
\]
Matching the structural exponent to that function is an internally consistent **definition** of \(S_{\mathrm{required}}\), but not a new fluid physics derivation. For accepted chronological Late300 at z=1090, \(p_{\rm chrono}=1.681668\), whereas \(p_{mr,\mathrm{required}}=1.619580\). At z=3400 the corresponding numbers are 1.808754 and 1.747461. At z=20 they are 1.509140 and 1.503007. The mismatch is not zero through equality.

One can define a different \(S_{\rm required}(z)/S_0=\sqrt{\rho_{mr}(0)/\rho_{mr}(z)}\); at a given fixed \(\kappa\) this would force \(\rho_G\propto\rho_m+\rho_r\) by construction. To match absolute normalization and keep the present reference \(\kappa\), the required \(S\) differs from chronological \(S\) by roughly:
- z=20, \(\sqrt{\rho_{G,\rm chrono}/\rho_{mr}}\approx0.805\);
- z=1090, \(\approx0.899\);
- radiation era z=1e7, \(\approx1.070\).
It would require a genuinely modified structural lapse \(\mathcal N=t_P\dot S\) history. Its variation must follow from the action, not just be fitted to the matter and radiation densities.

## 4. Full numerical comparison of what the coefficient describes

**Full source and physical conservation** are distinct:
- \(\rho_m=\rho_{m0}(1+z)^3\) with \(\omega_b+\omega_c=0.14507858642\);
- \(\rho_r=\rho_{r0}(1+z)^4\) with \(\omega_r=4.17998772\times10^{-5}\);
- accepted restricted Horndeski \(\rho_X^J=F(z)\rho_{\rm crit,0} E(z)^2-\rho_m-\rho_r\);
- \(\rho_G\) is the **hypothetical** mass-derived structural density from \(\kappa\).

| z | \(\rho_G/\rho_m\) | \(\rho_G/(\rho_m+\rho_r)\) | \(\rho_G/(\rho_m+\rho_r+\rho_X^J)\) |
|---:|---:|---:|---:|
| 0 | 1.00785 | 1.00756 | 0.29480 |
| 20 | 0.65198 | 0.64806 | 0.64378 |
| 1090 | 1.06107 | 0.80731 | 0.79959 |
| 3400 | 1.85051 | 0.93465 | 0.92501 |
| 10000 | 4.04461 | 1.04203 | 1.03052 |
| \(10^7\) | 3298.20 | 1.14434 | 1.13083 |

**Decisive result:** With \(\kappa\) fixed and chronological \(S(a)\), \(\rho_G\) is neither the conserved matter density nor the total matter+radiation density at all epochs. Nor is it the entire Jordan-frame Friedmann source. Since at z=20 \(\rho_m>\rho_G\), a literal decomposition \(\rho_G=\rho_b+\rho_c+\rho_\gamma+\dots\) with nonnegative components fails. It may instead be an *effective structural gravitational source quantity*, not the actual sum of particles and photons.

A required time-dependent coefficient for exact \(\rho_G=\rho_m+\rho_r\) would be
\[
\kappa_{\rm required}(z)=\kappa_0\frac{\rho_m(z)+\rho_r(z)}{\rho_G(z)}.
\]
Values:
- z=0: **1.63488** (matter+tiny radiation only);
- z=20: **2.54181**;
- z=1090: **2.04042**;
- z=3400: **1.76242**;
- z=1e7: **1.43948**.

Hence **one universal constant coefficient does not equal the full physical matter+radiation inventory across all those epochs** under the historical chronology.

The alternative proposal \(\rho_G=\rho_m\) at every epoch predicts the independent Part I matter stretch
\[
S_m=(4\pi/3\kappa)^{1/3}S^{2/3},\quad
K_m=S/S_m=(3\kappa/4\pi)^{1/3}S^{1/3},
\]
whereas a constant matter/structural coordinate ratio would require \(\dot K_m=0\). Under actual separately conserved matter \(S_m\propto a\), so \(d\ln K_m/d\ln a=p-1\). These constraints coincide exactly only in specific eras or under a separately derived structural lapse evolution. Present \(K_m\sim2.20\times10^{20}\) matches Part I normalization, but at recombination it is only \(\sim6.48\times10^{18}\) for independent conserved matter, not constant.

## 5. A useful Planck-to-Friedmann identity

The new result can also be tested *without treating the structural density as baryons*. From \(\rho_P=m_P/\ell_P^3\simeq c^2/(G\ell_P^2)\), define \(\mathcal N=\dot R_{\rm struct}/c\) and \(p=\dot R/(HR)\):
\[
\boxed{\frac{\rho_G}{\rho_{\mathrm{crit,ref}}}
=\frac{2\kappa Gm_P}{c^2\ell_P}\left(\frac{c}{H R_{\rm struct}}\right)^2
=\frac{2\kappa Gm_P}{c^2\ell_P}\frac{p^2}{\mathcal N^2}.}
\]
Here \(\rho_{\mathrm{crit,ref}}=3H^2/(8\pi G)\); the Planck-unit factor \(Gm_P/(c^2\ell_P)\) equals unity up to rounding. For Jordan gravity with \(F\), divide by \(F\) in the effective Friedmann source ratio:
\[
\frac{\rho_G}{F\rho_{\mathrm{crit,ref}}}\simeq\frac{2\kappa p^2}{F\mathcal N^2}.
\]
**This identity was independently checked against numerical integration across all epochs; maximum residual \(6.66\times10^{-16}\)**.

At present, accepted \(F_0\simeq1.020525\) gives \(\rho_G/(\rho_m+\rho_r+\rho_X^J)=0.294799\). A hypothetical **structural-only** Friedmann closure with the same \(R_0,F_0\) would predict \(H_0\simeq37.85\ \mathrm{km/s/Mpc}\), NOT Late300's accepted \(69.7148\). Reproducing total present Jordan source with \(\kappa\) alone would require \(\kappa_{J,\mathrm{total},0}\simeq5.588\), not 1.64. This comparison simply identifies the missing late scalar/vacuum sector; it is not a new model or a fit.

At z≈1e7, the structural density is 1.144× ordinary radiation density. Two deliberately separate **invalid replacements** illustrate why the bookkeeping matters:
- If \(\rho_G\) REPLACES all matter/radiation/field sources, the background H value would differ by about **+6.3%**, even before recombination changes.
- If \(\rho_G\) REPLACES ordinary matter **but radiation is also counted**, the H value would increase by about **+46.0%**, effectively double counting a radiationlike structural term. This strongly violates previously verified BBN expansion history (accepted tracker +0.596%). **Neither is a valid autonomous model.**

Do not assign these H ratios to real SDMC without a new covariant stress-energy decomposition.

## 6. Baryon loading, sound horizon and H0: no independent new prediction

Even if \(\rho_G\) is used as a common reference density to express baryon and photon concentrations,
\[
b(z)=\rho_b(z)/\rho_G(z),\qquad g(z)=\rho_\gamma(z)/\rho_G(z),
\]
conservation gives
\[
\frac{d\ln b}{d\ln a}=2p-3,\quad
\frac{d\ln g}{d\ln a}=2p-4,\quad
\boxed{R_b=\frac34\frac{b}{g}
=\frac34\frac{\omega_b}{\omega_\gamma}\frac1{1+z}.}
\]
The structural factor \(S^{-2}\) cancels **exactly at all redshifts** in the ratio for the same physical densities. The accepted \(R_b(1090)=0.61486222\), sound horizon \(\sim143.25\) Mpc and fitted \(H_0=69.71482083\) km/s/Mpc remain the correct Late300 model outputs. There is **no additional term to apply** to the CMB sound horizon solely due to a re-expression of density via Planck stretch.

If one really forces physical \(\rho_b\propto S^{-2}\), the continuity equation needs \(Q_b/(H\rho_b)=3-2p\); if forcing radiation, \(Q_\gamma/(H\rho_\gamma)=4-2p\). These are genuinely new interactions, in conflict with the accepted minimally coupled Jordan-frame matter action and, generally, precision BBN. They were separately falsified as simple substitutions in [the preceding density/acoustic-only audit](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-structural-mass-baryon-loading-h0-test-20261013/sdmc/reports/structural_mass_baryon_loading_H0_2026_10_13.md). Their non-self-consistent acoustic inversion \(H_0\sim91.6\) km/s/Mpc is **not** a candidate to retain.

The successful frozen-covariant Late300 fit already includes changing radiation/matter/tracker/gravitational response in \(H(z),F(z),D(z)\). A genuine new CMB H0 would require a new, independently motivated action source or propagation mechanism, with tight-coupling dynamics, visibility, light-element production, full Planck/BAO/FS constraints and independent SN comparison consistently recomputed.

## 7. What is genuinely learned, and the next gate

**Accepted**: the mass-to-density algebra, the independent matter-stretch \(K_m\) relation, the exact Friedmann source identity, and continuous era exponent \(n_G=2p\).

**Conditional**: the *universal* epoch-independent mass coefficient \(\kappa\simeq1.64\), the historical constant-lapse chronology \(S\propto t\), the interpretation of \(\rho_G\) as any particular matter/radiation source. These must not be counted as confirmed microscopic predictions.

**Fails as an ordinary physical-density replacement in accepted Late300**: one constant \(\kappa\) together with chronological \(S(t)\) cannot represent conserved matter, radiation and the complete Jordan source across all epochs while preserving the existing BBN/CMB fits. Nor can \(M_G(S)\propto S\) be the same conserved comoving mass \(M_{G0}\) across all epochs.

**Still possible, but currently unverified**: \(\rho_G\) may be a distinct scale-dependent *effective gravitational source quantity* with a source projection \(K_m(\sigma,\mathcal N)\), distinct from the physical baryon/photon densities; only a fully specified Jordan/structural metric plus covariant action could determine this unambiguously. This overlaps the paused Investigation G but no work on G was resumed.

**Scientific verdict:** The calculation establishes a transparent, physically testable interpretation boundary but **does not yet provide a new \(R_b\), acoustic horizon or \(H_0\), and does not solve the Hubble tension**. The accepted benchmark is unchanged. Do not promote a new candidate from the structural coefficient alone.

## Reproducibility and evidence

- [Independent numerical run](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38327737402) (all source/dimensional/conservation checks successful).
- [Full Python source](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-mass-stretch-era-density-full-closure-20261013/sdmc/run_mass_stretch_era_density_full_closure.py).
- Workflow artifact \`sdmc-universal-mass-stretch-era-closure\` has full 14-epoch JSON and three scientific graphs: \`density_closure.png\`, \`epoch_structural_clock.png\`, \`mass_coefficient_closure.png\`.
- Manuscript I §10 independent cubic matter stretch and \(K_m\); Manuscript II §5.3 modern global-mass diagnostic; accepted Late300 parameters and the prior proven frozen-covariant model.
