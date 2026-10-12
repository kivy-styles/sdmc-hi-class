# SDMC Investigation L — Primordial Lithium, \(^{7}\mathrm{Be}\) and the early structural tracker

**12 October 2026** · **Research status:** numerical audit COMPLETE, with a successfully executed historical nuclear-network cross-check; modern PRIMAT microscopic action calculation still outstanding. **Investigation G is ON HOLD**; the gravitational-history/structural-clock/Jacobian research and accepted Late300 action are unchanged.

## 1. Existing manuscript and observational basis

- SDMC **Part IV**, **Part VIII §5**, **Part IX §125** already identified the standard lithium excess and concluded the preferred small tracker cannot fix it through a uniform BBN-clock change without endangering D/H and \(^{4}\mathrm{He}\).
- **Part VIII's PRIMAT-calibrated *local sensitivity screen*** at \(\lambda_e=20\), \(\omega_b=0.02239952\), \(f_r=.01\) produced baseline Yp=.246990, D/H=2.43955e-5, He3/H=1.04049e-5, Li7/H=5.54707e-10; tracker Yp=.247820, D/H=2.45993e-5, He3/H=1.04344e-5, Li7/H=5.51555e-10. Percentage changes: +0.336%, +0.835%, +0.283%, **−0.568%**. These are NOT a direct full modified-background PRIMAT run, and \(\lambda=20\) is NOT accepted Late300's \(\lambda_e\).
- Accepted [Late300 candidate data](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/late300_candidate.json) fixes \(\lambda_e=18.40625\), \(\omega_b=0.0220832191946\), \(H_0=69.71482\), early-tracker midpoint \(z_t=16.19\). Since nuclear BBN occurs at \(z\sim10^8-10^{10}\), neither the *late* tracker shutoff near \(z=16\) nor Planck-mass midpoint \(z_c=3.93\) is occurring during nucleosynthesis.
- PDG 2025 Big Bang Nucleosynthesis review, §24.3/24.5: observed halo star Li/H \((1.45\pm0.25)\times10^{-10}\) should be viewed as a possible *lower bound* due to depletion; a representative CMB-baryon predicted Li/H \((4.72\pm0.70)\times10^{-10}\) exceeds it by a factor 3.3. PDG mean D/H \((2.508\pm0.029)\times10^{-5}\) and representative helium \(Y_p\sim0.2453\pm0.0034\). This is not proof the *primordial* Li abundance is actually lower; stellar depletion is a viable primary candidate.
  - https://pdg.lbl.gov/2025/reviews/rpp2025-rev-bbang-nucleosynthesis.pdf (published 1 December 2025).
- The dominant lithium-7 source is \({}^{3}{\rm He}(\alpha,\gamma){}^{7}{\rm Be}\) followed by later electron capture \({}^{7}{\rm Be}+e^{-}\to{}^{7}{\rm Li}+\nu_e\), with more direct Li produced and some destroyed earlier. When interpreting a BBN nuclear-network output, **count surviving direct \({}^{7}{\rm Li}\) PLUS surviving \({}^{7}{\rm Be}\)** for eventual Li-7/H. Selective destruction would involve \({}^{7}{\rm Be}(n,p){}^{7}{\rm Li}(p,\alpha)2\alpha\) or other actual reactions; neutron injection generally also perturbs D/H and He, and is not in the accepted SDMC action.

## 2. Recompute the accepted early expansion shift (not a new gravitational-history model)

At radiation tracking \(\rho_{\rm struct}/\rho_{\rm total}=f_r=4/\lambda_e^2\), hence
\[
H_{\rm SDMC}/H_{\rm std}=1/\sqrt{1-f_r}.
\]
Accepted Late300 yields:
- \(f_r=0.01180672256796\), a 1.18067% fractional early tracker.
- \(H_{\rm SDMC}/H_{\rm std}=1.005956155492\), **0.5956% speed-up** at fixed standard radiation density.
- Characteristic radiation-era duration factor \(\sqrt{1-f_r}=0.9940791\), about **0.5921% less nuclear processing time**.
- Radiation-equivalent local proxy near \(\sim1\) MeV, \(\Delta N_{\rm eff}^{\rm bg}\approx(43/7)f_r/(1-f_r)=0.07339355\). This is a **background expansion equivalence**, not evidence of an actual extra neutrino species, and electron–positron annihilation makes the conversion epoch dependent.
- \(\omega_b=0.0220832191946\) corresponds to \(\eta_{10}\simeq 6.04742\) under the PDG baryon-density conversion \(\eta_{10}=6.115\,\omega_b/0.02233\).

**No structural/comoving Jacobian directly rescales nuclear densities, cross sections or temperature**: the universal Jordan matter frame conserves baryons \(\rho_b\propto a^{-3}\), and the Jacobian \(\mathcal M=(S/S_0)/a\) is a change of coordinate dictionary. A new physical interaction would have to be derived to alter reaction rates.

## 3. Updated local-sensitivity estimate with transparent assumptions

Start from the earlier PRIMAT-calibrated reference and rescale the **small** tracker shifts by \(f_r/.01=1.1806723\), correcting the changed baryon density using local standard-BBN sensitivity \({\rm D/H}\propto \omega_b^{-1.6}\), eventual \({}^{7}{\rm Li/H}\propto\omega_b^{2}\). Approximate \(Y_p\) baryon slope \(dY_p/d\eta_{10}\simeq.0016\). These are **approximate response laws**, not exact nuclear-network predictions.

| Response quantity | Preliminary accepted Late300 | Observational reference |
|---|---:|---:|
| \(Y_p\) | **0.247831** | representative 0.2453±0.0034 |
| D/H | **2.52030×10⁻⁵** | 2.508±0.029 ×10⁻⁵ |
| eventual Li7/H | **5.35536×10⁻¹⁰** | stellar 1.45±0.25 ×10⁻¹⁰ |
| Li excess relative to inferred stellar central value | **3.693×** | should allow depletion |
| tracker-only Li reduction at fixed baryon density | **~0.671%** | ~73% total reduction needed after new baryon rescaling |

Under a purely multiplicative phenomenological selective-7Be-depletion scenario, reducing 5.355e-10 to the 1.45e-10 reference would require removing **72.92%** of surviving total Li7+Be7 progenitor inventory. The corresponding *hypothetical homogeneous selective destruction optical depth* is \(-\ln(1-0.72924)=1.30653\). **This is NOT a calculated SDMC reaction rate or mechanism.**

A pure baryon-density adjustment with Li\(\propto\eta^2\) would require \(\eta\) reduced to about **0.513×**, which in the same approximate model would enhance D/H by **2.91×**, strongly contradicting its high-precision observations. Uniform early expansion and baryon changes are therefore not a plausible lithium-only cure.

## 4. Direct nuclear-network stage (independent historical code)

A new branch specifically for this investigation is created: \`sdmc-primordial-lithium-investigation-L-20261012\`. The workflow
\`.github/workflows/sdmc_primordial_lithium_network.yml\` pins \`espensem/AlterBBN_master@8219dd4b578244e36d56ebda31541705ae382093\` (a legacy 2016 *modified AlterBBN* implementation), builds its **actual C nuclear network** from source, and executes seven controlled experiments:
1. local021 ordinary \(\omega_b\) and zero tracker;
2. Late300 ordinary \(\omega_b\) with zero tracker (same-baryon control);
3. Late300 accepted radiationlike tracker \(f_r=0.0118067\);
4. early tracker \(f_r=.01\);
5. \(f_r=.0177778\);
6. \(f_r=.04\);
7. \(f_r=.0625\).

All cases use identical network/neutron lifetime \(879.4\) s/\(N_\nu=3.046\)/nuclear rates. The final validated solver implements a **fixed fraction of the actual instantaneous physical radiation/plasma energy**, \(H^2/H_{\rm std}^2=1/(1-f_r)\); zero entropy transfer and zero direct nuclear couplings. Across electron–positron annihilation this is only an *approximate* representation of a constant fractional SDMC tracker. The \(\nu\)-decoupling and neutrino transport of a modern PRIMAT run may differ.

**Do not substitute the proxy for a genuine derived SDMC \(H(T)\)**, and do not claim that a historical alternative network is as precise as modern PRIMAT or PArthENoPE. Report direct Li7, Be7, eventual total, D, He3 and Yp separately; observe their mutual tradeoffs. The workflow must actually succeed before its numerical outputs are claimed as executed.

## 5. Independent selectivity pathways and required tests

**(A) Stellar depletion (presently strongest non-SDMC explanation).** Test observational plateau slope/scatter, \(^{6}\)Li constraints, binary and diffusion signatures, young metal-poor stars and interstellar Li. This does not require modifying accepted SDMC. The 2025 PDG review notes increasingly strong depletion evidence, but astrophysical systematics remain.

**(B) Modified baryon/nuclear binding dynamics.** SDMC's \(F(\sigma)\) changes the effective gravitational response in the cosmology; it does **not automatically** vary quark masses, binding energies, \(n\leftrightarrow p\) weak rates, or electromagnetic nuclear cross sections. To test, derive explicit Jordan-frame couplings that satisfy laboratory constraints, then run changed nuclear Q-values and reaction rates. Must compare D/H, Yp, CMB/BBN \(\eta\), not just Li7. These couplings are **new hypotheses**, not part of accepted Late300.

**(C) Neutron-assisted \(^{7}\)Be destruction.** Quantify actual neutron source, energy deposition and nonthermal spectra; D/H rises in many models, and helium can be modified. Accepted SDMC explicitly conserves matter/radiation after thermalization: there is **no independent post-BBN neutron injection**. No new matter-creation hypothesis is retained. This route requires an explicit physical interaction and is currently unmotivated for accepted Late300.

**(D) Primordial chemical or inhomogeneous baryon processing.** A pure comoving mapper cannot raise *local physical* baryon density. A real finite-scale perturbation in physical baryon density, potential or reaction rate would require action-level perturbation equations, CMB isocurvature and D/H consistency, and must not double count the Jacobian. No such primordial microphysics is presently validated.

## 6. Preliminary verdict and gate to finish

The existing manuscript negative verdict for a homogeneous tracker survives direct scaling to current Late300 values: the uniform BBN-clock adjustment is O(0.6%) and the desired lithium-specific correction O(73%). The physical signal needed is over two orders of magnitude larger *in impact on Li* than the already calibrated small tracker effect, which cannot be achieved merely by increasing \(H\) without confronting D/H and Yp. The large-parameter extrapolation would be invalid.

The independent historical BBN network is designed to quantify this tradeoff and validate the *direction* of the effect. Full closure requires a modern nuclear-network run with a real dynamically derived \(H_{\rm SDMC}(T)\), neutrino transport, matched \(\omega_b\) and rate uncertainties. No result should be called a solution unless the lithium central tendency moves decisively while D/H and helium remain jointly compatible with independent observations.

**Current hypothesis preference:** preserve the accepted successful SDMC light-element background; investigate stellar depletion and precision BBN astrophysics first; do not alter late300 to force lithium concordance.

## Sources

- SDMC Part IV §8, Part VIII §5, Part IX §125, plus 10 Oct covariant activation/handoff audit.
- [2025 Particle Data Group, *Big Bang Nucleosynthesis*](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-bbang-nucleosynthesis.pdf), §§24.3–24.6.
- [AlterBBN v2 code paper](https://arxiv.org/abs/1806.11095) and [pinned historical AlterBBN fork](https://github.com/espensem/AlterBBN_master/tree/8219dd4b578244e36d56ebda31541705ae382093).
- [2025–26 \(^{7}\)Be nuclear-reaction measurements](https://www.sciencedirect.com/science/article/pii/S0375947425002350).


---

## 7. COMPLETED: independently integrated matched nuclear network, with the correct H² profile

**Final successful execution:** [GitHub run 38291398803](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38291398803), \`success\`. Seven central-rate nuclear-network integrations were completed and validated, including all six original radiation-era density/expansion calculations of the pinned legacy AlterBBN engine.

### Source-level technical corrections verified

The first serial compilation succeeded but revealed that returned \`ratioH[8]\` already includes the \({}^{7}\mathrm{Be}\to{}^{7}\mathrm{Li}\) post-BBN conversion (\`bbn.c\` lines 1196 and 1995). **The correctly reported eventual Li7/H is returned \`ratioH[8]\`, not \`ratioH[8]+ratioH[9]\`.** Direct Li7 can be reported as total minus separate surviving Be7.

A second correction was essential: the upstream \`dark_density\` routine and Hubble equation use inconsistent historical energy-density unit conventions, and a fixed \(T^4\) source does NOT guarantee constant fractional radiation energy through electron–positron annihilation. The **final** version patches **six H evaluations** only, replacing dark-density additions with an instantaneous scaling
\[
H_{\rm SDMC}^2(T)=H_{\rm standard}^2(T)\,[1+\kappa],
\quad \kappa=f_r/(1-f_r).
\]
The code retains the actual time-dependent photon, electron–positron, neutrino and baryon energy density at each step. Nuclear reaction rates and nuclear abundance ODEs are **not** changed; no lithium-specific sink, matter creation or unknown gravitational projection is introduced. This is a source-pinned constant-fraction BBN background approximation, not an autonomous SDMC field equation or modern PRIMAT.

**The earlier workflow 38290688144**, before isotope/source correction, **double-counted Be7 and underweighted expansion and must NOT be used**. **Workflow 38291079141**, after correcting source units but using a fixed \(T^4\) relative density, was a *different* model with an artificially growing relative dark fraction at low T; its Li7 −3.4% is not the preferred SDMC constant-tracker result. Only **38291398803** passes the final physical-background gate.

### 7.1 Direct BBN nuclear-network abundances

Same baryon-density control uses accepted \(\omega_b=0.022083219194622913\), \(\eta_{10}=6.04741985558\), \(N_\nu=3.046\), \(\tau_n=879.4\) seconds, the same historical nuclear rates and initial conditions.

| BBN abundance | Control at same \(\omega_b\), zero tracker | Accepted Late300 tracker | Relative change |
|---|---:|---:|---:|
| Helium-4 mass fraction \(Y_p\) | 0.246937143567 | **0.247984924126** | **+0.424311%** |
| Deuterium D/H | \(2.486886577338\times10^{-5}\) | **\(2.515294512460\times10^{-5}\)** | **+1.142309%** |
| Helium-3/H | \(1.038285125463\times10^{-5}\) | \(1.041693434059\times10^{-5}\) | +0.328263% |
| Direct Li-7/H at BBN freeze-out | \(2.906474460896\times10^{-11}\) | \(2.937299944152\times10^{-11}\) | +1.060580% |
| Surviving Be-7/H precursor | \(5.001246594067\times10^{-10}\) | \(4.960320399345\times10^{-10}\) | −0.818320% |
| **Eventual total Li-7/H** | \(5.291894040157\times10^{-10}\) | **\(5.254050393761\times10^{-10}\)** | **−0.715125%** |

**Almost 94.4% of eventual Li7 was synthesized as Be7**, with the remaining ~5.6% being surviving direct Li7. As predicted, a faster expansion mildly decreases surviving Be7 while increasing direct Li7, and the net Li7 change is far too small to solve the cosmological lithium excess.

The direct-network lithium estimate is \(\sim3.62\times\) the PDG 2025 metal-poor stellar central value \((1.45\pm.25)\times10^{-10}\). A further decrease of \(\sim72.40\%\) would be needed to meet that central value if no stellar depletion occurs. A merely illustrative selective irreversible destruction optical depth would need \(\tau_{\rm destroy}=\ln(5.25405/1.45)=1.28744\). Neither such a source nor such a rate follows from accepted Late300.

At the chosen PDG observed D/H \((2.508\pm.029)\times10^{-5}\), the Late300 background-proxy central result differs by 0.25 quoted observational standard deviations. The helium result differs by about 0.79 of the *representative* PDG helium observational uncertainty \(0.0034\). These are **not a complete likelihood**: they omit modern network rate covariance, helium systematic pooling and BBN/neutrino solver differences. The historical-engine local021 D/H differs by ~2.8 observation-only standard deviations, demonstrating why **this code cannot be used for high-precision cosmological model selection** against a full modern nuclear network.

### 7.2 Radiation-fraction stress test at Late300's accepted \(\omega_b\)

| \(f_r\) used as a constant-fraction background proxy | H speed-up | \(Y_p\) | D/H (\(10^{-5}\)) | Eventual Li7/H (\(10^{-10}\)) |
|---|---:|---:|---:|---:|
| 0 | 0% | 0.246937 | 2.486887 | 5.291894 |
| 0.010000 (\(\lambda=20\)) | 0.504% | 0.247824 | 2.511416 | 5.260036 |
| **0.011807 (accepted \(\lambda=18.40625\))** | **0.596%** | **0.247985** | **2.515295** | **5.254050** |
| 0.017778 (\(\lambda=15\)) | 0.901% | 0.248520 | 2.527829 | 5.234477 |
| 0.040000 (\(\lambda=10\)) | 2.062% | 0.250541 | 2.599495 | 5.117408 |
| 0.062500 (\(\lambda=8\)) | 3.280% | 0.252638 | 2.647469 | 5.045230 |

At high tracker fractions, deuterium and helium are pushed progressively upward while lithium remains nowhere near the stellar central value. These are **stress-test parameterizations**, not observationally accepted SDMC alternatives.

### 7.3 Why the two independent numerical estimates differ

The original Part VIII PRIMAT-*calibrated local sensitivity* at \(\lambda=20\) predicted a **−0.568%** Li response and +0.835% D response. The pinned legacy direct nuclear-network at the *current* matched Late300 baryon density and \(\lambda=20\) gives **−0.602%** Li and **+0.986%** D. At accepted \(\lambda=18.40625\), this network gives **−0.715%** Li and **+1.142%** D, compared with extrapolated PRIMAT-calibrated response around **−0.671%** Li and **+0.986%** D (different reference baryon density). 

Thus the direction and scale of independent calculations agree, while differences reflect source code vintage, radiation modeling and baryon-density normalization. The appropriate scientific conclusion is not that SDMC solved lithium but that its small uniform BBN tracker remains compatible with basic D/He while making only percent-level adjustments to Li.

### 7.4 Exact reproducibility and preserved status

- [Completed source-pinned independent BBN workflow](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38291398803), passing all internal nuclear, isotope-accounting and parameter-matching checks.
- [Source-level instantaneous-background patch](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-primordial-lithium-investigation-L-20261012/sdmc/lithium/patch_alterbbn_instantaneous_tracker_H.py).
- [Nuclear-wrapper with exact post-capture accounting](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-primordial-lithium-investigation-L-20261012/sdmc/lithium/run_legacy_alterbbn_tracker.c).
- [Archived complete seven-run test sheet](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38291398803) in \`investigation-L-Late300-Li7-Be7-deuterium-helium-network\`.

**Final current scientific verdict: SDMC Late300's existing accepted uniform early tracker does NOT solve primordial lithium.** The nuclear network verifies the limitation quantitatively. Part IV and Part IX's verdict is still valid. Investigation G remains on hold. A separate matter-creation theory or a forced post-BBN neutron source is not introduced.
