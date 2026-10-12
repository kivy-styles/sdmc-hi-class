# SDMC Investigation L — Primordial Lithium, \(^{7}\mathrm{Be}\) and the early structural tracker

**12 October 2026** · **Research status:** early analytical audit plus independent historical nuclear-network cross-check. **Investigation G is ON HOLD**; the gravitational-history/structural-clock/Jacobian research and accepted Late300 action are unchanged.

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

All cases use identical network/neutron lifetime \(879.4\) s/\(N_\nu=3.046\)/nuclear rates. The radiationlike scalar fraction is implemented as a **separately conserved \(T^4\) dark-density proxy**, with \(\kappa_\rho=f_r/(1-f_r)\) and \(\rho_D\propto T^4\); zero entropy and zero direct nuclear couplings. Across electron–positron annihilation this is only an *approximate* representation of a constant fractional SDMC tracker. The \(\nu\)-decoupling and neutrino transport of a modern PRIMAT run may differ.

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
