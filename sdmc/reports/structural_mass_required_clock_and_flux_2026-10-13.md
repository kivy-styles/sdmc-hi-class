# SDMC continuation: derive the necessary structural clock from the Planck-to-global mass law

**Date:** 13 October 2026 · **Scientific classification:** independent mathematical/integration check of a proposed universal mass formula; accepted Late300 and local021 unchanged; **Investigation G still on hold** (this is a structural-mass/source consistency extension, not a new gravitational-field solution).

**Tested in GitHub Actions:** [Verified run #38349088042](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38349088042) (**SUCCESS**). Exact 16-epoch CSV-equivalent JSON, three high-resolution diagnostic graphs and source tests archived.

## 1. The missing physical interpretation: a moving structural sphere can enclose more conserved matter

Suppose the manuscript's *present-epoch* global diagnostic \(\kappa_0=M_{G,0}/(m_PS_0)\simeq1.6472447422\) is extended universally as a hypothesis \(M_G(S)=\kappa m_PS\), while \(R_{\rm struct}=\ell_PS\). Spherical volume implies \(\rho_G=\frac{3\kappa}{4\pi}\rho_P/S^2\).

For **separately conserved Jordan-frame dust**, \(\dot\rho_m=-3H\rho_m\), the ordinary *mass actually inside a time-dependent sphere of physical radius \(R(t)\)* is \(M_m=\frac{4\pi}{3}R^3\rho_m\). Its exact time derivative is
\[
\boxed{\dot M_m=4\pi R^2\rho_m(\dot R-HR).}
\]
This is the spherical-boundary flux of pre-existing matter relative to Hubble flow, **not matter creation**. The result is a bookkeeping identity; a structural "sphere" must actually be a physical domain with the Jordan-frame density and volume measure for this flux to be a physical mass flow.

Writing \(p=\dot R/(HR)=d\ln S/d\ln a\),
\[
\frac{\dot M_m}{HM_m}=3(p-1).
\]
If \(M_m\propto R\propto S\), it also requires \(\dot M_m/(HM_m)=p\). The two conditions coincide when **\(p=3/2\)**.

For mass-equivalent **radiation energy** in a sphere, \(\dot\rho_r=-4H\rho_r\). The energy changes by both boundary transport and photon redshift/pressure work:
\[
\boxed{\dot M_r=4\pi R^2\rho_r(\dot R-HR)-HM_r,}
\]
where \(M_r=\frac{4\pi}{3}R^3\rho_r\). Requiring \(M_r\propto S\) gives **\(p=2\)**. The extra \(-HM_r\) term is crucial; it would be wrong to equate photon energy with conserved baryon rest mass.

For the total ordinary matter+radiation energy,
\[
\frac{1}{H}\frac{\dot M_{mr}}{M_{mr}}=3p-(3+x_r),
\qquad
x_r=\frac{\rho_r}{\rho_m+\rho_r}.
\]
Therefore
\[
\boxed{M_{mr}\propto S\iff p_{mr}(z)=\tfrac12[3+x_r(z)]}.
\]
This is a **necessary and sufficient kinematic tracking condition** for a sphere enclosing homogeneous conserved matter+radiation to follow the proposed linear mass law, without creating matter. It is *not* an independent solution of the structural field equation.

## 2. Invert matter+radiation densities to reconstruct the *required* structural history

Accepted Late300 retains \(\rho_m=\rho_{m0}(1+z)^3\), \(\rho_r=\rho_{r0}(1+z)^4\). If its whole structural density is identified with **the matter+radiation inventory**—not the full scalar/vacuum Friedmann source—then, after matching the present normalization,
\[
\boxed{
\frac{S_{mr}(z)}{S_0}
=\sqrt{\frac{\rho_{m0}+\rho_{r0}}{\rho_m(z)+\rho_r(z)}}.
}
\]
This enforces **the same density slope** as matter+radiation at every era. Its actual structural lapse *required* by a given Late300 \(H(z)\) is
\[
\boxed{\mathcal N_{mr}=\frac{\dot R_{mr}}{c}
=\frac{p_{mr}H\ell_PS_{mr}}{c}.}
\]
The coefficient which makes the *absolute* present density equal to matter+radiation is \(\kappa_{mr}=1.6348819861\), versus the historically rounded Part II \(\kappa_0=1.6472447422\) (a 0.76% normalization difference). One must not claim the latter is *exactly* identical to the modern matter+radiation source.

| Redshift | Historical \(p\) | Required \(p_{mr}\) | Historical \(\mathcal N\) | Required \(\mathcal N_{mr}\) | \(R_{mr,\rm req}/R_{\rm historical}\) |
|---|---:|---:|---:|---:|---:|
| today | 1.031530 | **1.500144** | 3.413507 | **4.964232** | 1 |
| 20 | 1.509140 | **1.503007** | 3.413507 | **2.726490** | **0.801995** |
| 1090 | 1.681668 | **1.619580** | 3.413507 | **2.942705** | **0.895125** |
| 3400 | 1.808754 | **1.747461** | 3.413507 | **3.176270** | **0.963138** |
| 10000 | 1.909125 | **1.871183** | 3.413507 | **3.402409** | 1.016960 |
| \(10^7\) | 1.999885 | **1.999827** | 3.413507 | **3.637720** | 1.065715 |

**Key implication:** a dynamic lapse could in principle keep a constant structural-mass coefficient describing the moving matter+radiation content, but it must depart from the historical chronological constant-lapse approximation, especially today. It must be shown *independently* that the accepted field equations select this precise lapse.

At recombination the necessary \(R_{mr}\) is about **10.5% smaller** than the historical \(R_{\rm struct}\), and at z=20 approximately **19.8% smaller**. This is not an additional physical overdensity, because in this interpretation the domain mass changes with the moving boundary.

The **observable-universe comoving matter patch** does not obey \(M_{\mathrm{patch}}(z)\propto S\): for conserved dust its rest mass is constant, whereas \(M_{G,\rm structural}(z)/M_{G0}=S(z)/S_0\). Thus, the proposed evolving structural-sphere mass cannot simultaneously be the identical fixed inventory of the present observable progenitor patch at each redshift.

## 3. The full Jordan gravitational source demands a DIFFERENT history and normalization

The accepted covariant Late300 Friedmann constraint is
\[
\boxed{3M_{\rm ref}^2 F(z)H^2(z)=\rho_m+\rho_r+\rho_X^J\equiv\rho_J.}
\]
If the *entire* gravitational source \(\rho_J\) is instead equated to \(\rho_G=(3\kappa/4\pi)\rho_P/S^2\), the exact **required** structural evolution is
\[
\boxed{
\frac{S_J(z)}{S_0}=\frac{H_0}{H(z)}
\sqrt{\frac{F(0)}{F(z)}},
\qquad
p_J(z)=-\frac{d\ln H}{d\ln a}-\frac{\alpha_M}{2}.
}
\]
Use \(\mathcal N_J=p_J H\ell_PS_J/c\).

The present **absolute** source-normalization coefficient is **\(\kappa_J=5.5876836528\)**, over **3.4176 times** \(\kappa_{mr}=1.6348819861\). This makes the same fixed \(\kappa\simeq1.64\) impossible as the *entire* present Jordan source at the accepted radius/H0. No new fitted H0 follows.

| Epoch | Required \(p_J\) | Required \(\mathcal N_J\) | \(R_{J,\rm req}/R_{\rm historical}\) |
|---|---:|---:|---:|
| today | **0.472507** | **1.563606** | 1 |
| z=20 | 1.506172 | 5.034431 | 1.477762 |
| z=1090 | 1.619851 | 5.415095 | 1.646911 |
| z=3400 | 1.747834 | 5.842934 | 1.771371 |
| z=10⁷ | 1.999827 | 6.685337 | 1.958552 |

For a single structural field to be *both* \(\rho_m+\rho_r\) and the full Jordan source across epochs, their ratio would have to be constant after normalization:
\[
\frac{\rho_J(z)}{\rho_m(z)+\rho_r(z)}=\frac{\rho_J(0)}{\rho_m(0)+\rho_r(0)}
\quad \forall z.
\]
It demonstrably is not: at present this ratio is \(\kappa_J/\kappa_{mr}\approx3.4176\), whereas at recombination it is only \(\approx1.0097\). The two identifications cannot simultaneously be taken as universal physical equalities. In the accepted scalar-tensor action, the separate \(\rho_X^J\) sector is therefore essential.

## 4. A useful structural Hubble identity and why it does NOT independently solve H0

For any structural sphere with the linear mass-stretch law,
\[
\frac{8\pi G}{3}\rho_G=\frac{2\kappa Gm_P}{c^2\ell_P}\frac{c^2}{R_{\rm struct}^2}.
\]
If \(\rho_G\) is the **entire** Jordan source \(\rho_J\), then
\[
\boxed{
H^2=\frac{2\kappa Gm_P}{c^2\ell_P}\frac{c^2}{F R_{\rm struct}^2},
\quad
\frac{\mathcal N}{p}
=\sqrt{\frac{2\kappa Gm_P}{c^2\ell_P F}}.
}
\]
The numerical Planck factor \(Gm_P/(c^2\ell_P)\simeq1\). This is a **consistency identity**, not an independently obtained H(z): \(\kappa\) was inferred from the observed present density or chosen to fit the total Friedmann source.

At the accepted \(R_0\) and \(F_0\), using \(\kappa_0\simeq1.6472\) alone as the *complete* Jordan source would imply \(H_0\simeq37.85\,\mathrm{km\,s^{-1}\,Mpc^{-1}}\), not accepted **69.71482083**. Using \(\kappa_J=5.5877\) instead reconstructs the accepted H0 *by construction*, so it is not a new independent prediction or Hubble-tension solution.

Conversely, if \(\rho_G\) is **matter+radiation only**, the correct accepted equation is
\[
3M_{\rm ref}^2 F H^2=\rho_G+\rho_X^J
\]
with the genuine independent scalar/vacuum source retained. This is mathematically viable as a relabeling if \(\rho_G=\rho_m+\rho_r\) is enforced, but no new observable follows from the relabeling alone.

## 5. Baryon loading, CMB, BBN and SMBH

**Baryon loading:** in the matter+radiation identification, \(\rho_b,\rho_\gamma\) retain the independently conserved Jordan equations. Their common structural-reference normalization cancels exactly in \(R_b=3\rho_b/(4\rho_\gamma)\). Hence accepted \(R_b(1090)=0.6148622\), \(r_s\approx143.25\ \mathrm{Mpc}\) and **Late300 fitted \(H_0=69.71482083\)** are unchanged.

**BBN:** a newly *additional* \(\rho_G\propto S^{-2}\) radiation-like component would double count the early radiation energy; the accepted tracer is a small extra fraction \(\sim1.18\%\) and H speed-up \(+0.596\%\), not another \(O(1)\) background source. The inverse map to \(\rho_m+\rho_r\) remains a change of representation, not a new independently sourced BBN fluid.

**SMBH:** because \(M_G\propto S\) now has a possible moving-domain boundary-flux interpretation, a smaller global structural sphere does **not** necessarily contain the *same* conserved observable-patch matter mass at early epochs. It is invalid to cube global structural-vs-CDM patch radius ratios into halo overdensity/collision-rate gains unless a genuine local physical collapse solution separately establishes an equal-mass, smaller-volume halo. The existing nonlinear **local** SDMC mapper \(u\) remains the relevant candidate for such a mechanism.

## 6. Verdict and what is really needed next

The fixed-coefficient mass hypothesis has two coherent interpretations, but they cannot be identified without further physics:

1. **Moving matter+radiation structural domain:** \(\kappa_{mr}\simeq1.6349\), \(S_{mr}\) explicitly reconstructed from conserved fluid densities, \(p=(3+x_r)/2\), with a strongly changing required lapse. Physically possible as moving-domain bookkeeping, **not** a derivation of new CMB parameters.
2. **Full Jordan gravitational source:** \(\kappa_J\simeq5.5877\), \(S_J\propto [F H^2]^{-1/2}\), \(p_J=-d\ln H/d\ln a-\alpha_M/2\), with another drastically different lapse. This reconstruction already presupposes the observed/fitted H/F histories.

**Fundamental remaining gate:** derive from SDMC's independently specified action and structural initial conditions which (if either) dynamical clock is selected, with its matter/structural metric and energy-flux dictionary. The accepted restricted Horndeski reconstruction already reproduces Late300; merely reparameterizing \(\phi=\ln a\) as either required \(S(a)\) is *not* a fundamental prediction. The earliest Planck boundary is beyond this numerical Late300 extrapolation; neither required lapse can automatically be extended to \(S=1\), \(\mathcal N=1\).

**No accepted change to Late300, local021 or Investigation G. No new \(R_b\) or \(H_0\) is justified by these source inversions.**

### Reproducibility

- [Successful full GitHub numerical run](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38349088042)
- [Exact reproducible source](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-structural-mass-required-clock-20261013/sdmc/run_structural_mass_required_clock.py)
- [Prior mass-to-density audit](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-mass-stretch-era-density-full-closure-20261013/sdmc/reports/universal_mass_structural_density_full_closure_2026-10-13.md)
- Artifact contains full JSON, three figures and this report; use current run's artifact ID to download.
