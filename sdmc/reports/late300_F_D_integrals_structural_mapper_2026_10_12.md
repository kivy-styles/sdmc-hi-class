# Late300 gravitational and kinetic history integrals versus the SDMC structural clock and comoving mapper

**Date:** 12 October 2026  
**Scientific status:** Exact algebraic identities + independent numerical integration of the *accepted* Late300 target; separately labelled chronological structural diagnostic. No new physical parameter, automatic Jacobian force, new likelihood or independently solved UV gravitational action.

## 1. Three different time/coordinate variables must not be merged

Use \(n=\ln a\) for the conventional e-fold, \(\mathcal N=t_P\,dS/dt\) for the **structural lapse**, and \(S(t)\) for the structural stretch. The structural physical radius is \(R_{\rm struct}=\ell_PS\). Define \(\sigma=S/S_0\), \(p=d\ln S/dn\), and the homogeneous spatial comoving mapper \(\mathcal M=\sigma/a=R_{\rm struct}/(R_0a)\). Then

\[
\dot R_{\rm struct}=c\mathcal N,\quad
p={\mathcal N\over t_PSH},\quad
{d\ln\mathcal M\over dn}=p-1,\quad
S(t)=S(t_i)+{1\over t_P}\int_{t_i}^t\mathcal N(t')\,dt'.
\]

These relations are in Manuscript B and the structural-clock synthesis. Note \(\mathcal M\) is not the 3D Lagrangian mass-transport determinant, and \(\mathcal N\) is **not** \(n=\ln a\).

## 2. Exact accepted Late300 gravitational and kinetic target functions

Code \`sdmc/install_sdmc_v3_independent_kinetic.py\` explicitly sets

\[
U(n)={1\over2}[1+\tanh((n-n_c)/(2w_F))],\quad n_c=-\ln(1+z_c),
\]
\[
F(n)=\exp[A_FU(n)],\quad
\alpha_M={d\ln F\over dn}={A_F\over w_F}U(1-U),
\]
\[
D(n)=D_{\rm floor}+D_0U(n)^\nu,\quad
\alpha_B=-2\alpha_M,\quad D=\alpha_K+{3\over2}\alpha_B^2.
\]

For accepted \(\nu=1\), the exact accumulated-change integrals between arbitrary epochs are

\[
\ln{F(z)\over F(z_i)}
=\int_{\ln a_i}^{\ln a}\alpha_M(n')\,dn'
=\int_{z}^{z_i}{\alpha_M(z')\over1+z'}\,dz',
\]
\[
D(z)-D(z_i)={D_0\over A_F}\int_{z}^{z_i}
{\alpha_M(z')\over1+z'}\,dz'
={D_0\over A_F}\ln{F(z)\over F(z_i)}.
\]

Hence \(D=D_{\rm floor}+(D_0/A_F)\ln F\) with ratio \(D_0/A_F=16.713608919761782\). **This is a derived relation within the already fitted target; it does not independently predict that numeric ratio.** Nor is \(D\) an accumulated material energy: it is a dimensionless EFT kinetic combination.

The new numerical integrals from \(z_i=10^7\) to five selected redshifts agree with closed forms to \(\lesssim10^{-12}\). The accepted present values are \(F(0)=1.0205247753\), \(D(0)=0.3913854939\). Under the accepted No-Slip normalization, \(G_{\mathrm{eff}}(z)/G_{\mathrm{eff}}(0)=F(0)/F(z)\) in its stated linear response regime. This should not automatically be treated as a local laboratory measurement of Newton's constant.

## 3. Exact transformation into structural variables (not an independent gravitational law)

Since \(p=d\ln S/dn\) and \(d\ln\mathcal M/dn=p-1\),

\[
\alpha_M={d\ln F\over dn}
=p{d\ln F\over d\ln S}
=(p-1){d\ln F\over d\ln\mathcal M}
\]
along a well-defined monotonic trajectory (the last coefficient is undefined at \(p=1\) if \(\alpha_M\neq0\)).

The change integrals can therefore be rewritten
\[
\ln{F\over F_i}
=\int_{\ln S_i}^{\ln S}{\alpha_M\over p}\,d\ln S
=\int_{\ln\mathcal M_i}^{\ln\mathcal M}{\alpha_M\over p-1}\,d\ln\mathcal M.
\]
These are **coordinate changes of an already given gravitational response**, not a derivation from Jacobian kinematics. To *predict* \(F(S)\), an independent action must provide \(d\ln F/d\ln S\); to predict \(F(\mathcal M)\), it must provide \(d\ln F/d\ln\mathcal M\). The two rates need not be constant.

A structural field time integral exists independently:
\[
\sigma(t)=\sigma_i+(t_PS_0)^{-1}\int_{t_i}^t\mathcal N(t')\,dt'.
\]
Only if the action supplies \(G_4=M_{\rm Pl}^2F(\sigma)/2\) and \(\sigma\) is **physically** identified with \(S/S_0\) does its derivative give
\[
\alpha_M={\mathcal N\over Ht_PS_0}{F_\sigma\over F}
=p\sigma\,{F_\sigma\over F}.
\]
The accepted numerical linear-G3 reconstruction uses on-shell \(\phi=\ln a\); identifying that reconstruction coordinate with the underlying direct structural field requires justification.

## 4. Structural numerical diagnostic: chronological constant-lapse reference, NOT autonomous S dynamics

To explore the proposed coordinate connection, the new script uses the historical chronological approximation \(S/S_0=t/t_0\), hence a **constant \(\mathcal N\)** and \(p=1/(Ht)\), with \(t(a)=\int_0^a da/(aH_{\rm Late300})\) and \(\mathcal M=t/(t_0a)\). The full Late300 \(H(a)\) is retained across radiation, matter, scalar handoff and late expansion.

| z | F(z) | D(z) | alpha_M | approximate mapper M | p | beta_M=alpha_M/(p−1) |
|---:|---:|---:|---:|---:|---:|---:|
| 1090 | 1.00000000 | 0.0518155 | near 0 | 0.029506 | 1.68167 | near 0 |
| 30 | 1.00007902 | 0.0531361 | 0.00023770 | 0.222421 | 1.51159 | 0.00046 |
| 16.1932 | 1.00045997 | 0.0595015 | 0.00135755 | 0.300378 | 1.50852 | 0.00267 |
| 5 | 1.00731003 | 0.1735481 | 0.01417331 | 0.513197 | 1.50947 | 0.02782 |
| 3.92788 | 1.01029335 | 0.2229750 | 0.01546278 | 0.567153 | 1.50543 | 0.03059 |
| 2.2 | 1.01623873 | 0.3210432 | 0.01038656 | 0.701574 | 1.47334 | 0.02194 |
| 0 | 1.02052478 | 0.3913855 | 0.00049273 | 1 | 1.03153 | 0.01563 |

The mapper estimate near last scattering is close to, but *not identical* with, the previously reported Manuscript B \(\mathcal M_{1090}\approx0.0294675\), because the present calculation substitutes the Late300 background into the simplified chronology. Do not present it as the exact independently free-evolved structural field.

The required local response \(\beta_{\mathcal M}\) varies strongly from z=30 through z≈3.93 and today. Therefore a single constant \(F\propto\mathcal M^\beta\) is **not** implied by Late300. In this historical approximation \(\mathcal N\) is constant while \(F\) evolves by >2%, so varying structural lapse *alone* cannot be the necessary cause. The genuinely autonomous \(\mathcal N(t)\) history through all finite-era handoffs is still required for a universal clock-based force law.

## 5. The older N^{-3} gravitational-source projection is a distinct postulate

The *Part II Supplemental Closing Material* explicitly defines the **3D dynamical rate-space Jacobian** \(J^{(3)}_{\rm dyn}=\mathcal N^{-3}\) and a **local/global gravitational source projector** \(M_L=M_G K_m^{-1/2}\mathcal N^{-3}\).

That does **not** establish \(G_{\rm eff}\propto \mathcal N^{-3}\) or \(F\propto \mathcal N^3\). The source document calls the projection a structural postulate awaiting a covariant metric/action derivation. It is distinct from (i) the homogeneous comoving spatial mapper \(\mathcal M=\sigma/a\); (ii) the 3D local Lagrangian Jacobian; and (iii) the effective cosmological Planck mass \(F\).

As an explicit conditional falsification exercise only, identifying \(G_{\rm eff}\propto\mathcal N^{-3}\) with fixed \(K_m\) would require \(\mathcal N(0)/\mathcal N(z\gg z_c)=[F(0)/F(z\gg z_c)]^{1/3}\simeq 1.0067953\), a **0.68% rise** and \(d\ln\mathcal N/dn=\alpha_M/3\). This neither follows from the manuscript source projection nor agrees with its constant-lapse chronological submodel; do not elevate it to a new force law.

## 6. The genuine scalar-equation integration path

The accepted restricted Horndeski action is \(G_2=k_1(\sigma)X+k_2(\sigma)X^2-V(\sigma)\), \(G_3=g(\sigma)X\), \(G_4=M^2 F(\sigma)/2\), \(G_5=0\).

The already-derived conserved-current balance has the **genuine action-level** integral
\[
{d\over dt}(a^3\mathcal J_\sigma)=a^3P_\sigma
\quad\Rightarrow\quad
\mathcal J_\sigma(t)=a^{-3}(t)\Bigl[a_i^3\mathcal J_{\sigma,i}
+\int_{t_i}^t a^3(t')P_\sigma(t')\,dt'\Bigr].
\]
The current \(\mathcal J_\sigma=\dot\sigma[k_1+2k_2X-2g_\sigma X]+6HgX\). It is a **scalar current**, not a spatial Jacobian.

To determine a gravitational history rather than merely transform coordinates: independently supply the \(k_1,k_2,g,V,F(\sigma)\) functions and early structural \(\mathcal N\) boundary data, integrate this current equation together with the Friedmann/Raychaudhuri constraints and matter conservation, and compare resulting \(H,F,D,p,\mathcal M\) to accepted Late300. The on-shell reconstructed model was already shown to freely reproduce its selected target; the independent upstream selection of source functions remains the unresolved physics.

**Verdict:** The exact F/D redshift history integrals exist and pass independent checks. The structural clock and comoving mapper offer useful *indirect gravitational diagnostics*, with \(\alpha_M=p\,d\ln F/d\ln S=(p-1)\,d\ln F/d\ln\mathcal M\). Neither a varying lapse nor the spatial Jacobian **by itself** derives the cosmological gravitational running. The next decisive calculation is one physically independently specified \(F(\sigma)\)/kinetic/braiding action that predicts \(\alpha_M\), not one that reconstructs it from the already fitted \(F(z)\).
