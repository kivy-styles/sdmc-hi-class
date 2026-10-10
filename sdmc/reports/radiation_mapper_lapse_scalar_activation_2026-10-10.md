# SDMC: Does the structural/coasting radius deficit determine radiation-era scalar activation or matter-seed transfer?

**10 October 2026 — Mathematical continuation of Manuscript B and the existing early-galaxy/SMBH audit.**  
**Scientific status:** the mapper establishes the structural exponent, the manuscript's density-route identity fixes a *conditional* scalar activation when a scalar fraction and lapse are supplied, but no independent radiation-to-CDM transfer operator or new candidate is derived. **Frozen late300 unchanged.**

Sources: Manuscript B Part I *Causal Genesis and Thermal Handoff*, §§2, 10, 13–16; Manuscript B spatial map §§50–58; structural radius/Jacobian study §§2–3 and 10–12; independent activation/handoff report `sdmc/reports/activation_handoff_field_equations_2026-10-10.md`; late300 JSON and `sdmc/audit_radiation_mapper_activation.py`.

## 1. Three distinct ratios: compare like with like

Define `x=ln a`, `sigma=S/S0`, `a_R=sigma`, normalized physical structural radius `R_str=R0*sigma`, and `M=sigma/a`. Then

`chi_R=R_str/sigma=R0` (structural comoving),  
`chi_J=R_str/a=M R0` (Jordan/Jacobian comoving),  
`R_str/chi_J=a` (physical-to-comoving ratio),  
`R_str/(a R0)=chi_J/R0=M` (the matched, same-epoch coasting-reference ratio).

The last identity *validly* measures the relative structural stretch compared with `sigma_coast=a` normalized at today's radius, but is **not** a force or scalar-perturbation suppression law. At recombination `M*=0.02946748423`, `chi_J*=~420 Mpc`, and `R_str*=~0.385 Mpc`. Comparing `0.385/420` measures the ordinary scale factor `a*=1/1091`, not `M`. A physical / comoving comparison cannot identify an extra gravitational compression.

The normalized ratio `M` is fixed by historical `sigma` and `a`; it should not be double-counted in a Jordan-frame density already transformed correctly.

## 2. Insert the mapper in the manuscript's scalar activation equation

The retained relation is `rho_X=chi_act rho_v`, `rho_v∝sigma^-2`, and `p_s=dlnS/dln a`. Therefore

`**p_s = 1 + dlnM/dx**` and `dlnrho_v/dx = -2p_s`.

The manuscript's restricted density-route lapse equation is

`N_s = p_s Xi_v sqrt(chi_act/Omega_X), Xi_v=sqrt(8pi/3)`,

with `N_s=dt_s/dt_c`, **not** `x=ln a`. Hence the algebraic inverse is

`**chi_act=Omega_X (N_s/(Xi_v [1+dlnM/dx]))^2.**`

This is an identity within this specified density-route normalization. It does **not** imply `chi_act=M`, `chi_act=M²`, or that `chi_act` is the physical gravitational coupling in a generic Jordan-Horndeski model.

On the classical radiation tracker, `p_s=2`, `N_s=3.41350846` and manuscript `Omega_X=0.01`. Then

`chi_act = 0.0034771495577581`.

By contrast, the normalized *recombination* geometric mapper is `M*=0.02946748423`; this belongs to a **later epoch** and a different definition. Even in an illustrative identical-epoch comparison the two values are not related as a universal suppression. Keeping `p_s,N_s` constant but algebraically changing `Omega_X=0.005, 0.01, 0.02` gives `chi_act=0.0017385748,0.0034771496,0.0069542991`, demonstrating nonidentifiability from the mapper alone; varying fractions is **not** an independently evolved alternate SDMC solution.

The alternative steep radiation tracker fraction `Omega_X≈4/lambda_e²` with `lambda_e=18.40625` gives `Omega_X≈0.0118067` and `chi_act≈0.00410537` under the *same assumed* N_s and p_s. This 0.41% value is another conditional model-stage estimate, not a new measured radiation fraction.

## 3. Does the mapper determine the radiation transfer J?

Write `Achi=dlnchi_act/dx`. Differentiating the inverted lapse identity gives the exact kinematic identity

`**Achi = dlnOmega_X/dx +2 dlnN_s/dx -2 dlnp_s/dx.**`

Manuscript B's continuity/thermal bridge with `q_J=J/(H rho_X)` is

`**2p_s=Achi+3(1+w_X)+q_J.**`

Eliminating Achi yields a **constraint**, not an independent transfer prediction:

`**q_J=2p_s-3(1+w_X)-dlnOmega_X/dx-2dlnN_s/dx+2dlnp_s/dx.**`

Why it does not close: independent action equations are still required for `Omega_X(x),N_s(x),p_s(x)` and the bath decay current. Reconstructing q_J from the motion, when that motion is already calibrated to a handoff, is not a microphysical derivation of J.

Special exact-radiation tracker with `p_s=2,N_s=constant,Omega_X=constant,w_X=1/3` has `Achi=0` and `q_J=0`. Bare `rho_v∝a^-4` and active `rho_X∝a^-4` track ordinary radiation. This does **not** preclude the earlier *thermal-production pulse* before the radiation tracker became established. Manuscript B already supplies a candidate mediator/bath interaction during that *earlier* stage; it does not determine a later CDM-perturbation enhancement.

## 4. Perturbation transfer check

Manuscript B §16 already shows an **adiabatic transfer architecture** into *radiation*:

`zeta_A = -Psi - H delta_rho_A / dot(rho_A)`,  
`dot(zeta)=-H deltaP_nad/(rho+P)+O(k²/(a²H²))`,  
`delta J = dot(J) delta t` for a shared local clock.

When `deltaP_nad=0` and gradients are negligible, `dot(zeta)≈0`. This preserves a preexisting large-scale structural curvature perturbation into the thermal radiation bath: it is **not** an independent enhancement of CDM halo-scale power. CDM/radiation entropy `S_cr=3(zeta_c-zeta_r)`, the radiation pressure gradient source and subhorizon CDM transfer function require the coupled Einstein–Boltzmann/mapper perturbation equations.

A homogeneous background mapper `M(x)` by itself does not specify `delta M(x,k)`, its kinetic normalization, interaction vertex, or its coupling to Jordan-frame cold matter. Setting `delta M` equal to a background fraction or multiplying a power spectrum by `1/M` would be a new, unsupported physical assumption.

The original post-lock spatial trace mode satisfies the manuscript's heavy-envelope `|u|∝a^-3/2` (when applicable). A free mismatch left from recombination decays to about 0.0026705 by z=20, and from equality to about 0.0004714. A persistent halo-scale enhancement therefore requires a transfer to a genuinely conserved CDM/curvature degree of freedom *before lock* or a new local source after lock. The manuscript has shown adiabatic radiation capture, but **not the requested independent radiation→CDM amplification**.

## 5. What happens to the frozen late300 action F during radiation?

The accepted profile `F=exp(A_F U)`, `U=[1+exp(-(x-x_c)/width)]^-1`, `x_c=-ln(1+z_c)`, `A_F=0.020481464901`, `z_c=3.9278764`, `width=0.33114134`.

- At z=1090, `F-1≈1.6955e-9`, `alpha_M≈5.1201e-9`.
- At z=3466, `F-1≈5.1634e-11`, `alpha_M≈1.5593e-10`.

The specific late **Planck-mass transition** is essentially off at those epochs, so it does not supply a large hidden radiation-era *Planck-mass* activation. **This is not a bound on the early tracker scalar itself**, whose fractional background contribution can be ~1%; F and rho_X are different action effects.

## 6. Final closure decision

**Derived/verified:** `M` relative-to-coasting mapping, `p_s=1+dlnM/dx`, density-route `chi_act`, the exact activation-derivative and transfer *identities*, zero post-handoff J on the ideal radiation tracker, and negligible late300 F excursion at recombination/equality.

**Not determined from M,p,N alone:** an independent `J` beyond reconstruction, radiation-pressure sourcing of `delta M`, a CDM-relative entropy/transfer coefficient `T_c(k)`, and the nonlinear early galaxy/SMBH halo enhancement. Those require an off-shell interaction, specified kinetic/gradient sector, and coupled radiation+CDM+spatial-scalar perturbations.

**Scientific verdict:** preserve the mapper/coasting ratio as a valuable exact geometric and derivative diagnostic, but do not equate it with suppression of physical gravity or treat it as a substitute for radiation-to-matter coupling. No new candidate is promoted; late300 remains unchanged. The *next real microscopic target*, if pursued, is a variationally derived radiation-to-CDM perturbation transfer from a single Jordan-compatible action, with explicit BBN/CMB/isocurvature gates. 
