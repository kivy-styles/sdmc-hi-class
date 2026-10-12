# SDMC Late300: a deeper canonical–Horndeski transition, kinetic/gravity identity, and off-shell No-Slip test

**Date:** 12 October 2026  
**Scope:** Original analytic deductions, independently evaluated numeric tests, and a falsifiable microphysical closure strategy continuing the structural-parameter-independence audit.  
**Benchmark:** `sdmc/late300_candidate.json` (the accepted **late300/g022-A034**, unchanged).  
**Scientific status:** Verified necessary identities and counterexamples. **Not** a first-principles autonomous unified SDMC solution, complete nonlinear perturbation evolution, refitted candidate, or new observational likelihood.

## Executive conclusion

There is a genuine exact on-target connection between the accepted finite kinetic response `D(a)` and gravitational factor `F(a)`: their selected transition function is common. Under the accepted No-Slip constraint, the full Horndeski `alpha_K` law fixes an **additional necessary kinetic operator combination** `K_cov` as a function of redshift. This produces a quantitative microphysical target for deriving `D0` rather than just treating it as a fitted correction. Two serious obstacles are also demonstrated: (1) a two-exponential **adiabatically canonical** potential flattens the tracker slope, **increasing** the tracker density where Late300 prescribes a 50% early-component shutoff; and (2) exact No-Slip at **every** off-shell kinetic `X` with `F_phi != 0` is impossible for a regular `G3_X` at `X = 0`. The accepted reconstructed linear `G3=g(phi)X` is No-Slip **on its selected trajectory**, and is not invalidated by this second fact.

**No new matter generation is assumed or introduced.** Matter and radiation remain separately conserved in the Jordan frame. This work leaves the accepted Late300 and the optimized ΛCDM local021 models untouched and grants neither an automatic zero parameter-count exemption. It derives no new unique effective `Delta k` for model selection.

## 1. Starting facts and field conventions

For the *finite accepted* Late300 target, use `n=ln a` as the cosmological e-fold variable. In the direct structural manuscript, `N_lapse` may instead mean the structural lapse; never confuse the two. Let `sigma=S/S0` be the direct structural field, and let the fitted action reconstruction use a monotonic auxiliary field `phi=n` **on the target background only**. `phi` is not automatically `ln(S/S0)`. Any change from `phi` to `sigma` must transform kinetic couplings.

Accepted parameters:

| Input | Accepted value |
|---|---:|
| `A_F` | `0.02048146490100771` |
| `z_c` | `3.927876388467848` |
| `w_F=w_D` | `0.33114133956842123` |
| `D_floor` | `0.05181542627513409` |
| `D0` | `0.34231919445927034` |
| `p_D` | `1` |
| `lambda_e` | `18.40625` |
| `z_t` | `16.19317622240633` |

Define `n_c=-ln(1+z_c)` and the accepted shared sigmoid

    U(n) = {1+exp[-(n-n_c)/w]}^-1.
    F(n) = exp(A_F U(n)); D(n) = D_floor + D0 U(n).

Note this `U(n)` represents an **already selected profile**, not a dynamical derivation of the transition. The separate derived-*F* force-centre branch has **different** `F` timing and leaves the old `D` window, so the exact shared-U identity below is restricted to the canonical *accepted Late300* target.

## 2. Exact finite gravitational/kinetic identity and microscopic target

Eliminate `U`:

    U = ln(F)/A_F,
    D(F) = D_floor + C ln(F),
    C = D0/A_F = 16.713608919761782,
    D_n = C alpha_M; alpha_M=F_n/F.

The numerical identity was tested on 1,601 redshift points through `z=10^7`, with max absolute residual **2.03e-15** and max derivative residual **5.55e-17**. This is an exact *target-level* identity, not a new independent estimate of `C`.

The restricted covariant class is

    G2(phi,X) = k1(phi) X + k2(phi) X² - V(phi)
    G3(phi,X) = g(phi) X
    G4(phi) = M² F(phi)/2
    G5 = 0,      X = phidot²/2,

in the convention `L3=-G3 Box(phi)`. For arbitrary time-dependent field `phi`, the EFT kinetic coefficient of this class is

    alpha_K = [2X/(H² M² F)] [k1 + 6 k2 X - 4 g_phi X]
              + [12 phidot X g/(H M² F)].

The No-Slip condition `alpha_B = -2 alpha_M`, with

    alpha_M = phidot F_phi/(H F)
    alpha_B = [2 phidot/(H M² F)] [X g - M² F_phi/2],

implies **on the background solution**

    X g = -M² F_phi/2.

Substitution removes the explicit braiding contribution from `alpha_K`:

    12 phidot Xg/(H M² F) = -6 alpha_M.

Let

    K_cov := 2 X [k1 + 6 k2 X - 4 g_phi X]/(H² M² F).

Since `D=alpha_K+3 alpha_B²/2=alpha_K+6 alpha_M²`, a new exact **necessary on-shell identity** follows:

    **D = K_cov - 6 alpha_M + 6 alpha_M²**
    **K_cov = D + 6 alpha_M - 6 alpha_M².**

All factors refer to the same restricted Horndeski action convention. Unlike the equality `D(F)` alone, this specifies a *particular combination of microscopic kinetic coefficients* that an autonomously derived action must reproduce.

### Required kinetic-composite test sheet

| z | U | alpha_M | D | K_cov required |
|---:|---:|---:|---:|---:|
| 1090 | `0.000000083` | `0.000000005` | `0.051815` | `0.051815` |
| 30 | `0.003858` | `0.000238` | `0.053136` | `0.054562` |
| 16.19318 | `0.022453` | `0.001358` | `0.059501` | `0.067636` |
| 10 | `0.081293` | `0.004619` | `0.079644` | `0.107231` |
| 3.92788 | `0.500000` | `0.015463` | `0.222975` | **`0.314317`** |
| 2 | `0.817387` | `0.009232` | `0.331623` | `0.386505` |
| 1 | `0.938381` | `0.003576` | `0.373041` | `0.394423` |
| 0 | `0.991969` | `0.000493` | `0.391385` | `0.394340` |

At the central `z_c`, the direct term `6 alpha_M² = 0.001434586` is only **0.838%** of the `D0/2 = 0.171160` transition correction. So `D0` is **not** merely the braiding-squared contribution. The action-level composite `K_cov` grows significantly and must be explained through `k1`, `k2` and the field derivative `g_phi`. The present reconstruction determines these coefficients from an already selected `D(n)`; a true origin must predict them independently.

Early canonical radiation tracking with fixed `F` and no braiding gives `D_can,r=16/lambda_e²=0.0472268903`, versus finite accepted `D_floor=0.0518154263`; difference **0.004588536**. In the mature canonical infrared branch, a separately established action gives `D_IR=2`, whereas blindly extending the finite logistic gives only `D_floor+D0=0.394134621`. A future kinetic handoff remains necessary; it is not selected by the existing finite interpolation.

## 3. One scalar, two epochs: a two-exponential necessary-condition test

Consider, only as an explicit candidate to falsify,

    V(phi)=V_e exp(-lambda_e phi)+V_IR exp(-lambda_IR phi)
    lambda_e=18.40625, lambda_IR=sqrt(2), Delta_lambda=16.9920364376.

Define the fractional shallow-potential weight

    U_V(phi) = [1+(V_e/V_IR) exp(-Delta_lambda phi)]^-1.

Then **exactly**

    dU_V/dn = Delta_lambda phi_n U_V(1-U_V).

If we identify `U_V` with the accepted `U` locally, its transition width requires

    phi_n,required = 1/[Delta_lambda * w_F] = 0.1777219902.

A canonical scalar on the matter-era scaling solution instead has

    phi_n,tracker = 3/lambda_e = 0.1629881154
    w_predicted = lambda_e/[3(Delta_lambda)] = 0.3610760069.

The predicted width is **9.03985% wider** than the accepted `0.33114134`, and the required field speed is **9.03985% larger** than the canonical matter-tracker value. This is a conditional candidate, not a proof that the full kinetic/braided field cannot give the needed 9% departure. The field coordinate here is canonically normalized; equating it to `ln(S/S0)` would be incorrect without the field redefinition. With the manuscript's matter-era `p=dlnS/dn≈1.5`, a match would require `dphi/dlnS≈0.1184813` near transition.

The **two epochs are not the same**:

    n_t = -2.84451257385; n_c=-1.59490814241;
    n_c-n_t=1.24960443144 e-folds.

The accepted shared sigmoid has

    U(z_t)=0.0224528580,  U(z_c)=0.5, U(today)=0.9919691128.

If `U` is identified with `U_V`, the steep exponential remains **97.7547%** of potential at the early tracker-shutoff epoch. Therefore equal-potential crossover at `z_c` is not an explanation of the tracker shutoff `z_t`. The potential-amplitude ratio and field initial state must select the crossover epoch; the two slopes by themselves do **not** select `z_c`.

### Stronger sign/size test of the early handoff

The canonical matter-era scaling fraction is `f(phi)≈3/lambda_eff²`, where the two-exponential potential has

    lambda_eff = (1-U)lambda_e + U lambda_IR.

At `z_t`, `lambda_eff = 18.0247302183`, and the instantaneous adiabatic scaling approximation predicts

    f_steep = 3/lambda_e² = 0.00885504193;
    f_two_slope(z_t) = 3/lambda_eff² = 0.00923386899,
    a **4.2781% increase**.

Late300's independent early tracker window is **halfway shut off** at `z_t`; its matter-era reference has `f_Late300(z_t)≈0.00442752` (with the actual small radiation contribution, `0.00443480`). The direction is opposite: lowering the canonical potential slope increases the adiabatic tracking fraction instead of halving it.

**Scientific interpretation:** Ordinary *adiabatic* two-exponential canonical tracking cannot by itself reproduce the early component's shutdown. A viable common mechanism requires exiting the adiabatic tracker, a kinetic/braiding transfer to a different scalar behaviour, or a derived structural activation rule; it does **not** require creating new matter. `Gamma-1 = Delta_lambda² U(1-U)/lambda_eff² ≈ 0.019506` at `z_t`, giving an independent measure of how modest the potential-shape departure still is.

## 4. Exact No-Slip for a fixed microscopic action: an off-shell regularity theorem

For the restricted luminal Horndeski class with `G4=M² F(phi)/2` and `G5=0`, exact No-Slip means

    X G3_X = - (M²/2) F_phi.

If `F_phi != 0` and this equality must hold for *every kinetic state* `X>0` at the same `phi`, integrate with respect to `X`:

    G3(phi,X)=-(M²/2) F_phi(phi) ln(X/Xstar)+h(phi).

Thus, for a genuinely evolving gravitational coupling, **exact off-shell No-Slip for all `X>0` forces logarithmic dependence on `X`**. A `G3` analytic with bounded `G3_X` at `X=0` instead has `X G3_X -> 0`, while the right-hand side remains nonzero. The two conditions are incompatible at `X=0` unless `F_phi=0`. This is a simple mathematical regularity obstruction, **not a proof that the accepted on-shell No-Slip action fails**.

The accepted *linear* `G3=g(phi) X` avoids the off-shell requirement by selecting

    g(phi)=-M² F_phi(phi)/[2 X_bg(phi)]

along its chosen background. On an alternative trajectory at the same `phi`, its leading No-Slip residual is

    alpha_B+2alpha_M = alpha_M [1 - X/X_bg(phi)].

For example, changing `phidot` by +10% at fixed `phi,H,F` changes `X/Xbg` to 1.21 and produces residual `-0.21 alpha_M`; changing it by -10% gives `+0.19 alpha_M`. These are analytic instantaneous deviations, **not** results of re-evolving the original tested basin (whose trajectories can approach the accepted one at late times).

One possible *regularized* logarithmic demonstration is

    G3(phi,X) = -(M²/2)F_phi ln(1+X/Xstar) + h(phi),

for which

    (alpha_B+2alpha_M)/alpha_M = Xstar/(X+Xstar).

The residual can be small when `X >> Xstar` but is not identically zero; the scale `Xstar` is **new microphysics**, not a derived SDMC constant. The singular logarithmic option may also require careful `X->0` dynamics. Neither is proposed as a replacement for the accepted candidate without a full scalar stability and perturbation calculation.

**Physical consequence:** a valid microphysical unification cannot demand simultaneously (a) exact No-Slip for arbitrary field velocities, (b) analytic `G3_X` at `X=0`, and (c) nonzero `F_phi`. It must specify an operating branch, relaxing one assumption or employing a complete off-shell completion consistent with physical stability.

## 5. What a successful *parameter-reducing* action must now do

The following conditions are concrete and independently testable:

1. Specify `F(sigma)`, `k1(sigma)`, `k2(sigma)`, `g(sigma)`, `V(sigma)` and microscopic boundary conditions *before looking up* accepted `H(n),F(n),D(n)`; retain direct `sigma=S/S0` and calculate the field transformation to the auxiliary `phi` where required.
2. Solve Friedmann, Raychaudhuri and `dot J+3HJ=P_sigma` autonomously with separately conserved matter and radiation. Do **not** insert the fitted `z_t`, `z_c`, `A_late`, `B_late` or the accepted sigmoids as external forcing functions.
3. Derive the early tracker shutoff and later gravitational transition as **distinct dynamical observables**; a canonical two-exponential crossover alone fails the adiabatic sign/size test.
4. Obtain the noncanonical kinetic composite `K_cov = 2X(k1+6k2X-4g_phi X)/(M² F H²)` and compare it with the test sheet's required `D+6alpha_M-6alpha_M²`. Matching the table by reconstructing `k1,k2` *from the table itself* would be circular.
5. Derive the scalar trajectory and a genuinely physical No-Slip/near-No-Slip law, testing it for perturbed kinetic trajectories and `X->0` regularity; satisfy `F>0`, scalar no-ghost, positive `cs²`, tensor luminality, causality and local-gravity tests.
6. Derive the selected ratio `C=D0/A_F`, not merely insert its accepted value. Ensure the later canonical IR limit `D->2` is achieved by a physical future operator flow, not by extrapolating the finite sigmoid to `0.394`.
7. Then repeat Planck, DESI, SN, BBN, CMB lensing, shear and SH0ES tests against optimized local021. Count as independently fitted only the initial/action constants actually selected in a prospectively specified family. A frozen-run numerical parameter count of zero applies to *both* local021 and Late300, not to either full family.

**Decision gate:** A coherent fixed action that predicts the numerical late300 coefficients up to observational error, without target reconstruction, would legitimately reduce the count of independent SDMC-specific cosmological parameters. This report does not yet establish such a fixed action; it supplies quantitative **necessary conditions and a rejected simplistic mechanism**.

## 6. Source and reproducibility ledger

- Accepted Late300 values and results: [`sdmc/late300_candidate.json`](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/late300_candidate.json).
- Reconstruction ordering and definitions: [`sdmc/run_late300_linear_covariant.py`](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/run_late300_linear_covariant.py).
- Existing full Horndeski scalar equations and activation/handoff audit: [10 October field-equation report](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/reports/activation_handoff_field_equations_2026-10-10.md).
- Earlier restricted autonomous two-exponential test (does not match Late300): [autonomous countermodel report](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-activation-handoff-equations-audit-20261010/sdmc/reports/activation_autonomous_countermodel_results_2026-10-10.md).
- Mature action and kinetic limits: [D-sector final closure audit](https://github.com/kivy-styles/sdmc-hi-class/blob/sdmc-D-forward-canonical-attractor/sdmc/D_all_remaining_derivations_final_closure.md).
- Independent mathematical background/action jet nonuniqueness: [GitHub run 37833138647](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/37833138647).
- Independent physical No-Slip condition: [Linder, *No Slip Gravity*, 2018](https://arxiv.org/abs/1801.01503).
- Horndeski alpha identities: [hi_class paper](https://escholarship.org/uc/item/0qb774wk), Appendix A.32–A.34; Bellini & Sawicki 2014.
- Reproducible script: `sdmc/run_unified_transition_audit.py`; regression checks: `sdmc/test_unified_transition_audit.py`; plots: `sdmc/plot_unified_transition_audit.py`.

This audit supplements, **does not override**, the earlier frozen-action robustness and independent model-family parameter audits. The main reported new results are necessary mathematical constraints, and all accepted benchmark likelihoods remain unchanged.