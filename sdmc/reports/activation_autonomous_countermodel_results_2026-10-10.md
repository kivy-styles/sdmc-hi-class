# SDMC autonomous activation and handoff: verified forward experiment

**Date:** 2026-10-10  
**GitHub Actions verified success:** https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38042838049  
**Script:** [../run_activation_autonomous_countermodel.py](../run_activation_autonomous_countermodel.py)  
**Parent derivation:** [activation_handoff_field_equations_2026-10-10.md](activation_handoff_field_equations_2026-10-10.md)

## Question

Can a fixed scalar-tensor action evolve forward with ordinary matter/radiation conservation and generate the late handoff **without** imposing an observed or fitted transition redshift? This is a restricted countermodel and **not** a forward solution of the accepted late300 action.

## Fixed analytical action class

Jordan-frame reduced units, `8pi G_reference=1`:

```
G2 = X - V(phi)
G3 = 0
G4 = F(phi)/2
G5 = 0

F(phi) = exp(beta phi)
V(phi) = V_steep exp(-lambda_e phi) + V_IR exp(-sqrt(2) phi)
lambda_e = 18.40625
V_IR = 3 * 0.7014 = 2.1042, V_steep/V_IR in {0.01, 1}
beta in {0, +0.02048146490100771, -0.02048146490100771}
initial N = -12, phi_ini = -0.5, phi_N_ini = 0
rho_m(N) = 3*0.2985059603044 exp(-3N)
rho_r(N) = 3*[4.17998772e-5/(0.6971482083085)^2] exp(-4N)
```

The slopes and beta scale were motivated by manuscript/candidate values; potential amplitudes, initial field state and reference density normalizations remain **chosen**, not first-principles derived. No `H(z)`, `F(N)`, sigmoid, `z_t`, `z_c`, or handoff/activation function appears in the evolution equations.

## Exact autonomous equations integrated

With `y=dphi/dN`, `u=dlnH/dN`, `rho_m` and `rho_r` separately conserved, the canonical scalar equations are:

```
H² = [rho_m+rho_r+V]/[3F+3 F_phi y-y²/2]

y_N + (y-3F_phi) u = -3y -V_phi/H² +6F_phi

F_phi y_N +(2F+F_phi y)u
  = -(rho_m+4rho_r/3)/H² -(1+F_phiphi)y² +F_phi y
```

For `G3=0`, these are the scalar equation, Raychaudhuri equation, and exact Friedmann constraint under the model assumptions. Algebraic elimination gives the derivatives, which are integrated numerically. The independently checked maximum discrepancy between the numerical `dlnH²/dN` and `2u` was **3.3703e-9** across the sampled trajectories. This establishes internal background accuracy, not observational or perturbative stability.

## Verified six-case outcome

All redshifts are *outputs*, not inputs. The action-split equality `rho_X^J/rho_m=1` differs from the kinematic acceleration threshold `q=0`.

| beta | V_steep/V_IR | z(q=0) | z(rho_X^J/rho_m=1) | z(V_steep=V_IR) | q(today) | alpha_M(today) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.01 | 0.417136 | 0.223277 | 7.948426 | -0.156005 | 0 |
| 0 | 1 | 0.245488 | 0.079348 | 7.216594 | -0.110116 | 0 |
| +0.020481465 | 0.01 | 0.357679 | 0.197345 | 7.412945 | -0.123408 | +0.0153878 |
| +0.020481465 | 1 | 0.201933 | 0.059233 | 7.874110 | -0.083944 | +0.0136699 |
| -0.020481465 | 0.01 | 0.464027 | 0.249211 | 8.468704 | -0.184571 | -0.0142405 |
| -0.020481465 | 1 | 0.283161 | 0.097169 | 6.634331 | -0.135426 | -0.0125229 |

Representative `beta=+0.020481465`, ratio `0.01` has `F(today)=1.0087141`, `alpha_M(today)=0.0153878`, `q(today)=-0.123408`, `Q_J(today)=1.484212`, `H(today)/H_ref=sqrt(0.7352313)=0.85746`.

The accepted **late300 target** instead has `F(today)≈1.020525`, `alpha_M(today)≈0.0004927`, `q(today)≈-0.527`; its background-defined acceleration and action-split equality occur at `z≈0.707` and `z≈0.349`. The autonomous countermodel clearly does **not** reproduce late300.

## What was learned and what was not

1. A fixed action and initial state can **autonomously generate** a change of potential dominance, a deceleration/acceleration crossover, and equality of action-defined structural and matter contributions without a prescribed transition time.
2. **Nonuniqueness of timing is explicit:** varying the independent ratio of potential amplitudes and beta significantly moves the epochs. The late300 transition is not selected automatically by the general field equations.
3. The **sign of F running is dynamically possible** under nonadiabatic scalar rolling even while curvature falls. The simple stable moving-minimum sign obstruction does *not* apply to a rolling field out of equilibrium. This demonstrates logical viability, not realization by the accepted SDMC coefficients.
4. A potential made from slopes `lambda_e=18.40625` and `lambda_IR=sqrt(2)` is **not** a derived microscopic interpolating potential merely because those asymptotic slopes were borrowed from the SDMC program.
5. A nonminimal canonical scalar with `G3=0` and running `F` generally satisfies `alpha_B=-alpha_M`, **not** the accepted No-Slip target `alpha_B=-2alpha_M`. Therefore this is not a valid *replacement* for the accepted linear-G3 SDMC action.
6. `phi` in this canonical action is **not established to equal** the manuscript's direct structural clock `ln(S/S0)`; the matter-era structural target `p≈3/2` is not reproduced. Identifying the fields would require a justified field redefinition and consistent transformation of all coefficients.
7. Neither `H0` nor the absolute potential scale is predicted: `V_IR` is normalized by a reference density. There has been **no** Planck, DESI, SH0ES, BBN, growth, or Bayesian-evidence test of this countermodel.

## Stronger analytic consequence: No-Slip requires a kinetic-braiding law

For the restricted Horndeski class with `G4=M²F(phi)/2`, `G5=0`, and no `G4_X`, the exact linear-perturbation background identity is

```
alpha_M = dot(phi) F_phi /(H F)

alpha_B = [2 dot(phi)/(H M²F)] [X G3_X - (M²/2) F_phi].
```

Demanding the accepted **No-Slip branch** `alpha_B=-2 alpha_M` yields

```
X G3_X = -(M²/2) F_phi.
```

The canonical autonomous test has `G3=0` and therefore `alpha_B=-alpha_M` whenever `F_phi !=0`: the toy **fails the exact accepted No-Slip identity**, even though it passes the background Bianchi check. This is a precise physical reason not to promote it.

For the accepted restricted *linear-G3* action `G3=g(phi)X`, exact No-Slip **along a given trajectory** requires

```
g(phi) = - M² F_phi(phi)/[2 X_bg(phi)].
```

But the function `X_bg(phi)` must itself emerge from the actual field equations. Defining `g(phi)` from a preselected `X_bg(phi)` is a background reconstruction, not an independent derivation. One possible off-shell closure for `X>0` is `G3=-(M²/2)F_phi ln(X/Xstar)+C(phi)`, which satisfies the No-Slip identity algebraically, but its `X->0` behavior and kinetic stability must be tested; it is not identical to the accepted linear-G3 action.

## A second obstacle: accumulated F excursion versus tiny present drift

Canonical late300 has `ln F(today)≈0.020316` while `alpha_M(today)≈0.000493`, a ratio of order **41**. The representative autonomous exponential coupling instead gives `ln F(today)≈0.008676` and `alpha_M(today)≈0.015388`, ratio about **0.56**. In a pure `F=exp(beta phi)` law, `alpha_M=beta phi_N`, so any substantial ongoing field motion keeps the Planck-mass running appreciable. Reproducing the late300 combination of accumulated `F` excursion and tiny current drift therefore requires the field motion to slow sufficiently, or `F_phi/F` to become much smaller by the present epoch (for example near a derived plateau).

The history-dependent plateau is **not** supplied by merely setting a small present-day `dot G/G`; it must be generated from a fixed action with a physically determined characteristic field scale and initial state.


## Next necessary physical completion

Replace the illustrative `F=exp(beta phi), G3=0` with **independently specified** manuscript-consistent `F(sigma), G3=g(sigma)X, k1(sigma),k2(sigma),V(sigma)` that satisfy No-Slip, positivity, and the radiation/matter/stiff/mature limits. Then solve the *full* Friedmann + scalar + perturbation evolution without referencing an imposed late300 `H(N)` or `F(N)`, and test the emergent late handoff against all observables. In the absence of independent coefficient functions and UV initial-state selection, such a unique prediction cannot currently be claimed.

## Workflow provenance

GitHub workflow success: https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38042838049  
The separate historical `sdmc_reconstructed_covariant_evolution.yml` workflow continues to register validation failures with no jobs on this branch; they are **not** failures of the new autonomous forward script.
