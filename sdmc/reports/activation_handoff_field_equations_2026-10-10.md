# SDMC: Covariant activation and handoff from field equations — derivation and sign audit

**Date:** 2026-10-10  
**Branch status:** Independent theoretical diagnostic; no modification of accepted late300, no new cosmological likelihood claim.  
**Basis:** SDMC Manuscript B, Route A §§3–5 and 13; Route B §20; final retrospective §6; \`sdmc/late300_candidate.json\`, \`sdmc/apply_sdmc_full_background_patch.py\`, \`sdmc/run_late300_linear_covariant.py\`. Horndeski FLRW equations: Kobayashi, Yamaguchi & Yokoyama, *Generalized G-inflation* (2011), doi:10.1143/PTP.126.511.

## 1. Field definition and action

Use the Jordan-frame matter metric and a spatially flat background. Write \`N=ln a\`, \`sigma=S/S0\`, \`p=d ln S/dN=sigma_N/sigma\`, \`X=dot(sigma)^2/2\`, and Ricci scalar \`calR=6(2H^2+dot H)\`. **calR is curvature, not structural radius R_str=l_P S.** With \`M^2=M_Pl^2\` as a constant reference normalization, take the accepted *class* of action functions

\`\`\`
G2(sigma,X)=k1(sigma) X + k2(sigma) X^2 - V(sigma)
G3(sigma,X)=g(sigma) X
G4(sigma)=M^2 F(sigma)/2
G5=0
L3=-G3 Box(sigma).
\`\`\`

These are the restricted linear-G3 Horndeski conventions; all identities below require this exact convention. In the repository's on-shell reconstruction, \`phi=N\` is used as a convenient monotonic field coordinate. Identifying that field with the direct manuscript field \`sigma=S/S0\` requires a field redefinition and correct transformation of all kinetic/action coefficients. The two cannot silently be equated.

## 2. Exact flat-FLRW background identities

For separately conserved Jordan-frame matter and radiation,

\`\`\`
rho_m,N=-3 rho_m ; rho_r,N=-4 rho_r.
rho_m=rho_m0 exp(-3N); rho_r=rho_r0 exp(-4N).
dot(sigma)=H sigma_N=H p sigma ; X=(H p sigma)^2/2.
\`\`\`

The Hamiltonian equation reduces to

\`\`\`
3 M^2 F H^2 = rho_m+rho_r + rho_X^J

rho_X^J = k1 X + 3 k2 X^2 + V
          + 6 H g X dot(sigma) - 2 g_sigma X^2
          - 3 M^2 H dot(F).
\`\`\`

The acceleration equation is

\`\`\`
-M^2 F (2 dot(H)+3H^2) = rho_r/3 + p_X^J

p_X^J = k1 X+k2 X^2-V
        -2X [g_sigma X+g ddot(sigma)]
        +M^2[ddot(F)+2H dot(F)].
\`\`\`

The homogeneous scalar equation is

\`\`\`
dot(J)+3H J=P_sigma,

J = dot(sigma) [k1+2k2 X -2 g_sigma X] +6H g X,

P_sigma = k1_sigma X+k2_sigma X^2-V_sigma
          -2X[g_sigmasigma X+g_sigma ddot(sigma)]
          +(M^2/2) F_sigma calR.
\`\`\`

The Friedmann constraint, scalar evolution, and matter/radiation conservation form a coupled differential-algebraic problem once \`k1,k2,g,V,F\` are **independently supplied**. Reconstructing these functions from a target \`H(N)\` does not independently predict that history.

## 3. Activation from the action, rather than an imposed chi(N)

Manuscript B defines \`rho_v=rho_P S^-2 = rho_v0 sigma^-2\` and schematically \`rho_X=chi rho_v\`. In the action's Jordan-frame Friedmann split, define

\`\`\`
chi_J(N) := rho_X^J(N)/rho_v(N)
          = sigma^2 rho_X^J(N)/rho_v0.
\`\`\`

Thus the action supplies, *on any solved trajectory*, the effective activation:

\`\`\`
chi_J = sigma^2/rho_v0 [
    k1 X + 3k2 X^2 + V +6H g X dot(sigma)
    -2g_sigma X^2 -3M^2 H^2 F_N
].
\`\`\`

This is **not an independent evolution equation** for \`chi\`; it is a derived diagnostic once the true action and trajectory are known. It need not remain nonnegative. Also \`chi_J\` is not the same as the constant-G geometric residual

\`\`\`
chi_geom := [3 M^2 H^2-rho_m-rho_r]/rho_v,
chi_geom = chi_J/F +(1/F-1)(rho_m+rho_r)/rho_v.
\`\`\`

The latter can be negative even if the former is positive (and vice versa). Do not use \`ln chi\` across zeros or identify either with an observable vacuum density without specifying the gravitational frame and split.

## 4. Exact structural-vacuum handoff identity

The positive-density structural-to-matter ratio is \`Q=rho_X^J/rho_m=chi_J rho_v/rho_m\`. For intervals with \`chi_J>0\`, define \`gamma_J=d ln chi_J/dN\`. Then

\`\`\`
d ln Q/dN = gamma_J +3-2p.
Q(N)=Q(N_i) exp[ integral_(N_i)^N (gamma_J+3-2p) dN ].
\`\`\`

The non-log form valid even at \`chi_J=0\` is

\`\`\`
Q_N = Q(3-2p)+(rho_v/rho_m) (chi_J)_N.
\`\`\`

A structurally dominated handoff with \`Q=1\` is dynamically **predicted** only when this condition is encountered during autonomous evolution, not when a preset \`z_t\` or sigmoid crosses a chosen value. The onset of relative growth is \`gamma_J>2p-3\` where \`chi_J>0\`.

For pure matter tracking \`p=3/2\`, a constant \`chi_J\` leaves \`Q\` constant. For \`p<3/2\`, the raw structural density gains relative to matter even if \`chi_J\` is constant. That fact alone does not guarantee acceleration: \`rho_v∝a^-2\` in a \`p=1\` regime behaves like \`w=-1/3\` *if separately conserved*; additional action-dependent activation/dynamics are needed to obtain observed acceleration.

This late structural handoff is distinct from (i) radiation–matter equality near \`z~3466\`, (ii) the accepted early tracker shutoff at \`z_t=16.1932\`, and (iii) the selected Planck-mass \`F\` transition centre at \`z_c=3.9279\`. These epochs must not be conflated.

## 5. Gravitational drift and density relation

Under the manuscript's unscreened background identification \`G_eff ∝ F^-1\`,

\`\`\`
alpha_M=F_N/F;  dot(G_eff)/G_eff = -alpha_M H.
d ln G_eff / d ln rho_m = alpha_M/3.
\`\`\`

The manuscript's earlier accepted candidate reported \`-7.7e-14 per year\`, **not** the canonical late300 value. For canonical late300, the published F profile has \`A_F=0.02048146490100771\`, \`z_c=3.927876388467848\`, \`width=0.33114133956842123\`, and

\`\`\`
U=[1+exp(-(N-Nc)/width)]^-1; Nc=-ln(1+zc)
F=exp(A_F U), alpha_M=(A_F/width) U(1-U).
\`\`\`

At \`z=0\`: \`F=1.0205247753\`, \`alpha_M=0.0004927303\`, and with \`H0=69.71482083 km/s/Mpc\`, \`dot(G_eff)/G_eff=-3.51308e-14 per year\`.

At \`z=z_c\`: \`F=1.0102933482\`, \`alpha_M=0.0154627816\`, \`F_N=0.0156219454\`, and \`rho_m(z_c)/rho_m0=(1+z_c)^3=119.66838\`. The drift law is a *diagnostic of a selected F history*, not a unique field-law inversion.

## 6. Exact conditional curvature-minimum sign test

The manuscript's simplified slow-field curvature-force balance is

\`\`\`
U_sigma = V_sigma -(M^2/2) F_sigma calR = 0.
m_eff^2 = V_sigmasigma -(M^2/2) F_sigmasigma calR.
\`\`\`

Differentiating the minimum condition along \`N\`, holding \`V,F\` fixed functions of the field, yields

\`\`\`
F_N = (M^2/2) F_sigma^2 calR_N /m_eff^2.
\`\`\`

For a **stable adiabatic minimum** (\`m_eff²>0\`) when \`calR_N<0\`, this predicts **\`F_N<0\`**. The canonical late300 target has \`F_N>0\`; thus this simple equilibrium-tracking mechanism cannot produce that target.

Using the *published late300 target background formula* and \`calR=6(2H²+dot H)\` at \`z_c=3.927876\`:

\`\`\`
calR/H0²=109.91689,
d ln calR/dN=-2.70720,
calR_N/H0²=-297.56708,
F_N=+0.0156219454,
[m_eff²/(M² F_sigma²)]/H0²
   =calR_N/(2 F_N)
   ≈ -9524.0.
\`\`\`

The last line is the **negative inferred equilibrium curvature** and contradicts positive \`m_eff²\` for the simplified stable-minimum interpretation. It is **not** a measured physical tachyon mass, and it does not by itself invalidate the full kinetic/braided Horndeski action.

Generalized force balance includes a kinetic/braiding force \`K_sigma(sigma,N)\`. In an adiabatic generalized equilibrium,

\`\`\`
V_sigma -(M²/2)F_sigma calR+K_sigma(sigma,N)=0
m_eff²=V_sigmasigma -(M²/2)F_sigmasigma calR +partial_sigma K_sigma
F_N=F_sigma[(M²/2)F_sigma calR_N-partial_N K_sigma]/m_eff².
\`\`\`

With \`m_eff²>0\`, \`calR_N<0\`, and late300's \`F_N>0\`, the explicit time-dependent kinetic/braiding response must overcome the curvature-tracking term, **or** the field must leave adiabatic equilibrium. A delayed dynamical release/overshoot is another conditional possibility.

For the complete action the scalar source is not just \`U_sigma\`: it is the exact \`dot J+3HJ=P_sigma\` above. Solve that equation instead of fixing a desired \`F(N)\` by an external logistic gate.

## 7. Dynamical handoff algorithm and falsification criteria

1. Specify one microphysically motivated, covariant \`k1(sigma), k2(sigma), g(sigma), V(sigma), F(sigma)\` independent of the accepted target \`H(z)\`, \`F(N)\`, \`z_t\`, and late \`delta_H\`.
2. Set physically justified early initial \`sigma\`, \`dot sigma\`, \`rho_m\`, \`rho_r\`; solve the Hamiltonian constraint for positive \`H\`. No explicit matter-scalar energy exchange is required in the Jordan frame.
3. Integrate the exact Friedmann+Raychaudhuri+scalar system above. Monitor constraint conservation; identify the actual \`p=\dot sigma/(H sigma)\`.
4. Derive \`F(N)\`, \`alpha_M\`, \`rho_X^J\`, \`chi_J\`, \`Q\`, and **emergent** crossing times; do not insert the two-lobe deformation or a preselected handoff switch.
5. Track \`calR_N\`, \`m_eff²\` where meaningful, and signs of kinetic/braiding contributions. Compare against the stable-minimum sign audit.
6. Require \`F>0\`, \`Q_s>0\`, \`c_s²>0\`, acceptable No-Slip/growth, and recovered early BBN/CMB/BAO observables and calibrated local ladder. Fully recompute Planck and DESI only after physical background/perturbation checks pass.
7. Compare with baseline late300, leaving the benchmark unchanged.

**Scientific closure:** The governing equations and the activation/handoff identities are derived here. A *unique autonomous transition solution* is **not** derived from the currently supplied data: the current reconstruction was generated from a selected target and the microphysical coefficient/initial-state law remains open. No new observed likelihood or SH0ES-closure claim follows.

**Numerical reproducibility:** See the companion \`sdmc/audit_activation_handoff_sign.py\` on this branch for the branch-prescription curvature/F drift audit. 
