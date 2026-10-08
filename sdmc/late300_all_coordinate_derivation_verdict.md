# SDMC late300 — simultaneous derivation and independence audit of all ten structural coordinates

**Date:** 2026-10-08. **Scope:** The two agreed goals only: derive the remaining structural coordinates from SDMC dynamics; determine whether the resulting SDMC model surpasses the optimized local021 ΛCDM control under comparable observational data and parameter penalties.

**Sources:** Manuscript B Part I, secs. 7–10; Structural Radius / Jacobian Appending Continuation VI, secs. 122–137; \`sdmc/late300_candidate.json\`; \`sdmc/run_late300_full_structural_derivation.py\`; \`sdmc/run_late300_linear_covariant.py\`; completed GitHub runs cited below. This file does not postulate new microscopic values.

## Derivations and independent status

For N=ln(a), define S_F(N)=1/[1+exp(-(N-N_c)/w_F)] and F(N)=exp(A_F S_F).
Then

- alpha_M = d ln F/dN = A_F S_F(1-S_F)/w_F,
- N_c = -ln(1+z_c), S_F(N_c)=1/2,
- max(alpha_M) = A_F/(4w_F), and ln[F(+infinity)/F(-infinity)] = A_F.

These are exact identities *after* the F trajectory has been specified, not independent microscopic equations selecting A_F, N_c, or w_F. The prior force-decomposition root z_c=4.034077502899133 and width w_F=0.32925377073762596 are conditional candidates from the reconstructed action. The repeated-stability scan shows the derived-handoff lower amplitude near **A_F,min=0.020517949990009835**, while A_F=0.020520 is healthy. **Stability gives a one-sided condition; it does not uniquely choose the physical amplitude above threshold.** At the historical handoff z_t=16.19317622240633 the numeric minimum is A_F,min~0.020517355261795354.

For the independent kinetic support D(N)=D_floor+D_0 S_D(N)^{p_D}, with alpha_B=-2alpha_M (No-Slip),

**Authorial correction (8 October 2026):** D_0 was introduced as a **correction to D_floor**, not a second unrelated kinetic baseline or separate physical mechanism. In the present ansatz D_floor is the early-structural base support and D_0 S_D(N)^p_D is a time-dependent corrective increment. Consequently, if **independent** evolution of the same action predicts both asymptotic kinetic support levels, D_0 = D_late - D_floor is a conditional derived difference. For the accepted late300 window, the code yields S_D(0)=0.9919691128, D_correction(0)=0.3395700676 and D(0)=0.3913854939; at infinite late N, D_late=0.3941346207. The asymptotic full late normalization is not yet predicted without using the accepted trajectory. The physical relation reduces the number of **origins** that need explaining, but does not erase a statistically tunable normalization unless the same dynamics fix both endpoints.


D = alpha_K + (3/2) alpha_B^2 = alpha_K + 6 alpha_M^2.

If a *microscopic* action independently provided D(N), the asymptotic values would fix D_floor and D_0, and the logarithmic slope would fix p_D. Here the accepted action coefficients were reconstructed from an already selected target D(N). The linear-G3 inversion in \`sdmc/run_late300_linear_covariant.py\` explicitly gives

k2 = [F alphaK_target-C+4 g_phi X-6 g H^2]/(4X), k1=C-2k2 X,

and alphaK_reconstructed = (k1+6k2X-4g_phi X+6gH^2)/F.

Substitution recovers alphaK_target **for any prescribed D_target** in the admissible algebraic domain; this is an exact inverse-problem degeneracy, not evidence that all such actions are globally stable or observationally viable.

Canonical radiation tracker D_rad,can=16/lambda_e^2=0.047226890271848634 at accepted lambda_e=18.40625. The accepted D_floor=0.05181542627513409 has ratio 1.0971593932370478. Applying the manuscript's *reduced positive-k1,k2 kinetic* mapping (conditional) yields r_K=0.025530092455387544 and c_X^2=0.9114445960760634. This **inverts** the accepted floor; it does not independently derive r_K. The full reconstructed linear-G3 action has k1<0 below z~24.5; its physical k2X/k1 is -0.9037 today, +17.55 near z~100, +1.66 near z~1090, and +0.415 near z~10030. These are not the same universal reduced-kinetic parameter.

For the reduced homogeneous k1(sigma)X+k2(sigma)X^2 sector with negligible V,G3,G4 and transfer, define b_i=d ln k_i/d ln a. Then

d ln(r_K)/d ln(a) = [(1+3r_K)(b2-2b1)-6(1+2r_K)]/(1+6r_K).

If, *additionally*, k2/k1^2 scales as structural S^m and Manuscript B's no-transfer p=d ln S/d ln a=2+1/(1+3r_K) holds, the formula simplifies to

d ln(r_K)/d ln a = 3(1+2r_K)(m-2)/(1+6r_K).

Hence m>2 is necessary for growth **within these assumptions**. No accepted independent coefficient trajectory currently selects m. The reduced prethermal equations must not be substituted for the full late Horndeski equations.

An analytic small-S_D kinetic expansion with nonzero leading linear coefficient chooses p_D=1 as its leading power. If the linear coefficient vanishes, p_D=2 or higher is consistent. Even the positive D ansatz admits p_D=0.5,1,1.5,2 with positive D at representative epochs. Thus minimality does not prove p_D=1.

The tracker identity is f_r=4/lambda_e^2 (f_m=3/lambda_e^2). The accepted lambda_e=18.40625 returns f_r=0.011806722567962159. But f_r was inferred from a background already supplied lambda_e. Alternative positive lambda_e give alternative valid scaling fractions: an **independently predicted** f_r or scalar action slope is required.

The handoff W(z)=.5[1-tanh((N-N_t)/DeltaN_t)] automatically has W(z_t)=1/2 for *any* chosen z_t. An independent root of an action-derived activation-curvature diagnostic supplies the conditional candidate z_t~16.742. The eight helper-dependent root searches in run 37763535074 returned a central-root spread from **16.74133198 to 16.74314404** and multiple other curvature-event roots (~11.3 and ~29.6 at helper 16.7422). These other roots are **not** independently proven fixed points. The current investigation does not establish an eight-decimal unique handoff predicted without the helper ansatz, nor a robust initial-condition selection.

For the late two-lobe correction delta(z)=A_late z exp(-z/tauA)-B_late z^2 exp(-z/tauB), with tauA=.25,tauB=1.5:

- A_late=delta'(0).
- B_late= -delta''(0)/2 - A_late/tauA.
- If H(z)=H_flat(z)[1+delta(z)]/[1-f(z)]^1/2, then q0=q0_flat+A_late+f'(0)/[2(1-f(0))], using f' as derivative w.r.t. z.
- Two independent nonzero-z residual H/H_flat anchors algebraically fix A_late and B_late. But those anchors are observations or background inputs unless the action independently predicts H(z), and fixed tauA,tauB remain ansatz choices.

Finally Omega_x0=1-(omega_b+omega_cdm+omega_r)/(H0/100)^2=**0.7014079815036496** is algebraically closed under flat present-day bookkeeping and is not an independent structural fit coordinate.

## Status ledger

| Coordinate | Historical late300 | No-fit candidate | At present independent physical origin? |
|---|---:|---:|---|
| A_F | 0.020481465 | 0.020520 | No; one-sided health bound ~0.020517950 with derived handoff |
| z_c | 3.927876 | 4.0340775 | Force-balance *conditional on reconstructed trajectory* |
| Delta N_F | 0.33114134 | 0.32925377 | Force profile *conditional on reconstruction* |
| D_floor | 0.051815426 | 0.051815426 | No; rK calculation is inverse mapping |
| D_0 | 0.342319194 | 0.341993212 | Transition correction to D_floor, D_0=D_late-D_floor if endpoints independently derived; magnitude currently data-selected |
| p_D | 1 | 1 | Minimal analytic choice, no unique selector |
| lambda_e | 18.40625 | 18.40625 | No; tracker fraction not independently predicted |
| z_t | 16.193176 | 16.742290 | Conditional helper-dependent fixed-point candidate |
| A_late | 0.024822032 | unchanged | Algebraic if q0 independently predicted |
| B_late | 0.012647043 | unchanged | Algebraic if background curvature / acoustic distance independently predicted |

**Verdict:** All ten coordinates have been investigated in one audit. Their algebraic relationships, conditional constraints, and inverse-problem ambiguity are established. No new first-principles microscopic law fixing all ten is supplied by the presently reconstructed late300 action. This statement is about the currently supplied action/boundary data, not a universal impossibility theorem for every future SDMC completion.

## Observational/model-selection status

Derived-bundle HSC full TATT best available Delta chi^2 = -2.6668948 against local021 (optimizer iteration-limit caveat). For an illustrative independent HSC effective N=60, DeltaAIC=-2.6669+2Delta k and DeltaBIC=-2.6669+ln(60)Delta k.

Converged **conditional nuisance** evidence: Planck Delta ln Z=-2.12549±0.69167, HSC Delta ln Z=+2.0390447±0.6148358, sum -0.0864453±0.9254353, central BF=0.9172. This is not full SDMC-family marginalization over the structural-coordinate prior.

A "zero free parameters after action reconstruction" replay is not a valid reason to declare family Delta k=0: the chosen target history came from a data-informed search. The opposite naive choice Delta k=10 may also overcount if independent theory subsequently removes coordinates. Exact AIC/BIC and family evidence remain conditional until the model family, priors and real independence structure are declared.

## Reproducibility

- Source/test: \`sdmc/audit_late300_all_coordinate_closure.py\` and GitHub Actions \`sdmc_late300_all_coordinate_closure.yml\`.
- Initial all-coordinate pass: run 37801853530, artifact \`sdmc-late300-all-coordinate-closure\`.
- Handoff robustness input: run 37763535074.
- No-Slip action coefficient/spectral replay: run 37800426702.
- AF-D0 and paired transition thresholds: runs 37797826223, 37800023447.
- Independent reduced rK selector audit: run 37800037723.
- Derived-bundle Planck evidence: run 37764900113.
- Derived-bundle HSC evidence: run 37764410555.
- Latest HSC TATT: run 37764359939.

**Next solvable closure step:** Supply/derive independent coefficient functions and initial conditions for F and G2/G3/G4, plus the thermal transfer/activation sector and early structural fraction. Once these equations select the scalar trajectory and H(z) *without using late300 as a reconstruction input*, rerun all ten extraction formulas, cosmological gates, and a properly predeclared model-family Bayesian evidence calculation.


## New follow-up: D0 correction action-level and causal bound (8 October 2026)

### Physical meaning and distinct transition windows

The SDMC author clarified that **D0 was added as a correction to Dfloor**; it is not an unrelated kinetic base. Accordingly,

D(N)=Dfloor + D0 S_D(N)^pD, with S_D(N)=[1+exp(-(N-Nc,D)/wD)]^-1.

Accepted baseline Dfloor=0.05181542627513409, corrective amplitude D0=0.34231919445927034, pD=1, S_D(0)=0.9919691128046833. Hence correction today=0.33957006762377623, D(0)=0.3913854938989103, and the late asymptote Dlate=0.3941346207344044.

**Important hidden shape dependence:** In the "F-only force centre" derived-bundle workflow, the source patch \`sdmc/patch_F_only_force_center.py\` intentionally preserves the **historical D correction window** (z_c,D=3.927876388467848, wD=0.33114133956842123) while moving the **F force window** to z_c,F=4.034077502899133, wF=0.32925377073762596. Treating these two windows as interchangeable gives a wrong D0 causality formula. The derived F force profile therefore does not itself derive the independent kinetic window profile.

### Conditional correction extraction identities

If both endpoints arise from an independently evolved microscopic action, D0 = Dlate-Dfloor. At the D transition midpoint Nc,D=-ln(1+z_c,D), where S_D=1/2 and S_D'=1/(4wD),

D0 = [2^(pD+1) wD/pD] * [dD/dN at Nc,D].

The accepted pD=1 trajectory has midpoint slope=0.2584388851188267; the equation recovers D0=0.34231919445927034. **This is a diagnostic inversion of the chosen profile, not an independent prediction of its derivative.**

### Exact background-preserving kinetic degeneracy

Within the permitted linear-G3 reconstructed action class G2=k1(phi)X+k2(phi)X^2-V(phi), G3=g(phi)X, G4=F(phi)/2, introduce

deltaG2 = mu(phi)*[X-Xbg(phi)]^2.

On the background X=Xbg(phi), deltaG2=deltaG2_X=deltaG2_phi=0, but deltaG2_XX=2mu. Thus the background scalar/Friedmann solution and no-slip F,G3 are preserved while the scalar kinetic response changes. In the existing reconstructed-action conventions, deltaD=4mu Xbg/F. Choosing mu=F*(deltaD0*S_D^pD)/(4Xbg) changes only the correction amplitude along this chosen trajectory. The checked cases are deltaD0=-0.04,0,+0.04, with six tested redshifts each; all action-level identities passed in run **37803894563**. This is a **mathematical non-uniqueness construction**, not an independently motivated or observationally verified alternative theory.

### Subluminal kinetic correction lower bound from real hi_class

The native scalar propagation relation is c_s^2=N_s/[Dfloor+D0*S_D^pD], with N_s the action/background-determined sound-speed numerator. If N_s remains unchanged at fixed background and other EFT functions, gradient positivity only enforces N_s>=0; the condition 0<=c_s^2<=1 gives the *inequality*

D0 >= max_N [(N_s(N)-Dfloor)/S_D(N)^pD] for S_D>0, truncated at zero if the RHS is negative.

A full 23-point hi_class scan at force-derived A_F=0.020520, force F zc=4.0340775 and fixed Dfloor=0.0518154263 found that D0=0.01918 has max c_s^2=1.000116565185 (fails), while D0=0.01920 has max c_s^2=0.9998365574303 (passes); direct numerical interpolation gives a threshold of **D0≈0.019188325854**. The sound-speed numerator inferred by D*c_s^2 was invariant across the sampled z<=100 backgrounds at a maximum relative discrepancy of **8.935e-13**, supporting the fixed-numerator assumption in this specific parameterized calculation. The same pass/fail classification held over the full saved background redshift domain. These checks are in run **37804436088**, with corrected D-window source in follow-up **37804665278**.

For the true accepted D window, the analytic bound from the accepted reference numerator agrees with the numerical transition at **D0≈0.0191883**. An earlier preliminary analytic value ≈0.01917473 was an error caused by incorrectly using the *derived F window* instead of the *historical D correction window*; it must not be cited as the physical causal threshold.

The accepted correction D0=0.3423191945 is much larger than the minimum required by this subluminal criterion, and many smaller corrections are numerically healthy. The inequality cannot independently select the accepted D0.

### Revised complexity conclusion

Dfloor and D0 represent the early baseline and an added, time-dependent corrective response *within one kinetic sector*, but two physical interpretations do not imply two independent microscopic origin mechanisms. If the complete action independently fixes the early and late kinetic support, D0 is a computed difference and should not be separately charged. At present the accepted kinetic amplitude and the legacy D transition window were historically selected; the remaining statistical freedom cannot be removed merely by renaming D0 a correction. Model-family AIC/BIC and Bayesian evidence remain conditional until the kinetic action and boundary conditions predict D(N) without using late300's data-selected profile.

Source scripts: \`sdmc/audit_late300_D0_correction_action_degeneracy.py\`; \`sdmc/run_late300_D0_correction_causality_scan.py\`. Source runs: 37803894563, 37804436088 and 37804665278.
