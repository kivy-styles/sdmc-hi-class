# SDMC structural lapse N: chronology of N~57.89, N~50.13, and N~3.4135; late300 consistency corrigendum

**10 October 2026 — manuscript-source correction and independent late300 chronology audit**  
**Primary source authority:** Manuscript A Part III §17, Part IX §138, and Structural-Clock Lagrangian Report §§4–6; Manuscript B Part I §§8–14 and §28 (corrigendum-by-continuation). No change to late300.

## A. What N means

`N=dt_s/dt_c=t_P dS/dt_c=dot R/c`, and `p=dlnS/dlna`. Thus `HR/c=N/p`. Do not mistake `x=ln a` for structural N. The lapse is a *differential clock rate*, and an age ratio `t_s/t_c` equals it only on a constant-N trajectory with a compatible origin. `p` can evolve even if `N` remains constant.

## B. The two distinct closures in manuscript chronology

**Historical Part III canonical scalar, under the direct-density identification `rho_X=rho_v` (`chi=1`):**

`N=Xi_v p/sqrt(Omega_phi)`, `Xi_v=sqrt(8pi/3)`. For `lambda_e=20`, radiation `p=2,Omega_phi=4/lambda_e^2=0.01 -> N_r=Xi_v lambda_e=57.88810036`; matter `p=3/2,Omega_phi=3/lambda_e^2=0.0075 -> N_m=sqrt(2pi)lambda_e=50.13256549`. Both are exact fixed-point values within that **old** self-consistent model, and the lapse can change *between* the attractors. The old mature canonical scalar `N_infty=Xi_v=2.894405...` is another distinct asymptotic closure.

**Later accepted bare/active source split (Part IX, Structural-Clock Report, Manuscript B):**

`rho_X=chi_act rho_v` and `N=Xi_v p sqrt(chi_act/Omega_X)`. **The old N=57.89 is not compulsory.** Manuscript B adopts a *representative endpoint-matched reduced handoff* `N: 1 (Planck) -> 2/3 (kinetic crossover) -> 1.6155 (illustrative transfer midpoint) -> 3.41350846 (classical capture)`, with `chi_act=0.00347715, Omega_X=0.01, p=2` at classical radiation capture. It is an existence proof of a smooth bridge, not a unique UV solution. The value `N=50.13` does not survive automatically as the matter-era lapse in this revised branch.

**Chronological classical closure**, Part IX §138:
`S/S0=t_c/t0`, `N=t_P S0/t0=constant`, `p=1/(Ht_c)`. In the model used there, `t0=13.6129780 Gyr -> N=3.41350846`, `p0=1.03421566`; p evolves as radiation 2 -> matter 3/2 -> present ~1.034, **N remains constant on that adopted classical trajectory**. Manuscript B §28.4 states plainly: “On the chronological classical branch, N is approximately constant through radiation and matter tracking.” A microphysical action could allow deviations from that closure; they are not predicted by the numerical values 57/50 by themselves.

## C. Explicit correction to the 2026-10-10 failed radius extrapolation

The earlier `sdmc/audit_jacobian_gravity_collapse.py` took `p0=1.03421566`, `N0=3.41350846` **from Part IX's earlier trajectory** and combined them with late300's `H0=69.71482083` and `H(z)` to set `R0=(N0/p0)c/H0=14193.37 Mpc`. This differs from the earlier structural-clock normalization `R0=c N_earlier t0_earlier=14247.16 Mpc` by **53.79 Mpc**. Integrating with an erroneously constant copied `N=3.41350846` then gave an unphysical negative `R_star=-36.57 Mpc`. **This is a cross-branch consistency failure of the inputs, not evidence that N must vary during radiation/matter.** It is not an accepted physical SDMC prediction.

Repair: keep the source structural normalization constant across the comparison and use the later late300 cosmic age. For `t0_old=13.6129780 Gyr`, `N_old=3.41350846` we have `R0=14247.16 Mpc`. Integrating the frozen late300 background gives `t0_late300≈13.5968940 Gyr`. The same chronological scaling then gives

`N_late300=N_old t0_old/t0_late300≈3.417546` (with small numerical differences from rounded ages).

`p0_late300=1/(H0 t0_late300)≈1.0315297`.

At z*=1090, `t_star≈0.0003677291 Gyr` (367729 years), hence `R_star=c N_late300 t_star≈0.3853157 Mpc`, `M_star=(t_star/t0_late300)(1+1090)≈0.02950618`, and `p_star≈1.6816677`. These agree with manuscript B `R_star≈0.3848 Mpc`, `M_star=0.02946748423` within ~0.1–0.2%. The latter manuscript's `p_star≈1.6814565` also closely agrees. Such agreement validates the chronology **within the adopted late300 background and normalization**, not its fundamental origin.

The chronological-lapse scenario is still an *adopted closure*: the actual covariant kinetic equation must be solved jointly before claiming an all-epoch unique N(a). In particular, prethermal N=1 and 2/3 belongs to a *separate* Planck-to-classical matching stage and should not be extrapolated into the classical age integral without explicit initial conditions.

## D. Implications for structural gravity and early galaxies

The identity `H=cN/(p a M R0)` and derivative `q=p-1+(dlnp/dlna)-(dlnN/dlna)` remain correct. Under *constant classical N*, `q=p-1+p'/p`, **not** always `q=p-1` across radiation-matter transitions. Exact individual trackers with p'=0 give q_r=1 and q_m=1/2. At recombination where the transition is ongoing, q~0.61985 from late300 while p~1.68167; its p' contribution makes up the difference.

Neither historical high N=57/50 nor the newer classical N~3.41 is by itself the strength of local gravitational attraction: `mu(k,a)` needs perturbation dynamics. The previous approximate pressureless growth comparison uses fitted late300 H and conditional mu=1/F; it is not a complete original-triad nonlinear collapse prediction and does not prove an early-galaxy explanation.

## E. Scientific verdict

**Correct:** original 57.89 radiation and 50.13 matter *within the early conditional canonical branch*; newer 3.41350846 classical chronological branch with evolving p *and approximately fixed N*; Planck and prethermal values 1 and 2/3 from a candidate kinetic bridge. These are not all simultaneous epoch values of one solved late300 action.

**Retracted as inference:** earlier unphysical negative-radius test as evidence against constant-N classical closure. Its input normalizations were mismatched. Reconstructed late300 classical chrono gives positive recombination radius ~0.3853 Mpc and mapper ~0.02951.

**Not yet established:** a fundamental, unique Planck->radiation->matter->vacuum N(a) solved from full accepted action and its spatial triad; local gravity or galaxy collapse enhancement from N or M alone.

**Reproducibility:** `sdmc/audit_structural_lapse_chronology_correction.py` and associated GitHub Actions workflow; late300 candidate remains unchanged.
