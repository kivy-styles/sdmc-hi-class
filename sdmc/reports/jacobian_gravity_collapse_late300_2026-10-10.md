# SDMC structural/Jacobian reconstruction, gravitational influence, and early-halo collapse

**10 October 2026 — completed continuation; late300 unchanged.**  
**Scientific status:** verified **background identity**, nonunique inverse lapse reconstruction, and controlled **conditional** linear-growth/free-fall proxies. **No independently derived scalar gravitational force or new early-galaxy/SMBH candidate**.

## 1. Sources and exact definitions

Manuscript B §2 and the structural-clock Lagrangian define the physical structural radius `R(a)=R0 sigma(a)`, `sigma=S/S0`; the ordinary/Jacobian-mapped comoving scale `chi_J=R/a=M R0`, with `M=sigma/a`; and `p=dlnS/dln a=1+dlnM/dln a`. The structural-clock differential law `dot R=c N_s` uses lapse `N_s=dt_s/dt_c`; do **not** confuse it with the e-fold independent variable `x=ln a`.

These definitions yield
`p H R = c N_s`, so
`H=c N_s/[p a M R0]`,
and, after differentiation with respect to x,
`q=-1-dlnH/dx = p-1 + dlnp/dx - dlnN_s/dx`.

That is an algebraic identity. **It does not independently determine H, N or an effective gravitational force.** It shows precisely where each of those quantities would enter a genuinely independent inverse-gravity reconstruction if one were separately measured/predicted.

The raw physical/comoving ratio `R/chi_J=a` cannot be a force suppression. The matched coasting-reference comparison `R/(a R0)=M` encodes **relative structural expansion**, not a direct force multiplier.

## 2. Frozen late300 expansion reconstructed

The background is **exactly the frozen effective family** implemented in `sdmc/apply_sdmc_full_background_patch.py` and recovered in `sdmc/run_late300_full_structural_derivation.py`:

`H(z)=H0 sqrt([Omega_m(1+z)^3+Omega_r(1+z)^4+Omega_X0]) (1+delta_H(z))/sqrt(1-f_tr(z))`,

`delta_H(z)=A_late z exp(-z/0.25) - B_late z² exp(-z/1.5)`,

`f_tr=0.5[1-tanh((ln a-ln a_t)/0.5)] * [3(1+w_b)/lambda_e²]`,

with all parameters fixed from `sdmc/late300_candidate.json` and `Omega_r h²=4.17998772e-5`.

The analytic `q(z)=-1-dlnH/dln a` gives:

| z | H (km/s/Mpc) | q | Omega_m(z) | F(z) |
|---:|---:|---:|---:|---:|
| 0 | 69.7148244 | -0.527246735 | 0.29850598 | 1.020524775 |
| 1 | 121.812 | +0.13966343 | 0.78203805 | 1.019405305 |
| 3 | 305.630 | +0.44753757 | 0.99387282 | 1.013453562 |
| 7 | 859.999 | +0.51293016 | 1.00443127 | 1.003857424 |
| 20 | 3688.2 | +0.50655100 | 0.98766774 | 1.000253990 |
| 1090 | 1.581e6 | +0.61985140 | 0.75356512 | 1.0000000017 |

The effective `Omega_m(z)` can be slightly above 1 around z7 because `delta_H<0` and the effective modified sector can be negative under that diagnostic split. This is not a negative *physical matter density*. q is a **background deceleration diagnostic**, not an enhancement of local collapse by itself.

## 3. The crucial inverse-lapse consistency issue

For comparison with the manuscript values, use `N_s0=3.41350846`, `p0=1.03421566`; they infer `R0=(N_s0/p0)c/H0=14193.36794 Mpc`. The earlier recombination mapper is `M(z*=1090)=0.02946748423`, implying the conditional matched target
`Rstar=Mstar a_star R0=0.3833573291 Mpc`, within about 0.4% of the earlier manuscript's rounded `Rstar≈0.3848 Mpc`. **These are not claimed to be an independently co-evolved exact full-action trajectory; they are cross-manuscript endpoint constraints.**

Integrate `dR/dx=c N_s/H` **backward** using the late300 H(z). An erroneously constant lapse `N_s=N_s0` gives

`R_constN(z*) = R0 -(c N_s0/H0)∫_{x*}^0 dx/E(x) = -36.57225 Mpc`.

A negative radius is unacceptable. The interval-weighted mean lapse required to hit the same two geometry endpoints is

`<N_s>_dt = (R0-Rstar)/(c[t0-tstar]) = 3.40464347`

versus `N_s0=3.41350846`; a mean reduction of **~0.260%** is enough to remove an apparent **37 Mpc** accumulated discrepancy. This is a highly cancellation-sensitive integrated-radius consistency condition, not a new measured physical radius.

To **demonstrate nonuniqueness**, two smooth lapse-family inverse reconstructions were calibrated with the same H(z), same `N_s(0)=N_s0`, same `Rstar`, and near-the-same early endpoint:

- Family 1: `N_s(a)=N_s0 [1-0.00362524425 * 4a(1-a)]`.
- Family 2: `N_s(a)=N_s0 [1-0.00451433597 * (27/4) a(1-a)^2]`.

Both pass the frozen q and geometry endpoints, but yield distinct `M(a)` and `p(a)` in between. **They are not derived SDMC lapse solutions.** Their only scientific purpose is to show the need for an additional independent lapse evolution equation.

The finite-difference identity check of `q=p-1+p'/p-N_s'/N_s` at z=0,3,7,20 shows maximum absolute error `1.1031e-8`, numerical rather than theoretical.

## 4. Conditional subhorizon pressureless growth test

With `x=ln a`, the controlled growth equation is

`delta_xx+(1-q)delta_x -(3/2) Omega_m(a) mu(k,a) delta=0`.

Identical approximate growing initial conditions `delta=delta_x=1` were imposed at `z_i=50`. Three **comparators** were integrated:

1. Same late300 H(z), with `mu=1`.
2. Same H(z), with `mu=1/F(z)`, used as a **No-Slip-inspired scalar coupling proxy**, *not* the full scale-dependent Horndeski response `G_eff(k,z)`.
3. Flat LambdaCDM H(z) with same physical densities and `mu=1`.

An arbitrary `mu=1.1` under the same late300 H was run **only to prove nonidentifiability** (different gravity -> different growth with identical Jacobian background); it is not a candidate.

| z | D_late300(mu=1) | D_late300(mu=1/F) | D_LCDM |
|---:|---:|---:|---:|
| 20 | 2.41722424 | 2.41717718 | 2.42302147 |
| 10 | 4.60137261 | 4.60066398 | 4.61724516 |
| 7 | 6.32659216 | 6.32410576 | 6.34212353 |
| 0 | 40.15205236 | 39.65883949 | 39.46579805 |

All D values here are normalized to D(z50)=1; they are not sigma8 directly. At z7, `D(mu=1/F)/D(mu=1)=0.99960699` (only a **0.0393%** suppression). Compared to matched-density LambdaCDM, `D_late300(mu=1/F)/D_LCDM≈0.997159` at z7, i.e. approximately **0.284% less** growth. At z0, the No-Slip proxy gives 1.23% less normalized growth than mu=1 on the *same* late300 H.

**Interpretation:** the already accepted homogeneous mapping and simple `1/F` matter-coupling proxy do **not** deliver the 20%-plus rare-halo amplitude enhancement invoked in the phenomenological FRESCO/JWST treatment. This is not a full action-derived perturbative likelihood: k-dependent braiding, metric potentials, radiation transfer, neutrinos, nonlinear mode coupling and galaxy formation are absent.

## 5. Direct collapse-timescale translation (toy Newtonian bound)

For a uniform pressureless top-hat at actual local density `rho_loc=Delta * rho_bar_m`, a Newtonian free-fall time is

`t_ff=sqrt(3pi/[32 G mu rho_loc])`.

If `rho_bar_m=3 H² Omega_m/(8pi G)`, then

`**t_ff H = pi/[2 sqrt(Delta Omega_m mu)].**`

This is the simplest way to use the reconstructed expansion and **independently supplied** gravitational force to predict a collapse timescale; `Delta` must be physical, not a coordinate-Jacobian volume power.

For **illustrative Delta=201** (overdensity delta=200), `mu=1/F` and the late300 background:

| z | Hubble time (Gyr) | naive pressureless t_ff (Myr) | shift due only to F relative to mu=1 |
|---:|---:|---:|---:|
| 20 | 0.26511 | 29.56 | +0.0127% longer |
| 7 | 1.13703 | 125.94 | +0.1927% longer |

The times are simple gravity-from-rest spherical **free-fall diagnostics**, **not** actual formation/collapse ages: they omit gas pressure, feedback, accretion, non-sphericity, virialisation, radiation and nonlinear scalar terms. The percentage change is `sqrt(F)-1` at fixed rho. They emphasize that **a 2.95% early mapper cannot be read as 97% suppression or enhancement of Newtonian G**.

## 6. What actually follows for the early-galaxy/SMBH theory

**Derivation level:**
- The coasting-reference `M=R/(aR0)` and lapse supply a valid identity for background H and q.
- The constant-lapse test exposes an integrated-radius consistency requirement.
- The reference matter growth and free-fall equations make the missing gravity degree of freedom explicit.

**Not derived:**
- A unique `N_s(a)` from the action independent of the fitted H.
- A unique `mu(k,a)` purely from `M`, without perturbing the accepted action.
- The prior nonlinear `g3` or high-k `alphaK_spatial`, a radiation-to-CDM seed transfer, an early JWST halo-mass function, or SMBH seed/accretion microphysics.

**Next theory-level closure:** jointly evolve the independent structural lapse with H from the accepted action, and extract the scale-dependent dynamical `G_eff(k,a)` from its full linear perturbation equations rather than replacing it with `1/F`. Only then build non-linear halo collapse and compare to JWST/SMBH population distributions. Lapse windows and illustrative mu=1.1 were **not** promoted.

**Reproducibility:** [successful GitHub Actions 38053451173](https://github.com/kivy-styles/sdmc-hi-class/actions/runs/38053451173), with `sdmc/audit_jacobian_gravity_collapse.py` and uploaded `sdmc-jacobian-gravity-collapse-diagnostics` JSON/CSV. An initial run 38053404892 failed a test tolerance for a rounded manuscript value; corrected using a 1e-6 tolerance and reran successfully. Adjacent inherited `sdmc_reconstructed_covariant_evolution.yml` entries show zero executable jobs under this investigation branch (not evidence that the present mathematical diagnostic failed).

**Verdict:** proposed Jacobian/coasting reference is a valid **gravity-expansion background diagnostic**, not a standalone measure of scalar attraction or local collapse. No new candidate; late300 preserved.
