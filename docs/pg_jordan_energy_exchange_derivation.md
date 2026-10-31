# SDMC Part III — Jordan-time trace energy exchange and pressurelessness

## 1. Physical accounting

The conserved material inventory M_u is not the quadratic energy of the
supporting space V_s. The special nested-sphere M_u/M_inferred ratio does
not determine the stress of the supporting geometry.

All observable times and the reference spatial volume use Jordan t_c and
a_J. The structural clock gives a_s=S/S0 and

    H_s = (d ln S/dt_c) = p_v H_J.

The clock ratio N=dt_s/dt_c is not the ADM lapse N_ADM. Nothing below
identifies local curvature dilution with a global Friedmann source.

## 2. Restricted action and mode definitions

For one scalar trace mode in the quadratic supporting-space sector of
Manuscript B Part II Section 65,

    K=B f_L^2/3 (>0),
    omega_k²=m_L²+c_L² k²/a_J²,
    u=delta ln M - Phi_J,
    q=ell/3+u,
    C_k=4 W eta_K k^4/a_s^4.

The reduced action density is

    L_k = K/2 dot(ell)^2
          - K omega_k²/2 ell² - C_k q² .

The Jordan a_J^3 volume prefactor in the action gives the trace Euler
Lagrange equation

    ddot(ell) + (3 H_J + dot(K)/K) dot(ell)
      + omega_k² ell + (2 C_k/(3 K))q = 0.                 (A)

The prescribed driver u must eventually be solved by the combined
Jordan lapse, shift, metric and structural scalar constraints.
The equations below are exact **within this reduced mode**. They do
not amount to a fully covariant \rho_g, P_g, momentum or lensing solver.

## 3. The on-shell Jordan energy identity

The reduced mode energy and restricted isotropic Jordan pressure are

    E=K/2 [dot(ell)²+omega_k² ell²] + C_k q²,

    P_red = K/2 [dot(ell)²-m_L² ell²
                    -(c_L² k²/(3a_J²))ell²] - C_k q².

This P_red is obtained by varying the Jordan homogeneous scale at
fixed S, fixed reduced scalar fields and action coefficients.
It does not include the full off-shell Jordan metric coupling.

At time-dependent coefficient functions, the exact on-shell derivative is

    dot(E)=-3H_J K dot(ell)²
          +dot(K)/2 [omega_k² ell²-dot(ell)²]
          +(K/2) dot(omega_k²) ell²
          +dot(C_k) q² + 2 C_k q dot(u),

where

    dot(omega_k²) = dot(m_L²)
                  +k²/a_J² [dot(c_L²)-2H_J c_L²],

    dot(C_k) = (4 k^4/a_s^4)[dot(W eta_K)-4 H_s W eta_K].

Crucially, all time derivatives are Jordan derivatives, but the
curvature coefficient still knows the structural expansion H_s.

Now add 3H_J(E+P_red) and simplify *on shell*. The Hubble-friction
and Jordan spatial-gradient terms cancel, giving

    J_red = dot(E)+3H_J(E+P_red)

          = dot(K)/2 [omega_k² ell²-dot(ell)²]
            +K/2 [dot(m_L²)
                   +dot(c_L²) k²/a_J²] ell²
            +dot(C_k)q²+2C_k q dot(u).                 (B)

This is the **reduced sector energy-exchange residual**, not yet the
physical transfer J_g in a fully constrained FLRW system.

For coefficients K,c_L²,m_L²,W eta_K constant and u stationary,

    J_red = -4 H_s C_k q²
          = -4 p_v H_J E_curvature.

Consequently the **reduced curvature term is generally not a
separately conserved pressureless source**, even if its energy
density is positive. The exchange term reflects structural-scale
variation. In the complete SDMC system this energy work would need
to be balanced by the active scalar, Jordan metric, and/or actual
matter sources, as dictated by the full action.

This is not a statement about the actual global rho_g dilution
law, which must be derived from the full averaged physical stress.

## 4. Independent conservation and pressureless conditions

For a Jordan effective cosmological source, the desired conservation
equation is

    dot(rho_g)+3H_J(rho_g+P_g)=J_g.

A genuinely separately conserved pressureless source satisfies
P_g=0 and J_g=0 simultaneously; neither follows merely from the
existence of a long-lived local gravitational well.

Within the **reduced trace sector**, the *instantaneous* restricted
pressure condition reads

    E_kin=E_mass+E_gradient/3+E_curvature.

The second necessary condition for zero reduced exchange is

    J_red=0,

with J_red from (B). Even if P_red=0 at one epoch, J_red can be
nonzero because a_s and/or couplings evolve. The required condition
on dot(u) can be solved algebraically for nonzero C_k q, but is NOT
an autonomous physical equation for u.

If the external driver and action couplings are constant, and the
curvature term is absent (C_k=0), J_red=0 in the reduced equations.
The pressure generally remains nonzero, however.

This distinguishes:
  (1) an action-derived energy and its restricted spatial stress,
  (2) the exchange of structural-scale/metric energy,
  (3) global Jordan FLRW pressurelessness,
  (4) observable cold-like perturbations.

## 5. Numerical verification

The file sdmc/pg_energy_exchange.py implements the exact energy
derivative, restricted pressure and residual under arbitrary running
K,c_L²,m_L²,W eta_K and H_s/H_J. The tests independently:

- verify the equation of motion and energy terms;
- compute a direct central finite difference of the evolving
  energy along the trace equation, with time-varying coefficients;
- check the full energy-exchange residual against its independently
  simplified expression;
- verify J_red=0 for the constant-parameter, zero-curvature case;
- find nonzero exchange from the structural clock if curvature is active;
- exhibit an instantaneous P_red=0 state with J_red != 0;
- verify that running couplings can exchange energy independently;
- solve *formally* for dot(u) needed to cancel exchange and label
  it a condition, not a derived metric driver.

The demonstration uses H_J=.05 and H_s=.075, so p_v=1.5,
with ell=dot(ell)=0 but nonzero prescribed u=.3 at one instant
and W eta_K=.04, k=2. In that configuration,

    J_red/(H_J E) = -4 p_v = -6.

This is a **single-mode normalized reduced result**, not an
observable transfer rate or fitted \rho_g density.

Run:

    python3 -m unittest discover -s tests -p 'test_pg_energy_exchange.py' -v
    python3 -m sdmc.pg_energy_exchange

## 6. Boundary of the present result

The original quadratic structural action does not supply a unique
physical homogeneous \rho_g from an exactly unexcited flat branch.
A coarse-grained inhomogeneous supporting-space state can have
positive quadratic energy and long-scale modulation, as verified in
the previous 3D tests, but the physical background stress tensor,
action-derived driver u, constraint algebra and energy transfer
with X and ordinary matter are still to be calculated.

No true no-CDM CMB/BAO/weak-lensing fit is claimed.
