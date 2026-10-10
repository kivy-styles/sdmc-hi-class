#!/usr/bin/env python3
"""
SDMC late300 curvature-release and F drift *diagnostic*, not an
autonomous scalar-field or cosmological-likelihood run.

Reads the canonical late300 candidate, evaluates the preselected
'sdmc_full' background and preselected F window exactly as defined by
apply_sdmc_full_background_patch.py and the late300 reconstruction.
Computes the conditional stable-curvature-minimum sign obstruction.

Usage: python3 sdmc/audit_activation_handoff_sign.py
Only Python standard-library dependencies are required.
"""
import json
import math
from pathlib import Path

MPC_KM = 3.0856775814913673e19
SECONDS_PER_YEAR = 365.25 * 86400.0
OMEGA_R_PHYSICAL = 4.17998772e-5


def load():
    path = Path(__file__).with_name("late300_candidate.json")
    return json.loads(path.read_text(encoding="utf-8"))["parameters"]


def make_model(p):
    h = p["H0"] / 100.0
    om = (p["omega_b"] + p["omega_cdm"]) / h**2
    ore = OMEGA_R_PHYSICAL / h**2
    ox = 1.0 - om - ore
    lam = p["lambda_e"]
    nt = -math.log1p(p["z_t"])
    width_tracker = 0.5
    tau_a, tau_b = 0.25, 1.5
    nc = -math.log1p(p["z_c"])
    width_f = p["width"]
    af = p["A_F"]

    def z_of_n(n):
        return math.expm1(-n)

    def e2(n):
        z = z_of_n(n)
        a_inv = 1.0 + z
        rho_m = om * a_inv**3
        rho_r = ore * a_inv**4
        x_r = rho_r / (rho_m + rho_r)
        w = 0.5 * (1.0 - math.tanh((n - nt) / width_tracker))
        f_tracker = w * (3.0 + x_r) / (lam * lam)
        delta = (
            p["A_late"] * z * math.exp(-z / tau_a)
            - p["B_late"] * z**2 * math.exp(-z / tau_b)
        )
        return (rho_m + rho_r + ox) * (1.0 + delta)**2 / (1.0 - f_tracker)

    def f_and_alpha(n):
        # Numerically stable logistic in the inspected domain.
        x = (n - nc) / width_f
        u = 1.0 / (1.0 + math.exp(-x)) if x >= 0 else math.exp(x) / (1.0 + math.exp(x))
        return math.exp(af * u), af / width_f * u * (1.0 - u)

    def dlnh_dn(n):
        eps = 1e-5
        return (math.log(e2(n + eps)) - math.log(e2(n - eps))) / (4 * eps)

    def curvature(n):
        # Dimensionless flat-FLRW Ricci scalar Rcal/H0^2.
        return 6.0 * e2(n) * (2.0 + dlnh_dn(n))

    def curvature_n(n):
        eps = 1e-4
        return (curvature(n + eps) - curvature(n - eps)) / (2 * eps)

    def structural_ratio(z):
        n = -math.log1p(z)
        a_inv = 1.0 + z
        rho_m = om * a_inv**3
        rho_r = ore * a_inv**4
        f, _ = f_and_alpha(n)
        rho_x_j = f * e2(n) - rho_m - rho_r
        return rho_x_j / rho_m

    def deceleration(z):
        n = -math.log1p(z)
        return -1.0 - dlnh_dn(n)

    def bisect_z(fn, lo, hi):
        flo = fn(lo)
        fhi = fn(hi)
        if flo * fhi > 0:
            raise RuntimeError('redshift interval does not bracket a crossing')
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            fm = fn(mid)
            if flo * fm <= 0:
                hi = mid
                fhi = fm
            else:
                lo = mid
                flo = fm
        return 0.5 * (lo + hi)

    def rho_x_j_over_3M2H02(z):
        n = -math.log1p(z)
        f, _ = f_and_alpha(n)
        return f * e2(n) - om * (1.0 + z)**3 - ore * (1.0 + z)**4

    records = []
    for label, z in [
        ("last_scattering", 1090.0),
        ("tracker_handoff_input", p["z_t"]),
        ("F_transition_center_input", p["z_c"]),
        ("today", 0.0),
    ]:
        n = -math.log1p(z)
        fv, alpha = f_and_alpha(n)
        rcal = curvature(n)
        rcal_n = curvature_n(n)
        f_n = fv * alpha
        annual_drift = -alpha * math.sqrt(e2(n)) * p["H0"] * SECONDS_PER_YEAR / MPC_KM
        records.append({
            "label": label, "z": z,
            "F": fv, "alpha_M": alpha, "F_N": f_n,
            "Rcal_over_H0_squared": rcal,
            "dRcal_dN_over_H0_squared": rcal_n,
            "dlnRcal_dN": rcal_n / rcal,
            "background_G_drift_per_year_unscreened": annual_drift,
            "curvature_minimum_implied_meff2_over_M2_Fsigma2_H02":
                (rcal_n / (2.0 * f_n)) if f_n != 0 else None,
        })
    return {
        "scope": "selected-late300-background sign audit, not autonomous action solution",
        "candidate_name": "late300",
        "H0_input_km_s_Mpc": p["H0"],
        "Omega_m0": om, "Omega_r0": ore, "Omega_x0": ox,
        "matter_dilution_factor_at_F_center": (1.0 + p["z_c"])**3,
        "action_split_Q_equals_one_z_from_selected_background": bisect_z(lambda z: structural_ratio(z) - 1.0, 0.3, 0.4),
        "acceleration_q_equals_zero_z_from_selected_background": bisect_z(deceleration, 0.7, 0.8),
        "action_split_rhoX_J_zero_crossings_z_from_selected_background": [
        bisect_z(rho_x_j_over_3M2H02, 4.0, 5.0),
        bisect_z(rho_x_j_over_3M2H02, 8.0, 9.0),
        ],
        "Q_J_today": structural_ratio(0.0),
        "results": records,
        "stable_adiabatic_curvature_minimum_contradicted_at_F_center":
            records[2]["F_N"] > 0 and records[2]["dRcal_dN_over_H0_squared"] < 0,
        "qualification": (
            "Requires independent k1,k2,g,V,F and early initial data to test true "
            "covariant activation and handoff. Nonadiabatic evolution or braiding "
            "can evade this simplified equilibrium sign constraint."
        ),
    }


if __name__ == "__main__":
    results = make_model(load())
    at_center = results["results"][2]
    assert at_center["F_N"] > 0
    assert at_center["dRcal_dN_over_H0_squared"] < 0
    assert abs(at_center["Rcal_over_H0_squared"] - 109.9169) < 0.1
    assert abs(at_center["alpha_M"] - 0.0154628) < 1e-5
    assert abs(results["action_split_Q_equals_one_z_from_selected_background"] - 0.34925) < 0.002
    assert abs(results["acceleration_q_equals_zero_z_from_selected_background"] - 0.70723) < 0.003
    print(json.dumps(results, indent=2))
