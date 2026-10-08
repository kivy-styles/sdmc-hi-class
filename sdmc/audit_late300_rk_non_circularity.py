#!/usr/bin/env python3
"""SDMC late300 r_K selector audit (no observational fitting).

Manuscript B Part I Eq. 8.5: cs_X^2=(1+2r_K)/(1+6r_K).
Continuation VI Sec.125 conditionally identifies D_floor/D_canonical = 1/cs_X^2,
with D_canonical = 16/lambda_e^2.
This audit distinguishes algebraic inversion from an independent microscopic
selection of r_K. The kinetic identity is a reduced early G2 diagnostic, not
a proof that it continues unchanged through late braiding and scalar-tensor
domains.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

LAMBDA_ACCEPTED = 18.40625
D_FLOOR_ACCEPTED = 0.05181542627513409


def kinetic_cs2(r: float) -> float:
    if r < 0:
        raise ValueError("r_K must be nonnegative in the positive G2 limit")
    return (1.0 + 2.0*r)/(1.0 + 6.0*r)


def invert_r_from_ratio(q: float) -> float:
    # q = D_floor / D_can = (1+6*r)/(1+2*r), conditional match.
    if not (1.0 <= q < 3.0):
        raise ValueError("positive-r_K mapping requires 1 <= q < 3")
    return (q-1.0)/(6.0-2.0*q)


def do_audit():
    d_can = 16.0/(LAMBDA_ACCEPTED**2)
    q = D_FLOOR_ACCEPTED/d_can
    r = invert_r_from_ratio(q)
    cs2 = kinetic_cs2(r)
    assert math.isclose(d_can/cs2, D_FLOOR_ACCEPTED, rel_tol=1e-12)
    assert math.isclose((1.0+6.0*r)/(1.0+2.0*r), q, rel_tol=1e-12)
    cs2_bound = 0.9
    r_bound = (1.0-cs2_bound)/(6.0*cs2_bound-2.0)
    d_bound = d_can/cs2_bound
    assert r <= r_bound and d_can <= D_FLOOR_ACCEPTED <= d_bound
    alternatives=[]
    for assumed_cs2 in [1.0,0.975,0.95,0.925,0.90]:
        alternative_d = d_can/assumed_cs2
        inferred_r = invert_r_from_ratio(alternative_d/d_can)
        assert math.isclose(kinetic_cs2(inferred_r), assumed_cs2, rel_tol=1e-12)
        alternatives.append({
            "postulated_cs2_X":assumed_cs2,
            "D_floor_compatible":alternative_d,
            "rK_compatible":inferred_r,
            "stable_reduced_G2":True,
        })
    return {
        "status":"PASS: algebraic map is invertible but rK is not independently selected",
        "assumptions":[
            "positive k1,k2,Z with reduced prethermal G2 = k1 Z + k2 Z^2",
            "conditional identification D_floor/D_can = 1/cs_X^2",
            "accepted lambda_e is input, not derived from an independent early fraction",
            "c_s,X^2 >= 0.9 is a guide/inequality, not a value selector",
        ],
        "accepted":{
            "lambda_e_input":LAMBDA_ACCEPTED,
            "D_floor_input":D_FLOOR_ACCEPTED,
            "D_r_canonical_from_lambda_e":d_can,
            "ratio_Dfloor_over_Dcanonical":q,
            "rK_inverted_from_inputs":r,
            "cs2_X_implied":cs2,
        },
        "causal_guide":{
            "cs2_X_min":cs2_bound,
            "rK_max":r_bound,
            "allowed_Dfloor_interval_at_fixed_lambda":[d_can,d_bound],
        },
        "continuous_counterexamples":alternatives,
        "scientific_conclusion":{
            "confirmed":"Accepted D_floor is consistent with a positive-G2 rK approximately 0.02553 and the causal guide.",
            "not_confirmed":"No independent coefficient ratio k2(Z) Z/k1(Z), early trajectory Z(N), or independent early structural fraction has been derived here.",
            "model_selection":"Do not remove D_floor or lambda_e from full-family complexity solely because their accepted inputs invert to rK. Count them as conditional until action evolution supplies an independent selector.",
            "domain_warning":"The prethermal rK limit, the classical tracker floor, and the late designer low sound-speed sector are not interchangeable.",
        },
    }


if __name__=="__main__":
    result=do_audit()
    path=Path("output/rk_non_circularity_audit.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("SDMC_RK_NONCIRCULARITY",json.dumps(result,sort_keys=True))
