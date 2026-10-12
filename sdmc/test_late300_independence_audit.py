#!/usr/bin/env python3
"""Regression tests for the source-constrained Late300 target independence audit."""
import unittest
import numpy as np
from run_late300_independence_audit import audit,flat_density,P,target,jacobian,THETA

class TestLate300Independence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=audit()

    def test_exact_flat_closure(self):
        om,orr,ox=flat_density()
        self.assertLess(abs(ox-P["Omega_x0"]),1e-12)
        self.assertAlmostEqual(om+orr+ox,1.,places=13)

    def test_all_four_target_sensitivity_ranks(self):
        r=self.result["ranks"]
        self.assertEqual((r["H"]["rank"],r["F"]["rank"],r["D_conditional"]["rank"],r["H_F_D"]["rank"]),
                         (4,3,3,10))
        self.assertGreater(r["H_F_D"]["smallest_over_largest"],0.02)

    def test_mixed_redshift_conditioning(self):
        r=self.result["redshift_conditioning"]
        self.assertEqual(r["1.0"]["rank"],10)
        self.assertEqual(r["30.0"]["rank"],10)
        self.assertGreater(r["30.0"]["smallest_over_largest"],1e-2)
        self.assertLess(r["0.5"]["smallest_over_largest"],1e-5)

    def test_canonical_floor_not_equal_accepted(self):
        d=self.result["D_analytic_checks"]
        self.assertAlmostEqual(d["canonical_radiation_Dfloor"],0.047226890271848634,places=13)
        self.assertGreater(abs(d["fitted_Dfloor"]-d["canonical_radiation_Dfloor"]),0.004)

    def test_2pi_slope_candidate_not_exact(self):
        d=self.result["D_analytic_checks"]
        self.assertAlmostEqual(d["candidate_2pi_over_lambda"],0.34136151074659893,places=13)
        self.assertGreater(abs(d["candidate_2pi_over_lambda"]-d["fitted_D0"]),0.0009)

    def test_target_logarithms_finite(self):
        z=np.expm1(np.linspace(0,np.log1p(1e5),200))
        t=target(z)
        for field in ("lnE","lnF","lnD"):
            self.assertTrue(np.isfinite(t[field]).all())

    def test_finite_difference_derivative(self):
        z=np.array([0.,.2,1.,3.9,16.,40.])
        j=jacobian(z)
        eps=1e-7;p1=dict(P);p2=dict(P)
        p1["A_late"]+=eps;p2["A_late"]-=eps
        numeric=(target(z,p1)["lnE"]-target(z,p2)["lnE"])/(2*eps)
        other=j["lnE"][:,2]/P["A_late"]
        self.assertLess(float(np.max(np.abs(numeric-other))),1e-7)

if __name__=="__main__":unittest.main()
