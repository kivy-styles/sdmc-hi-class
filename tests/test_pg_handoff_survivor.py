"""Transition-induced surviving structural trace tests, without external CDM."""
import unittest
import numpy as np
from sdmc.pg_handoff_survivor import (
    HandoffScenario,background,rhs,energies,evolve,diagnostic)

class HandoffSurvivor(unittest.TestCase):
    def setUp(self):
        self.s=HandoffScenario()

    def test_jordan_radiation_deceleration_and_structural_coasting(self):
        for t in (0.,3.,5.,10.,30.,60.):
            b=background(t,self.s)
            self.assertAlmostEqual(b["aJ"]**2,1+t/self.s.t0,
                                   places=12)
            self.assertAlmostEqual(b["aS"],1+t/self.s.t0,
                                   places=12)
            self.assertEqual(b["jordan_deceleration"],1.)

    def test_handoff_B_and_W_are_regular_and_selector_switches_off(self):
        a=background(0.,self.s)
        b=background(5.,self.s)
        c=background(60.,self.s)
        self.assertGreater(a["B"],0.)
        self.assertLess(a["B"],.01)
        self.assertAlmostEqual(b["B"],.5,places=12)
        self.assertAlmostEqual(b["W"],1.,places=12)
        self.assertGreater(c["B"],.9999)
        self.assertLess(c["W"],1e-10)
        self.assertLess(c["Lambda"],1e-12)

    def test_dynamic_mode_is_excited_from_zero_initial_state(self):
        r=evolve(self.s)
        self.assertAlmostEqual(r["records"][0]["ell"],0.,places=13)
        self.assertAlmostEqual(r["records"][0]["velocity"],0.,places=13)
        self.assertGreater(abs(r["records"][2]["ell"]),.01)
        self.assertGreater(r["records"][10]["free"],0.)
        self.assertTrue(r["free_tail_nonzero"])

    def test_selector_only_with_zero_external_metric_perturbation_produces_no_mode(self):
        s=HandoffScenario(source_on=False)
        r=evolve(s)
        self.assertEqual(max(abs(a["ell"]) for a in r["records"]),0.)

    def test_no_selector_no_mode_excited_from_zero_state(self):
        s=HandoffScenario(selector_on=False)
        r=evolve(s)
        self.assertEqual(max(abs(a["ell"]) for a in r["records"]),0.)

    def test_surviving_tail_has_approximately_dustlike_comoving_energy(self):
        r=evolve(self.s)
        self.assertLess(r["late_comoving_action_spread_fraction"],.09)
        t20=next(x for x in r["records"] if x["t"]==20.)
        t60=r["records"][-1]
        self.assertGreater(t60["free"],0.)
        self.assertLess(t60["free"],t20["free"])
        self.assertLess(t60["selector"],1e-16)
        self.assertAlmostEqual(t60["comoving_free_action"],
                               t20["comoving_free_action"],
                               delta=.1*t20["comoving_free_action"])

    def test_doubling_external_source_quadruples_excitation_energy(self):
        r1=evolve(HandoffScenario(u0=.01))
        r2=evolve(HandoffScenario(u0=.02))
        self.assertAlmostEqual(
             r2["records"][-1]["ell"]/r1["records"][-1]["ell"],
             2.,places=10)
        self.assertAlmostEqual(
             r2["records"][-1]["free"]/r1["records"][-1]["free"],
             4.,places=10)

    def test_finite_late_response_not_instantaneous_forced_stationarity(self):
        r=evolve(self.s)
        t60=r["records"][-1]
        self.assertLess(t60["W"],1e-10)
        self.assertGreater(abs(t60["ell"]),1e-10)
        self.assertGreater(t60["free"],0.)

    def test_action_energy_components_are_positive_on_stable_branch(self):
        r=evolve(self.s)
        for item in r["records"]:
            self.assertGreaterEqual(item["free"],0.)
            self.assertGreaterEqual(item["selector"],0.)
            self.assertAlmostEqual(item["total"],
                            item["free"]+item["selector"])

    def test_invalid_histories_and_time_bounds_rejected(self):
        with self.assertRaises(ValueError):
            HandoffScenario(width=0.)
        with self.assertRaises(ValueError):
            background(-1.,self.s)
        with self.assertRaises(ValueError):
            evolve(self.s,t_end=-1.)
        with self.assertRaises(ValueError):
            evolve(self.s,t_end=2.,sample=[0.,3.])

    def test_diagnostic_keeps_gravity_observation_unproven(self):
        d=diagnostic()
        self.assertEqual(d["no_perturbation_max_abs_ell"],0.)
        self.assertEqual(d["no_selector_max_abs_ell"],0.)
        self.assertAlmostEqual(d["double_source_free_energy_ratio"],4.,
                               places=10)
        self.assertFalse(d["main"]["observational_gravitational_extra_source_verified"])
        self.assertIn("actual Planck-to-radiation kinetic",str(d["missing"]))

if __name__=="__main__":
    unittest.main()
