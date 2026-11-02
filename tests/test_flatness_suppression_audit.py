"""Verify conditional flatness suppression and ordinary FLRW kinematics."""
import unittest
import math
from sdmc.flatness_suppression_audit import (
 N_RAD,structural_history,conditional_curvature,required_extra_at_present,
 archival_growth_factors,curvature_mapping_check,nonflat_FLRW_consistency,
 diagnostic)

class SDMCFlatnessSuppressionAudit(unittest.TestCase):
    def test_full_transition_not_just_thermal_handoff(self):
        r=structural_history()
        self.assertAlmostEqual(r["phase1"],2.25)
        self.assertAlmostEqual(r["phase2"],.0381431,delta=2e-7)
        self.assertAlmostEqual(r["total"],.0858219,delta=2e-7)
        self.assertAlmostEqual(r["percent_reduction"],91.4178,delta=.001)

    def test_early_small_curvature_can_grow_later(self):
        x=conditional_curvature(initial=1.)
        self.assertAlmostEqual(x["freezeout"],6.35e-17)
        self.assertAlmostEqual(x["recombination"],.000210696)
        self.assertAlmostEqual(x["present"],.0917953)
        self.assertTrue(x["freezeout"]<x["recombination"]<x["present"])

    def test_target_requires_small_initial_or_c_ratio(self):
        r=required_extra_at_present(initial=1.,target=.001)
        self.assertAlmostEqual(r["maximum_C0"],.010893803931138087,
                               delta=1e-13)
        self.assertAlmostEqual(r["suppression_needed_from_C1"],91.7953,
                               delta=1e-12)
        self.assertFalse(r["adequate_with_C1"])
        q=required_extra_at_present(initial=.01,target=.001)
        self.assertTrue(q["adequate_with_C1"])
        self.assertLess(conditional_curvature(.01)["present"],.001)

    def test_constant_C1_conditional_curvature_linearity(self):
        a=conditional_curvature(.1)
        b=conditional_curvature(.2)
        for key in ("present","recombination","freezeout"):
            self.assertAlmostEqual(b[key]/a[key],2.)

    def test_growth_recombination_to_present(self):
        r=archival_growth_factors()
        self.assertAlmostEqual(r["recombination_to_present_growth"],
                               435.67651972510157,places=9)
        self.assertAlmostEqual(r["inferred_Hstar_over_H0"],
                               22772.2964054,delta=.0001)

    def test_exact_time_only_conformal_bridge_has_Cequal1(self):
        for a,M,H,p in [(1.,1.,1.,1.),(.02,1e-6,2e4,2.),
                        (.1,2.,1.3,1.5),(1e-4,.0008,4e7,2.3)]:
            r=curvature_mapping_check(a,M,H,p,kJ=.01)
            self.assertAlmostEqual(r["omegaJ"],
                         r["omegaJ_from_bridge"],delta=1e-11*
                         max(1.,abs(r["omegaJ"])))
            self.assertAlmostEqual(r["algebraic_ratio"],
                         r["bridge_ratio"],delta=1e-11*
                         max(1.,abs(r["bridge_ratio"])))
            self.assertEqual(r["effective_Ck"],1.)

    def test_nonzero_flatness_not_required_by_locked_bridge(self):
        r=curvature_mapping_check(aJ=.3,M=.3,HJ=2.,p=1.6,kJ=1e-3)
        self.assertNotEqual(r["omegaJ"],0.)
        self.assertAlmostEqual(r["omegaJ"],r["omegaJ_from_bridge"])

    def test_two_constant_k_cannot_give_time_variable_C(self):
        x=nonflat_FLRW_consistency(1.,1.,.5,.5)
        self.assertTrue(x["fixed_structural_k"])
        self.assertTrue(x["fixed_jordan_k"])
        self.assertFalse(x["C_changed"])
        self.assertFalse(x["both_FLRW_fixed_and_C_changed"])

    def test_changing_C_requires_change_of_at_least_one_k(self):
        a=nonflat_FLRW_consistency(1.,1.,1.,.01)
        self.assertTrue(a["C_changed"])
        self.assertFalse(a["fixed_jordan_k"])
        self.assertFalse(a["both_FLRW_fixed_and_C_changed"])
        b=nonflat_FLRW_consistency(1.,.1,1.,1.)
        self.assertTrue(b["C_changed"])
        self.assertFalse(b["fixed_structural_k"])

    def test_empty_geometry_or_unsupported_parameters_rejected(self):
        with self.assertRaises(ValueError):
            nonflat_FLRW_consistency(0.,1.,1.,1.)
        with self.assertRaises(ValueError):
            required_extra_at_present(0.)
        with self.assertRaises(ValueError):
            conditional_curvature(initial=-1.)
        with self.assertRaises(ValueError):
            curvature_mapping_check(-1.,1.,1.,1.,1.)

    def test_physical_nonlinear_evolution_not_claimed(self):
        r=diagnostic()
        self.assertIn("NOT a solved nonflat SDMC",r["scope"])
        self.assertTrue(any("not integrated" in x for x in r["limitations"]))

if __name__=="__main__":
    unittest.main()
