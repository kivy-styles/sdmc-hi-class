"""Nonzero background curved S3 geometry tests; no flat initial metric."""
import unittest
from sdmc.flatness_nonflat_structural_slice import (
 closed_s3_curvature,diagnostic)

class NonflatStructuralSlice(unittest.TestCase):
    def test_nonflat_reference_at_zero_deformation(self):
        r=closed_s3_curvature(0.)
        self.assertAlmostEqual(r.physical_average_R3,6.)
        self.assertAlmostEqual(r.curvature_ratio_to_initial,1.)
        self.assertAlmostEqual(r.reduction_percent,0.,delta=1e-12)
        self.assertFalse(r.global_FLRW_curvature_parameter_selected)

    def test_exact_integration_by_parts_on_compact_sphere(self):
        for e in (0.,.1,.5,1.,2.,3.,4.,5.):
            r=closed_s3_curvature(e)
            self.assertAlmostEqual(r.physical_average_R3,
                r.physical_average_gradient_identity,
                delta=2e-12*max(1,r.physical_average_R3))
            self.assertGreater(r.physical_average_R3,0.)

    def test_weak_modulation_of_positive_curvature(self):
        r=closed_s3_curvature(.1)
        self.assertAlmostEqual(r.physical_average_R3,
             5.955267164522356,places=11)
        self.assertAlmostEqual(r.reduction_percent,
             100*(1-5.955267164522356/6),places=10)

    def test_substantial_averaged_suppression_without_Kzero(self):
        r=closed_s3_curvature(3.)
        self.assertAlmostEqual(r.physical_average_R3,
             .18663244545942548,places=11)
        self.assertGreater(r.reduction_percent,96)
        self.assertLess(r.physical_average_R3,1)
        self.assertFalse(r.global_FLRW_curvature_parameter_selected)

    def test_large_modulation_approaches_small_mean_without_exact_zero(self):
        r=closed_s3_curvature(5.)
        self.assertAlmostEqual(r.physical_average_R3,
             .006136525361956885,places=11)
        self.assertLess(r.curvature_ratio_to_initial,.0011)
        self.assertGreater(r.curvature_ratio_to_initial,0.)
        self.assertFalse(r.actual_thermal_density_spectrum_derived)

    def test_local_sign_change_coexists_positive_mean(self):
        r=closed_s3_curvature(2.)
        self.assertLess(r.min_local_R3,0)
        self.assertGreater(r.max_local_R3,0)
        self.assertGreater(r.physical_average_R3,0)

    def test_opposite_epsilon_same_nonflat_average(self):
        for e in (.1,1,3):
            a=closed_s3_curvature(e)
            b=closed_s3_curvature(-e)
            self.assertAlmostEqual(a.physical_average_R3,
               b.physical_average_R3,places=11)

    def test_dimensional_structural_scale_inverse_square(self):
        r1=closed_s3_curvature(2.,a_s=1.)
        r2=closed_s3_curvature(2.,a_s=4.)
        self.assertAlmostEqual(r2.physical_average_R3/
                               r1.physical_average_R3,1/16,places=12)
        self.assertAlmostEqual(r2.curvature_ratio_to_initial,
                               r1.curvature_ratio_to_initial,places=12)

    def test_deformed_volume_changes_while_global_topology_does_not(self):
        r=closed_s3_curvature(3.)
        self.assertGreater(r.physical_volume_ratio_to_unperturbed,1.)
        self.assertFalse(r.global_FLRW_curvature_parameter_selected)

    def test_constraints_not_assumed_at_planck_boundary(self):
        x=diagnostic()
        self.assertIn("prescribed q",x["scope"])
        self.assertTrue(any("prethermal" in y for y in x["missing"]))
        self.assertTrue(all(not r["global_FLRW_curvature_parameter_selected"]
                     for r in x["records"]))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            closed_s3_curvature(1.,a_s=0)
        with self.assertRaises(ValueError):
            closed_s3_curvature(1.,n_quadrature=4)

if __name__=="__main__":
    unittest.main()
