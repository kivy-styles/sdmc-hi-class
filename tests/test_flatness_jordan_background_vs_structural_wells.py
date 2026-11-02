"""Regression: physical Jordan perturbations vs structural homogeneous spatial metric."""
import unittest
from math import isclose
from scipy.special import iv
from sdmc.flatness_jordan_background_vs_structural_wells import (
 ScalarFrameMode,frame_curvature,diagnostic)

class FrameSeparation(unittest.TestCase):
    def test_jordan_exact_homogeneous_and_structural_nonuniform_requires_relative_strain(self):
        s=frame_curvature(ScalarFrameMode(mapper_amplitude=0,
             strain_trace_amplitude=.3,jordan_potential_amplitude=0))
        self.assertTrue(s["jordan_full_exactly_homogeneous"])
        self.assertTrue(s["structural_local_geometry_nontrivial"])
        self.assertFalse(s["lock_L_is_zero"])
        self.assertEqual(s["local_jordan_3R_mean"],0.)
        self.assertGreater(s["local_structural_3R_mean"],0.)

    def test_exact_lock_uniform_mapper_and_jordan_exact_flat_gives_no_well(self):
        s=frame_curvature(ScalarFrameMode(strain_trace_amplitude=0))
        self.assertTrue(s["jordan_full_exactly_homogeneous"])
        self.assertTrue(s["lock_L_is_zero"])
        self.assertTrue(s["mapper_spatially_uniform"])
        self.assertFalse(s["structural_local_geometry_nontrivial"])
        self.assertEqual(s["local_structural_3R_mean"],0.)

    def test_flat_jordan_background_with_local_jordan_potential_remains_perturbed(self):
        s=frame_curvature(ScalarFrameMode(strain_trace_amplitude=0,
                                         jordan_potential_amplitude=.1))
        self.assertEqual(s["background_jordan_K"],0)
        self.assertFalse(s["jordan_full_exactly_homogeneous"])
        self.assertGreater(s["local_jordan_3R_mean"],0)
        self.assertGreater(s["local_structural_3R_mean"],0)

    def test_local_mapper_can_generate_strainless_structural_curvature(self):
        s=frame_curvature(ScalarFrameMode(mapper_amplitude=.1,
                              strain_trace_amplitude=0,
                              jordan_potential_amplitude=0))
        self.assertTrue(s["lock_L_is_zero"])
        self.assertFalse(s["mapper_spatially_uniform"])
        self.assertTrue(s["jordan_full_exactly_homogeneous"])
        self.assertGreater(s["local_structural_3R_mean"],0)
        self.assertFalse(s["physical_Jordan_constraint_solved"])

    def test_oppositely_signed_mapper_strain_can_cancel_local_structural_mode(self):
        s=frame_curvature(ScalarFrameMode(mapper_amplitude=.1,
                     strain_trace_amplitude=-.3,
                     jordan_potential_amplitude=0))
        self.assertFalse(s["structural_local_geometry_nontrivial"])
        self.assertAlmostEqual(s["local_structural_3R_mean"],0,delta=1e-15)

    def test_averaged_metric_3R_is_not_averaged_3R(self):
        s=frame_curvature(ScalarFrameMode())
        self.assertEqual(s["averaged_metric_scalar_curvature"],0.)
        self.assertGreater(s["local_structural_3R_mean"],.009)
        self.assertTrue(s["physical_average_curvature_not_equal_to_metric_average"])

    def test_averaged_metric_normalization_is_bessel(self):
        s=frame_curvature(ScalarFrameMode())
        self.assertAlmostEqual(s["coordinate_metric_mean_factor"],
                               iv(0,.2),places=13)
        self.assertAlmostEqual(s["coordinate_coarse_jordan_scale_ratio"],
                               1.0050000138284305,places=13)
        self.assertAlmostEqual(s["volume_coarse_jordan_scale_ratio"],
                               1.0074861114079738,places=13)

    def test_coarse_grained_rulers_not_unique_at_finite_well_amplitude(self):
        s=frame_curvature(ScalarFrameMode())
        self.assertNotAlmostEqual(s["coordinate_coarse_jordan_scale_ratio"],
                                  s["volume_coarse_jordan_scale_ratio"],places=4)

    def test_weak_deformation_homogeneous_limit_recovers_identity(self):
        s=frame_curvature(ScalarFrameMode(strain_trace_amplitude=3e-5))
        self.assertAlmostEqual(s["coordinate_coarse_jordan_scale_ratio"],
                               1.,delta=1e-10)
        self.assertAlmostEqual(s["volume_coarse_jordan_scale_ratio"],
                               1.,delta=1e-10)

    def test_curvature_scaling_with_structural_ruler(self):
        s=frame_curvature(ScalarFrameMode())
        t=frame_curvature(ScalarFrameMode(Mbar=3))
        self.assertAlmostEqual(t["local_structural_3R_mean"]/
                s["local_structural_3R_mean"],1/9,places=12)

    def test_full_jordan_perturbations_cannot_be_dismissed_by_background_flatness(self):
        s=frame_curvature(ScalarFrameMode(
             jordan_potential_amplitude=.03,
             strain_trace_amplitude=.3))
        self.assertEqual(s["background_jordan_K"],0)
        self.assertGreater(s["local_jordan_linear_R_rms"],0)
        self.assertFalse(s["physical_Jordan_constraint_solved"])

    def test_unsupported_4D_flatness_prediction_disclaimed(self):
        d=diagnostic()
        self.assertIn("not 4D cosmological",d["scope"])
        self.assertIn("full Jordan and structural",str(d["missing"]))

    def test_invalid_scale_rejected(self):
        with self.assertRaises(ValueError):
            ScalarFrameMode(aJ=0)
        with self.assertRaises(ValueError):
            ScalarFrameMode(Mbar=-1)
        with self.assertRaises(ValueError):
            ScalarFrameMode(k=0)
if __name__=="__main__":
    unittest.main()
