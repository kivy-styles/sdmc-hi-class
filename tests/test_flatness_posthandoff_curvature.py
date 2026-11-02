"""Tests: conservatively redistributed M_u and structural curvature at handoff."""
import unittest
import math
import numpy as np
from sdmc.flatness_posthandoff_curvature import (
 exact_conformal_curvature_1d,harmonic_example,
 conserved_gaussian_material,structural_snapshot_from_delta,
 material_fragmentation_case,ensemble_report,diagnostic)

class StructuralCurvatureFromRedistribution(unittest.TestCase):
    def test_constant_structural_q_is_flat_for_zero_background_K(self):
        zeros=np.zeros(1024)
        r=exact_conformal_curvature_1d(zeros,zeros,zeros)
        for key in ("mean","positive_identity","rms",
                    "coordinate_linear_mean"):
            self.assertEqual(r[key],0.)

    def test_exact_physical_volume_mean_equals_positive_gradient_identity(self):
        rng=np.random.default_rng(93)
        x=np.arange(8192)*2*np.pi/8192
        a=.01;b=.015
        q=a*np.cos(3*x)+b*np.sin(7*x)
        dq=-3*a*np.sin(3*x)+7*b*np.cos(7*x)
        lap=-9*a*np.cos(3*x)-49*b*np.sin(7*x)
        r=exact_conformal_curvature_1d(q,dq,lap)
        self.assertGreater(r["mean"],0.)
        self.assertAlmostEqual(r["mean"],r["positive_identity"],
                               delta=2e-14)
        self.assertAlmostEqual(r["coordinate_linear_mean"],0.,
                               delta=5e-15)
        self.assertGreater(r["rms"],r["mean"])

    def test_one_mode_bessel_exact_mean_matches_3curvature_integral(self):
        for k in (1,2,4,8,16):
            r=harmonic_example(k,C=.02)
            self.assertAlmostEqual(r["mean"],r["mean_bessel"],delta=3e-14)
            self.assertAlmostEqual(r["positive_identity"],r["mean"],
                                   delta=3e-14)
            self.assertAlmostEqual(r["source_amplitude"],.02)

    def test_fixed_density_source_fragment_scale_reduces_mean_not_local_RMS(self):
        ss=[harmonic_example(k,C=.02) for k in (1,2,4,8,16)]
        self.assertTrue(all(ss[i]["mean"]>ss[i+1]["mean"]
                       for i in range(len(ss)-1)))
        self.assertAlmostEqual(ss[0]["mean"]/ss[-1]["mean"],255.8,
                               delta=1)
        for s in ss:
            self.assertAlmostEqual(s["rms"],.05656854,delta=2e-5)

    def test_fixed_metric_amplitude_does_not_imply_suppression_by_fragmentation(self):
        # At fixed A (not fixed density source), nonlinear mean ~ A² k²,
        # so adding shorter-wavelength wells can INCREASE mean curvature.
        a=.002
        x=np.linspace(0,2*np.pi,8192,endpoint=False)
        results=[]
        for k in (1,2,4):
            q=a*np.cos(k*x)
            dq=-a*k*np.sin(k*x)
            lap=-a*k*k*np.cos(k*x)
            results.append(exact_conformal_curvature_1d(
                 q,dq,lap)["mean"])
        self.assertLess(results[0],results[1])
        self.assertLess(results[1],results[2])

    def test_nonuniform_gaussian_distributions_have_exact_conserved_Mu(self):
        rng=np.random.default_rng(17)
        positions=rng.uniform(0,2*np.pi,(16,3))
        d=conserved_gaussian_material(positions,.4,N=32)
        self.assertAlmostEqual(d.mean(),1.,delta=4e-15)
        b=conserved_gaussian_material(positions,.4,N=32,
                                      core_fraction=.3)
        self.assertAlmostEqual(b.mean(),1.,delta=4e-15)
        self.assertGreater((b-d).std(),.01)

    def test_three_dimensional_poisson_constraint_and_exact_curvature(self):
        case=material_fragmentation_case(n=8,seed=7,grid=48)
        self.assertAlmostEqual(case["mean_density"],1.,delta=3e-15)
        self.assertLess(case["max_poisson_residual"],2e-14)
        self.assertLess(case["mean_density_contrast"],1e-14)
        self.assertAlmostEqual(case["mean_curvature"],
                   case["integrated_identity"],delta=2e-14)
        self.assertGreater(case["rms_curvature"],
                           case["mean_curvature"]*100)
        self.assertLess(case["max_potential"],.05)

    def test_zero_density_contrast_gives_zero_local_curvature(self):
        r=structural_snapshot_from_delta(np.zeros((32,32,32)))
        self.assertEqual(r["mean_curvature"],0.)
        self.assertEqual(r["rms_curvature"],0.)
        self.assertFalse(r["background_K_R_determined"])
        self.assertFalse(r["jordan_Omega_k_determined"])

    def test_independent_material_densities_give_small_mean_but_finite_wells(self):
        for n in (8,32,128):
            r=material_fragmentation_case(n,seed=7,grid=48)
            self.assertAlmostEqual(r["total_mass_fraction"],1.)
            self.assertEqual(r["redistributed_fraction"],1.)
            self.assertGreater(r["rms_curvature"],.02)
            self.assertLess(r["mean_curvature"],.0001)
            self.assertGreater(r["mean_curvature"],0.)

    def test_ensemble_averaging_suppresses_mean_with_more_smaller_blobs(self):
        r=ensemble_report(ns=(8,32,128),seeds=range(1,7),grid=48)
        xs=r["groups"]
        self.assertGreater(xs[0]["mean_structural_curvature"],
                           xs[1]["mean_structural_curvature"])
        self.assertGreater(xs[1]["mean_structural_curvature"],
                           xs[2]["mean_structural_curvature"])
        self.assertGreater(xs[0]["mean_structural_curvature"]/
                           xs[-1]["mean_structural_curvature"],3.)
        self.assertLess(max(x["mean_local_rms_curvature"] for x in xs)/
                       min(x["mean_local_rms_curvature"] for x in xs),1.2)
        for x in xs:
            self.assertLess(x["max_density_mass_error"],4e-15)
            self.assertLess(x["max_integrated_identity_residual"],2e-13)

    def test_no_primordial_background_K_is_derived_by_flat_reference(self):
        r=material_fragmentation_case(n=8,seed=7,grid=32)
        self.assertFalse(r["background_K_R_determined"])
        self.assertFalse(r["jordan_Omega_k_determined"])
        d=diagnostic()
        self.assertIn("was imposed",d["no_flatness_claim"])
        self.assertIn("NOT primordial",d["gaussian_ensemble"]["scope"])

    def test_invalid_source_rejected(self):
        with self.assertRaises(ValueError):
            harmonic_example(0)
        with self.assertRaises(ValueError):
            harmonic_example(8,resolution=16)
        with self.assertRaises(ValueError):
            conserved_gaussian_material(np.zeros((4,2)),.1,N=32)
        with self.assertRaises(ValueError):
            conserved_gaussian_material(np.zeros((4,3)),.1,N=32,
                                         core_fraction=1.)
        with self.assertRaises(ValueError):
            structural_snapshot_from_delta(np.zeros((4,4,4)))

if __name__=="__main__":
    unittest.main()
