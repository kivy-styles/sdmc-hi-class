"""Exact post-handoff geometric curvature checks, not FLRW selection."""
import unittest
import numpy as np
from sdmc.flatness_posthandoff_curvature import (
    sinusoid_exact,sinusoid_bessel_mean,periodic_fourier_metrics,example
)
class StructuralPosthandoffCurvature(unittest.TestCase):
    def test_uniform_relative_geometry_has_zero_intrinsic_curvature(self):
        d=sinusoid_exact(0.)
        self.assertEqual(d.physical_mean_exact,0.)
        self.assertEqual(d.linear_rms,0.)

    def test_laplacian_linear_mean_zero_under_periodic_boundaries(self):
        for eps in (.001,.01,.1,-.08):
            d=sinusoid_exact(eps,mode=3)
            self.assertAlmostEqual(d.coordinate_mean_linear,0.,delta=1e-14)

    def test_full_mean_positive_even_if_linear_mean_zero(self):
        for eps in (.001,.02,.1,-.1):
            d=sinusoid_exact(eps,mode=3)
            self.assertGreater(d.physical_mean_exact,0.)
            self.assertAlmostEqual(d.physical_mean_exact,
                 d.physical_mean_integral_identity,delta=1e-11)

    def test_exact_analytic_bessel_value(self):
        for eps in (.001,.01,.03,.1,.2):
            for mode in (1,2,7):
                d=sinusoid_exact(eps,mode=mode)
                theory=sinusoid_bessel_mean(eps,mode=mode)
                self.assertAlmostEqual(d.physical_mean_exact,theory,
                          delta=5e-12*max(1,abs(theory)))

    def test_nonzero_curvature_at_well_center_and_trough(self):
        d=sinusoid_exact(.1)
        self.assertAlmostEqual(d.center_curvature,4*.1*np.exp(-.2),
                               places=12)
        self.assertAlmostEqual(d.trough_curvature,-4*.1*np.exp(.2),
                               places=12)
        self.assertGreater(d.center_curvature,0.)
        self.assertLess(d.trough_curvature,0.)

    def test_weak_limit_mean_is_second_order(self):
        eps=1e-3
        d=sinusoid_exact(eps)
        self.assertAlmostEqual(d.physical_mean_exact/eps**2,1,
                               delta=1e-5)
        self.assertAlmostEqual(d.linear_rms/eps,np.sqrt(8),
                               delta=1e-11)

    def test_curvature_scaling_under_structural_expansion(self):
        d=sinusoid_exact(.1,mode=3,a_s=1.)
        expanded=sinusoid_exact(.1,mode=3,a_s=2.5)
        self.assertAlmostEqual(expanded.physical_mean_exact/
                              d.physical_mean_exact,1/2.5**2,
                              places=12)
        self.assertAlmostEqual(expanded.linear_rms/d.linear_rms,
                               1/2.5**2,places=12)

    def test_local_well_curvature_can_dominate_average(self):
        d=sinusoid_exact(.01)
        self.assertGreater(d.linear_rms,200*d.physical_mean_exact)
        self.assertGreater(d.nonlinear_variance,0)

    def test_three_mode_periodic_field_integral_identity(self):
        d=periodic_fourier_metrics([.03,.01,.007],[1,3,7],
                                  [.2,1.1,2.3])
        self.assertAlmostEqual(d["coordinate_mean_laplacian"],0.,
                               delta=1e-13)
        self.assertAlmostEqual(d["exact_volume_mean_curvature"],
                              d["integration_by_parts_mean"],delta=1e-12)
        self.assertGreater(d["exact_volume_mean_curvature"],0)
        self.assertTrue(d["homogeneous_FLRW_K_not_inferred"])

    def test_weak_spectral_mode_identity(self):
        a=[.003,.001,.0007]
        k=[1,3,7]
        d=periodic_fourier_metrics(a,k,[.2,1.1,2.3])
        want=sum((x*y)**2 for x,y in zip(a,k))
        self.assertAlmostEqual(d["weak_mode_power_mean"],want,places=14)
        self.assertAlmostEqual(d["weak_linear_curvature_variance"],
                    8*sum((x*y*y)**2 for x,y in zip(a,k)),places=14)
        self.assertAlmostEqual(d["exact_volume_mean_curvature"],want,
                  delta=.02*want)

    def test_duplicate_harmonic_does_not_fake_independent_spectrum(self):
        d=periodic_fourier_metrics([.03,.03],[2,2],[0.,0.])
        self.assertIsNone(d["weak_mode_power_mean"])

    def test_zero_well_strength_does_not_prove_cosmological_curvature(self):
        d=periodic_fourier_metrics([0.],[1])
        self.assertEqual(d["exact_volume_mean_curvature"],0.)
        self.assertTrue(d["homogeneous_FLRW_K_not_inferred"])

    def test_validation_rejects_invalid_input(self):
        with self.assertRaises(ValueError):
            sinusoid_exact(.1,a_s=0.)
        with self.assertRaises(ValueError):
            sinusoid_exact(.1,mode=0)
        with self.assertRaises(ValueError):
            periodic_fourier_metrics([.1,.2],[1])
        with self.assertRaises(ValueError):
            periodic_fourier_metrics([.1],[1],a_s=0.)

    def test_demo_qualifies_unknown_primordial_spectrum(self):
        e=example()
        self.assertIn("not dynamical Planck-to-radiation",e["scope"])
        self.assertTrue(any("primordial q_k" in s for s in e["missing"]))

if __name__=="__main__":
    unittest.main()
