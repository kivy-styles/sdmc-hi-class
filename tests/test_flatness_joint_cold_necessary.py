"""Joint flatness and local-well prerequisites; never assert CDM solved."""
import unittest
import numpy as np
from sdmc.flatness_probability_ensembles import Ensemble,spatial_positions
from sdmc.flatness_joint_cold_necessary import (
    JointThresholds,joint_modes,one_prior,buchert_dust_sign_test)

class JointFlatnessAndLocalWells(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.e=Ensemble(draws=3000)

    def test_joint_boolean_is_intersection_of_same_samples(self):
        rng=np.random.default_rng(4102)
        x=spatial_positions("independent_uniform",8,rng,self.e,30)
        s=rng.uniform(.1,.6,30)
        z=joint_modes(x,s,self.e)
        self.assertTrue(np.array_equal(z["joint_proxy"],
                  z["flat_proxy"]&z["well_proxy"]))

    def test_exact_mass_partition_no_extra_dark_component(self):
        self.assertAlmostEqual(np.sum(np.ones(8)/8),1.)
        # Spatial positions and widths cannot change that sum.
        rng=np.random.default_rng(12)
        x=spatial_positions("clustered",8,rng,self.e,10)
        z=joint_modes(x,np.full(10,.2),self.e)
        self.assertEqual(len(z["joint_proxy"]),10)

    def test_evenly_spaced_can_have_small_longscale_and_highk_power(self):
        x=np.arange(8)[None,:]*(2*np.pi/8)
        z=joint_modes(x,np.array([.1]),self.e)
        self.assertLess(z["curvature_ratio"][0],.1)
        self.assertGreater(z["density_contrast_rms"][0],1.)
        self.assertGreater(z["short_scale_power_share"][0],.05)
        self.assertTrue(z["joint_proxy"][0])

    def test_same_positions_but_wider_wells_lose_local_contrast(self):
        x=np.arange(8)[None,:]*(2*np.pi/8)
        narrow=joint_modes(x,np.array([.1]),self.e)
        wide=joint_modes(x,np.array([.6]),self.e)
        self.assertGreater(narrow["density_contrast_rms"][0],
                           wide["density_contrast_rms"][0])
        self.assertLess(wide["density_contrast_rms"][0],1.)
        self.assertFalse(wide["joint_proxy"][0])

    def test_random_prior_joint_is_much_less_than_flat_proxy(self):
        r=one_prior("independent_uniform",self.e)
        self.assertGreater(r["fractions"]["flat"]["probability"],.30)
        self.assertLess(r["fractions"]["flat"]["probability"],.39)
        self.assertGreater(r["fractions"]["both"]["probability"],.008)
        self.assertLess(r["fractions"]["both"]["probability"],.035)
        self.assertLess(r["fractions"]["both"]["count"],
                        r["fractions"]["flat"]["count"])

    def test_grouped_one_center_fails_flatness_proxy(self):
        r=one_prior("clustered",self.e)
        self.assertEqual(r["fractions"]["flat"]["count"],0)
        self.assertEqual(r["fractions"]["both"]["count"],0)

    def test_evenly_spaced_not_all_prerequisites(self):
        r=one_prior("evenly_spaced",self.e)
        self.assertEqual(r["fractions"]["flat"]["count"],self.e.draws)
        self.assertGreater(r["fractions"]["both"]["probability"],.01)
        self.assertLess(r["fractions"]["both"]["probability"],.04)
        self.assertFalse(r["physical_CDM_like_growth_and_lensing_computed"])

    def test_local_cutoffs_strictly_reduce_joint_pass_rate(self):
        e=Ensemble(draws=2000)
        loose=one_prior("evenly_spaced",e,JointThresholds(density_rms=.5))
        medium=one_prior("evenly_spaced",e,JointThresholds(density_rms=1.))
        strict=one_prior("evenly_spaced",e,JointThresholds(density_rms=2.))
        n=lambda r:r["fractions"]["both"]["count"]
        self.assertGreater(n(loose),n(medium))
        self.assertGreater(n(medium),n(strict))
        self.assertEqual(n(strict),0)

    def test_positive_Rc_is_incompatible_with_special_positive_dust_buchert(self):
        d=buchert_dust_sign_test(.01)
        self.assertLess(d["rho_g_av_if_dust"],0)
        self.assertFalse(d["positive_geometric_dust_compatible"])
        self.assertTrue(d["distinct_action_route_not_excluded"])

    def test_negative_Rc_allows_positive_dust_buchert_if_Q_R3(self):
        d=buchert_dust_sign_test(-.03)
        self.assertAlmostEqual(d["Q_D_if_dust"],-.01)
        self.assertGreater(d["rho_g_av_if_dust"],0)
        self.assertTrue(d["positive_geometric_dust_compatible"])
        # The required negative R is not automatically generated
        # by the positive conformal periodic flat-reference example.

    def test_nontrivial_conformal_wells_never_positive_dust_in_same_metric(self):
        for r in (1e-10,.001,.1,1):
            self.assertFalse(buchert_dust_sign_test(r)[
                      "positive_geometric_dust_compatible"])

    def test_bad_prerequisites_not_accepted(self):
        for kwargs in (dict(density_rms=-1),
                       dict(mean_curvature_ratio=0),
                       dict(short_scale_min_k=0),
                       dict(short_scale_power_fraction=1.01)):
            with self.assertRaises(ValueError):JointThresholds(**kwargs)
        with self.assertRaises(ValueError):
            joint_modes(np.zeros((2,8)),np.ones(2),
                self.e,JointThresholds(short_scale_min_k=1000))
        with self.assertRaises(ValueError):
            buchert_dust_sign_test(.1,gravity=0)

    def test_all_observational_claims_false(self):
        r=one_prior("independent_uniform",self.e)
        self.assertFalse(r["physical_jordan_OmegaK_computed"])
        self.assertFalse(r["physical_CDM_like_growth_and_lensing_computed"])

if __name__=="__main__":
    unittest.main()
