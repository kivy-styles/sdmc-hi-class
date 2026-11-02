"""Conditional near-flatness and redistribution ensemble tests.

An ensemble fraction is defined only relative to its stated prior.
No SDMC dynamical primordial probability is asserted.
"""
import unittest
import numpy as np
from sdmc.flatness_probability_ensembles import (
    Ensemble,conditional_observable_probability,weak_curvature_ratio,
    spatial_positions,wilson_interval,evaluate,fixed_width_mean_ratio,
    diagnostic)

class PrimordialConditionalProbability(unittest.TestCase):
    def test_archived_jordan_mapping_uniform_prior(self):
        self.assertAlmostEqual(
            conditional_observable_probability(1,.004),
            .04357521572455235,places=13)
        self.assertAlmostEqual(
            conditional_observable_probability(1,.001),
            .010893803931138087,places=13)

    def test_narrow_prior_can_be_all_nearflat(self):
        self.assertEqual(conditional_observable_probability(.01,.004),1.)
        self.assertEqual(conditional_observable_probability(.01,.001),1.)

    def test_probability_depends_on_prior_scale_and_mapper(self):
        p1=conditional_observable_probability(1,.004)
        p2=conditional_observable_probability(.1,.004)
        p3=conditional_observable_probability(1,.004,Ck=.1)
        self.assertAlmostEqual(p2,10*p1)
        self.assertAlmostEqual(p3,p2)
        self.assertLess(p1,p2)

    def test_curvature_prior_is_not_cmb_measure(self):
        # This is exactly an assumed distribution of the structural
        # provisional Planck curvature fraction, not Planck/BAO posterior.
        d=diagnostic()
        self.assertFalse(d["user_requested_probability_physically_determined"])
        self.assertIn("not physical SDMC",d["scope"])

    def test_wilson_confidence_is_bounded_with_zero_successes(self):
        for count in (0,10,1000):
            a,b=wilson_interval(count,1000)
            self.assertGreaterEqual(a,0.)
            self.assertLessEqual(b,1.)
            self.assertLessEqual(a,b)
        self.assertGreater(wilson_interval(0,20000)[1],0.)

    def test_equal_center_profiles_do_not_automatically_flatten(self):
        pos=np.zeros((1,8))
        widths=np.array([.1])
        ratio=weak_curvature_ratio(pos,widths)
        self.assertGreater(ratio[0],1.)

    def test_perfectly_even_distribution_cancels_longest_modes(self):
        pos=np.arange(8)[None,:]*(2*np.pi/8)
        ratio=weak_curvature_ratio(pos,np.array([.3]))
        self.assertLess(ratio[0],.01)

    def test_variations_only_positions_and_width_no_new_inventory(self):
        e=Ensemble(draws=1000)
        rng=np.random.default_rng(123)
        for kind in ("independent_uniform","clustered","evenly_spaced"):
            x=spatial_positions(kind,8,rng,e,1000)
            self.assertEqual(x.shape,(1000,8))
            self.assertGreaterEqual(x.min(),0.)
            self.assertLess(x.max(),2*np.pi)
            self.assertEqual(np.sum(np.ones(8)/8),1.)

    def test_weak_ensemble_curvature_ratio_is_nonnegative(self):
        e=Ensemble(draws=500)
        rng=np.random.default_rng(3)
        x=spatial_positions("clustered",8,rng,e,500)
        vals=weak_curvature_ratio(x,np.full(500,.3))
        self.assertTrue(np.isfinite(vals).all())
        self.assertTrue(np.all(vals>=0))

    def test_uniform_random_expected_ratio_exact_formula(self):
        e=Ensemble(draws=2000,min_width=.3,max_width=.3)
        out=evaluate("independent_uniform",e)
        expected=fixed_width_mean_ratio(.3,count=8)
        self.assertAlmostEqual(out["mean_ratio"],expected,
                               delta=.15*expected)

    def test_geometric_near_flatness_is_prior_dependent(self):
        e=Ensemble(draws=3500)
        a=evaluate("independent_uniform",e)
        b=evaluate("clustered",e)
        c=evaluate("evenly_spaced",e)
        pa=a["estimates"]["90pct_reduction"]["probability"]
        pb=b["estimates"]["90pct_reduction"]["probability"]
        pc=c["estimates"]["90pct_reduction"]["probability"]
        self.assertGreater(pa,.27)
        self.assertLess(pa,.42)
        self.assertLess(pb,.01)
        self.assertGreater(pc,.98)

    def test_invalid_prior_inputs_rejected(self):
        with self.assertRaises(ValueError):
            Ensemble(min_width=0)
        with self.assertRaises(ValueError):
            conditional_observable_probability(0,.004)
        with self.assertRaises(ValueError):
            conditional_observable_probability(1,.004,Ck=-1)
        with self.assertRaises(ValueError):
            weak_curvature_ratio(np.zeros((2,8)),np.ones(3))
        with self.assertRaises(ValueError):
            spatial_positions("nonsense",8,np.random.default_rng(),Ensemble(),10)

if __name__=="__main__":
    unittest.main()
