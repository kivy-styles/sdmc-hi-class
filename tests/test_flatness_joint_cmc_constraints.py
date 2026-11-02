"""Exact GR Hamiltonian and momentum initial-data filtering of toy ensemble."""
import unittest
import numpy as np
from sdmc.flatness_probability_ensembles import Ensemble,spatial_positions
from sdmc.flatness_joint_cmc_constraints import (
    CMCSlice, build_gr_cmc_slices,one_prior,positive_R_buchert_bound)

class GRComparisonCMCConstraints(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.e=Ensemble(draws=600)
        cls.rng=np.random.default_rng(100)
        cls.positions=spatial_positions("independent_uniform",8,cls.rng,cls.e,64)
        cls.widths=cls.rng.uniform(.1,.6,64)
        cls.data=build_gr_cmc_slices(cls.positions,cls.widths,cls.e)

    def test_fixed_total_material_inventory(self):
        self.assertAlmostEqual(
           max(abs(self.data["conserved_physical_Mu"]-1)),0.,delta=1e-12)

    def test_gr_hamiltonian_constraint_exactly_satisfied(self):
        self.assertLess(
            max(self.data["Hamiltonian_constraint_max_abs_residual"]),1e-14)

    def test_gr_momentum_constraint_exact_in_isotropic_constant_mean_curvature(self):
        self.assertTrue(self.data["Momentum_constraint_identically_zero"])

    def test_exact_proper_volume_curvature_average_identity(self):
        self.assertLess(max(abs(self.data["geometry_R3_mean"]-
                                self.data["gradient_identity"])),1e-12)

    def test_positive_curvature_mean_for_nontrivial_wells(self):
        self.assertTrue(np.all(self.data["geometry_R3_mean"]>0))

    def test_not_all_template_profiles_are_positive_gr_initial_dust(self):
        self.assertTrue(self.data["GR_CMC_initial_dust_admissible"].any())
        self.assertFalse(self.data["GR_CMC_initial_dust_admissible"].all())

    def test_CMC_Q_exactly_zero(self):
        self.assertTrue(np.all(self.data["Buchert_Q"]==0))
        self.assertTrue(np.all(self.data["Buchert_geometric_rho"]<0))

    def test_none_are_positive_pressureless_Buchert_dust(self):
        self.assertEqual(np.sum(self.data["positive_pressureless_Buchert_dust"]),0)

    def test_if_Rpositive_and_geometric_energy_positive_stress_is_stiff(self):
        for R,Q in [(.01,-.03),(1,-10),(.01,-.011)]:
            r=positive_R_buchert_bound(R,Q)
            self.assertTrue(r["positive_energy"])
            self.assertGreater(r["w"],1.)
            self.assertFalse(r["near_pressureless_positive"])

    def test_Buchert_exact_dust_is_possible_for_negative_R_not_our_family(self):
        d=positive_R_buchert_bound(-.03,-.01)
        self.assertTrue(d["positive_energy"])
        self.assertAlmostEqual(d["w"],0.,delta=1e-14)
        self.assertTrue(d["near_pressureless_positive"])

    def test_counts_are_intersection_on_same_configurations(self):
        r=one_prior("independent_uniform",self.e)
        count=r["counts"]
        self.assertLessEqual(count["flat_wells_gr_cmc"]["count"],
                         count["flat_and_wells"]["count"])
        self.assertLessEqual(count["flat_wells_gr_cmc"]["count"],
                         count["gr_cmc_admissible"]["count"])
        self.assertLessEqual(count["flat_and_gr_cmc"]["count"],
                         count["gr_cmc_admissible"]["count"])
        self.assertLess(r["max_mass_conservation_residual"],1e-11)
        self.assertLess(r["max_hamiltonian_constraint_residual"],1e-13)

    def test_nearly_spaced_constraint_filter_rejects_some(self):
        e=Ensemble(draws=2000)
        r=one_prior("evenly_spaced",e)
        self.assertEqual(r["counts"]["flat"]["count"],2000)
        self.assertLess(r["counts"]["flat_wells_gr_cmc"]["count"],
                        r["counts"]["flat_and_wells"]["count"])
        self.assertLess(r["counts"]["positive_buchert_dust"]["count"],1)
        self.assertGreater(r["counts"]["gr_cmc_admissible"]["count"],1500)

    def test_clustering_template_fails_many_positive_physical_rho(self):
        e=Ensemble(draws=600)
        r=one_prior("clustered",e)
        self.assertEqual(r["counts"]["flat_and_wells"]["count"],0)
        self.assertEqual(r["counts"]["flat_wells_gr_cmc"]["count"],0)

    def test_rigorous_scope_restriction(self):
        r=one_prior("independent_uniform",self.e)
        self.assertTrue(r["GR_CMC_initial_data_not_FULL_SDMC"])
        self.assertTrue(r["observational_Jordan_flatness_not_computed"])
        self.assertTrue(r["physical_CDM_like_gravity_not_computed"])

    def test_bad_inputs_rejected(self):
        with self.assertRaises(ValueError):
            CMCSlice(grid=40)
        with self.assertRaises(ValueError):
            CMCSlice(gravity=0)
        with self.assertRaises(ValueError):
            build_gr_cmc_slices(np.zeros((2,3)),np.ones(2),self.e)
        with self.assertRaises(ValueError):
            positive_R_buchert_bound(.01,-.03,G=-1.)

if __name__=="__main__":
    unittest.main()
