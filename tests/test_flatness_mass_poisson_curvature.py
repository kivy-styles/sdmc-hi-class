"""Mass-conserving 3D curvature of a redistributed post-handoff geometry."""
import unittest
import numpy as np
from sdmc.flatness_mass_poisson_curvature import (
 MassGeometry,periodic_gaussian_density,solve_structural_curvature,
 equal_spaced_centers,redistribution_comparison,diagnostic)

class MassPoissonCurvature(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=MassGeometry()
        cls.one=solve_structural_curvature(
            periodic_gaussian_density([(np.pi,)*3],[1.],.6,cls.p),cls.p)
        cls.eight=solve_structural_curvature(
            periodic_gaussian_density(equal_spaced_centers(cls.p),
                                      [1/8]*8,.3,cls.p),cls.p)
        cls.narrow=solve_structural_curvature(
            periodic_gaussian_density(equal_spaced_centers(cls.p),
                                      [1/8]*8,.1,cls.p),cls.p)

    def test_mass_preserved_by_one_to_eight(self):
        for x in (self.one,self.eight,self.narrow):
            self.assertAlmostEqual(x["M_u"],1.,places=12)

    def test_average_curvature_reduces_for_moderate_redistribution(self):
        self.assertAlmostEqual(self.one["physical_mean_R3"],
                           .00013018751559455962,delta=5e-10)
        self.assertAlmostEqual(self.eight["physical_mean_R3"],
                           .00003221578677563654,delta=5e-10)
        reduction=1-self.eight["physical_mean_R3"] / self.one["physical_mean_R3"]
        self.assertAlmostEqual(reduction,.7525431941111349,delta=1e-5)

    def test_max_local_curvature_survives_mean_suppression(self):
        self.assertGreater(self.eight["maximum_local_R3"],
                           self.one["maximum_local_R3"])
        self.assertAlmostEqual(self.one["maximum_local_R3"],
                               .6667786749516493,delta=2e-6)
        self.assertAlmostEqual(self.eight["maximum_local_R3"],
                               .7126635508344074,delta=2e-6)

    def test_overly_compact_lumps_can_increase_average_curvature(self):
        # This rejects any universal law that "breaking the primordial well"
        # must reduce spatial curvature at fixed total material inventory.
        self.assertGreater(self.narrow["physical_mean_R3"],
                           self.one["physical_mean_R3"])
        self.assertGreater(self.narrow["maximum_local_R3"],
                           15*self.one["maximum_local_R3"])

    def test_periodic_laplacian_has_zero_linear_coordinate_mean(self):
        for x in (self.one,self.eight,self.narrow):
            self.assertAlmostEqual(x["linear_coordinate_mean_R3"],
                                   0.,delta=2e-15)

    def test_potential_yields_structural_R3_linear_density_constraint(self):
        for x in (self.one,self.eight,self.narrow):
            self.assertLess(x["poisson_linear_max_error"],2e-12)

    def test_full_proper_mean_matches_positive_gradient_identity(self):
        for x in (self.one,self.eight):
            self.assertGreater(x["physical_mean_R3"],0.)
            self.assertAlmostEqual(x["physical_mean_R3"],
                            x["gradient_identity_R3"],delta=1e-11)

    def test_narrow_field_high_k_aliasing_is_exposed(self):
        # E-fold tiny wells approach grid limits. No false exact equality.
        relative=abs(self.narrow["physical_mean_R3"]-
               self.narrow["gradient_identity_R3"]) / self.narrow["physical_mean_R3"]
        self.assertLess(relative,.01)

    def test_nonuniform_curvature_has_substantial_variance(self):
        for x in (self.one,self.eight):
            self.assertGreater(x["nonlinear_R3_variance"],100*
                               x["physical_mean_R3"]**2)
            self.assertLess(x["q_max"],.1)

    def test_eight_centers_are_distinct_and_all_inside_box(self):
        c=equal_spaced_centers(self.p)
        self.assertEqual(len(set(c)),8)
        self.assertTrue(all(0<=coord<self.p.box_length
                            for center in c for coord in center))

    def test_mass_partition_requires_no_extra_species(self):
        c=equal_spaced_centers(self.p)
        masses=[.5/8]*8
        rho=periodic_gaussian_density(c,masses,.3,self.p)
        self.assertAlmostEqual(rho.sum()*self.p.cell_volume,.5)
        self.assertAlmostEqual(.5+rho.sum()*self.p.cell_volume,1.)
        self.assertAlmostEqual(.5+sum(masses),1.)

    def test_linear_scaling_for_weak_mass_source(self):
        a=periodic_gaussian_density([(np.pi,)*3],[1.],.6,self.p)
        original=solve_structural_curvature(a,self.p)
        doubled=solve_structural_curvature(2*a,self.p)
        self.assertAlmostEqual(doubled["RMS_local_linear_R3"]/
                               original["RMS_local_linear_R3"],2.)
        self.assertGreater(doubled["physical_mean_R3"],
                           original["physical_mean_R3"])

    def test_structural_scale_and_jordan_scale_are_distinct(self):
        c=equal_spaced_centers(self.p)
        density=periodic_gaussian_density(c,[1/8]*8,.3,self.p)
        alt=MassGeometry(aS=2.)
        result=solve_structural_curvature(density,alt)
        self.assertAlmostEqual(result["physical_mean_R3"]/
                        self.eight["physical_mean_R3"],.25,delta=5e-12)

    def test_wrong_grid_bad_mass_invalid(self):
        with self.assertRaises(ValueError):
            MassGeometry(grid=3)
        with self.assertRaises(ValueError):
            periodic_gaussian_density([[1,2]],[1],.3,self.p)
        with self.assertRaises(ValueError):
            periodic_gaussian_density([[1,2,3]],[-1],.3,self.p)
        with self.assertRaises(ValueError):
            solve_structural_curvature(np.zeros((8,8,8)),self.p)

    def test_no_observational_flatness_without_posthandoff_initial_conditions(self):
        d=diagnostic()
        self.assertFalse(d["moderate_redistribution"][
            "actual_planck_to_radiation_physics_solved"])
        self.assertFalse(d["moderate_redistribution"]["after"][
            "primordial_K_J_predicted"])
        self.assertGreater(d["moderate_redistribution"][
            "average_curvature_reduction_percent"],70)
        self.assertLess(d["overly_compact_redistribution"][
            "average_curvature_reduction_percent"],0)

if __name__=="__main__":
    unittest.main()
