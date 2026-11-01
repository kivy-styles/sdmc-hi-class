"""SDMC occupied mode: conditional scalar-to-strain coupling and exact Ward tests."""
import unittest
import numpy as np
from sdmc.pg_occupied_scalar_curvature_transfer import (
    OccupiedWell,Hubble,omega_L2,omega_sigma2,imposed_Ricci,
    scalar_rhs,l_friction,integrated_linear_response,
    integrate_nonlinear_perturbed_strain,centered_fd_variation,
    occupied_mode_continuity,response_comparison)

class OccupiedStructuralWellScalarLoading(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.e=OccupiedWell()
        cls.o=integrated_linear_response(cls.e,n=301)

    def test_occupied_scalar_response_nonzero(self):
        d=self.o
        self.assertGreater(np.max(abs(d["s"])),1e-4)
        self.assertGreater(np.max(abs(d["delta_A"])),1e-6)
        self.assertGreater(np.max(abs(d["delta_E"])),1e-7)

    def test_exactly_unexcited_mode_has_zero_first_order_response(self):
        e=OccupiedWell(initial_strain=0.,initial_strain_x=0.)
        d=integrated_linear_response(e,n=251)
        self.assertGreater(np.max(abs(d["s"])),1e-4)
        self.assertTrue(np.array_equal(d["A0"],np.zeros_like(d["A0"])))
        self.assertTrue(np.array_equal(d["delta_A"],np.zeros_like(d["delta_A"])))
        self.assertTrue(np.array_equal(d["delta_E"],np.zeros_like(d["delta_E"])))

    def test_beta_zero_removes_sigma_to_strain_forcing(self):
        e=OccupiedWell(beta=0.)
        d=integrated_linear_response(e,n=251)
        self.assertGreater(np.max(abs(d["s"])),1e-4)
        self.assertGreater(np.max(abs(d["A0"])),1e-4)
        self.assertTrue(np.array_equal(d["delta_A"],np.zeros_like(d["delta_A"])))
        self.assertTrue(np.array_equal(d["delta_E"],np.zeros_like(d["delta_E"])))

    def test_Fsigma_zero_removes_curvature_loading(self):
        e=OccupiedWell(F_sigma=0.)
        d=integrated_linear_response(e,n=251)
        self.assertTrue(np.array_equal(d["s"],np.zeros_like(d["s"])))
        self.assertTrue(np.array_equal(d["delta_A"],np.zeros_like(d["delta_A"])))

    def test_external_Ricci_zero_has_no_induced_signal(self):
        e=OccupiedWell(Ricci_pulse=0.)
        d=integrated_linear_response(e,n=251)
        self.assertTrue(np.array_equal(d["s"],np.zeros_like(d["s"])))
        self.assertTrue(np.array_equal(d["delta_A"],np.zeros_like(d["delta_A"])))

    def test_finite_difference_strain_tangent_matches_full_exponential_K(self):
        d=centered_fd_variation(self.e,epsilon=.002,n=301)
        self.assertLess(d["max_abs_strain_linearization_error"],4e-5)
        self.assertGreater(d["max_induced_deltaA"],1e-5)
        self.assertTrue(d["requires_nonzero_background_occupation"])

    def test_finite_difference_energy_tangent_matches_full_exponential_K(self):
        d=centered_fd_variation(self.e,epsilon=.002,n=301)
        self.assertLess(d["max_abs_energy_linearization_error"],4e-5)
        self.assertGreater(d["max_energy_response"],1e-7)

    def test_zero_strain_response_to_any_scalar_forcing_at_full_nonlinear_order(self):
        e=OccupiedWell(initial_strain=0.,initial_strain_x=0.)
        d=integrated_linear_response(e,n=151)
        nonlin=integrate_nonlinear_perturbed_strain(e,d["s"],d["s_x"],
                                             .2,d["x"])
        self.assertTrue(np.array_equal(nonlin["A"],np.zeros_like(nonlin["A"])))
        self.assertTrue(np.array_equal(nonlin["E"],np.zeros_like(nonlin["E"])))

    def test_scalar_amplitude_linear_in_external_curvature_source(self):
        a=integrated_linear_response(OccupiedWell(Ricci_pulse=.12),n=251)
        b=integrated_linear_response(OccupiedWell(Ricci_pulse=.24),n=251)
        self.assertTrue(np.allclose(b["s"],2*a["s"],rtol=3e-6,atol=1e-10))
        self.assertTrue(np.allclose(b["delta_A"],2*a["delta_A"],rtol=3e-6,atol=1e-10))

    def test_induced_structural_response_linear_in_primordial_occupation(self):
        a=integrated_linear_response(OccupiedWell(initial_strain=.12),n=251)
        b=integrated_linear_response(OccupiedWell(initial_strain=.24),n=251)
        self.assertTrue(np.allclose(b["delta_A"],2*a["delta_A"],rtol=1e-6,atol=1e-9))
        self.assertTrue(np.allclose(b["E0"],4*a["E0"],rtol=1e-6,atol=1e-9))

    def test_mode_energy_first_variation_splits_into_K_and_dynamical(self):
        d=self.o
        self.assertTrue(np.allclose(d["delta_E"],
                d["delta_E_direct_K"]+d["delta_E_dynamical"],atol=2e-13))
        self.assertGreater(np.max(abs(d["delta_E_direct_K"])),1e-7)
        self.assertGreater(np.max(abs(d["delta_E_dynamical"])),1e-7)

    def test_exact_W_zero_adiabatic_continuity_identity(self):
        e=self.e
        for x in [-.4,-.15,.1,.5,.8]:
            for A,Ax in [(.1,.3),(-.2,.06),(.04,-.2)]:
                d=occupied_mode_continuity(e,x,A,Ax)
                self.assertLess(abs(d["balance_residual"]),5e-13)
                self.assertGreater(d["E"],0.)

    def test_extra_scalar_perturbation_changes_Kdot_and_exchange(self):
        e=self.e
        bare=occupied_mode_continuity(e,.1,.11,.03)
        excited=occupied_mode_continuity(e,.1,.11,.03,s=.04,sx=.06)
        self.assertNotEqual(bare["K_log_derivative"],
                            excited["K_log_derivative"])
        self.assertNotEqual(bare["J_over_H"],excited["J_over_H"])
        self.assertLess(abs(excited["balance_residual"]),5e-13)

    def test_no_K_run_restores_zero_scalar_exchange(self):
        e=OccupiedWell(beta=0.)
        for x in [-.5,0,.7]:
            d=occupied_mode_continuity(e,x,.1,.13,s=.3,sx=.5)
            self.assertEqual(d["J_over_H"],0.)
            self.assertLess(abs(d["balance_residual"]),2e-13)

    def test_exact_continuity_gradient_scaling_keeps_J_zero(self):
        e=OccupiedWell(beta=0.,c_L2=.7,k_L=2.)
        d=occupied_mode_continuity(e,.2,.1,.24)
        self.assertEqual(d["J_over_H"],0.)
        self.assertLess(abs(d["balance_residual"]),5e-13)

    def test_isotropic_pressure_differs_from_mass_mode_pressure_when_gradient(self):
        e=OccupiedWell(c_L2=1.,k_L=2.)
        d=occupied_mode_continuity(e,.1,.11,.0)
        self.assertLess(d["Pstar"],d["P"])
        self.assertNotEqual(d["P"],d["Pstar"])

    def test_energy_of_empty_mode_remains_zero_in_continuity_identity(self):
        e=self.e
        d=occupied_mode_continuity(e,.1,0.,0.,s=.3,sx=.1)
        self.assertEqual(d["E"],0.)
        self.assertEqual(d["P"],0.)
        self.assertEqual(d["J_over_H"],0.)

    def test_all_scales_in_Jordan_clock(self):
        e=self.e
        self.assertAlmostEqual(Hubble(e,0),e.H0)
        self.assertAlmostEqual(Hubble(e,1),e.H0*np.exp(-e.h))
        self.assertAlmostEqual(omega_L2(e,0),e.m_L2+e.c_L2*e.k_L**2)

    def test_response_summary_does_not_claim_closure(self):
        d=response_comparison()
        self.assertGreater(d["occupied"]["max_perturbed_trace"],0)
        self.assertEqual(d["empty"]["max_perturbed_trace"],0)
        self.assertEqual(d["beta_zero"]["max_perturbed_trace"],0)
        self.assertEqual(d["Ricci_zero"]["max_perturbed_trace"],0)
        self.assertEqual(d["Fsigma_zero"]["max_perturbed_trace"],0)
        self.assertLess(abs(d["continuity_check"]["balance_residual"]),5e-13)
        self.assertTrue(d["finite_difference_crosscheck"]["source_is_prescribed_not_physical_constraint"])
        self.assertIn("derive deltaR_J",str(d["physical_constraints_not_yet_solved"]))

    def test_parameter_validation(self):
        for kw in [{"H0":0},{"K0":0},{"m_L2":0},{"pulse_width":0},
                  {"m_sigma2":0},{"x_end":-.6},{"k_L":-.1}]:
            with self.assertRaises(ValueError):
                OccupiedWell(**kw)
        with self.assertRaises(ValueError):
            integrated_linear_response(self.e,n=4)
        with self.assertRaises(ValueError):
            centered_fd_variation(self.e,epsilon=0.)

if __name__=="__main__":
    unittest.main()
