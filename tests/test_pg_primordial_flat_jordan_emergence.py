"""Tests of SDMC primordial-flat to matter-curved Jordan geometric handoff."""
import unittest
import numpy as np
from sdmc.pg_primordial_flat_jordan_emergence import (
 FlatHandoff,coefficient,energy,source_and_hessian,
 gamma_curvature_mode,split_strain_forces,matched_flat_initial_shape,
 static_positive_B,right_limit_B_zero,
 activation_energy_derivative,constraint_rank_check,diagnostic)

class FlatStructuralJordanHandoff(unittest.TestCase):
    def setUp(self):
        self.e=FlatHandoff()

    def test_primordial_flat_structural_metric_can_coexist_with_curved_jordan(self):
        e=self.e
        ell=matched_flat_initial_shape(e)
        g=gamma_curvature_mode(e,ell)
        self.assertAlmostEqual(ell,-3*e.u)
        self.assertTrue(g["exact_linear_flat"])
        self.assertEqual(g["delta_R3"],0.)
        self.assertNotEqual(e.u,0.)

    def test_flat_initial_strain_depends_on_actual_Jordan_potential(self):
        e=self.e
        for u in [-.12,-.02,.05,.11]:
            ell=matched_flat_initial_shape(e,u=u)
            self.assertAlmostEqual(gamma_curvature_mode(e,ell,u=u)["q"],0.,places=14)

    def test_exact_flat_initial_metric_has_zero_selector_force_but_active_restoring(self):
        e=self.e
        ell=matched_flat_initial_shape(e)
        for B in [.01,.1,.5,.9,1.]:
            f=split_strain_forces(e,B,ell)
            self.assertTrue(f["no_force_from_selector_when_q_zero"])
            self.assertEqual(f["curvature_selector_gradient"],0.)
            self.assertNotEqual(f["restoring_gradient"],0.)
            self.assertAlmostEqual(f["total_gradient"],
                                   f["restoring_gradient"])
        dead=split_strain_forces(e,0,ell)
        self.assertEqual(dead["total_gradient"],0.)

    def test_locking_restoring_drive_makes_gamma_inherit_Jordan_curvature(self):
        e=self.e
        flat=matched_flat_initial_shape(e)
        self.assertEqual(gamma_curvature_mode(e,flat)["q"],0.)
        for B in [.01,.2,.5,.9,1.]:
            eq=static_positive_B(e,B)
            self.assertGreater(eq["q"],0.)
            self.assertLessEqual(eq["q"],e.u)
            self.assertLess(abs(eq["ell"]),abs(flat))
        self.assertAlmostEqual(static_positive_B(e,1)["q"],e.u)

    def test_B_zero_switches_off_both_kinetic_and_selector(self):
        e=self.e;c=coefficient(e,0)
        self.assertEqual(c["K"],0.)
        self.assertEqual(c["W"],0.)
        self.assertEqual(c["C"],0.)
        self.assertIsNone(c["Lambda"])
        self.assertFalse(c["kinetic_rank_nonzero"])
        self.assertTrue(constraint_rank_check(e,0)["no_classical_differential_equation_at_B0"])

    def test_at_B_zero_no_geometric_state_selected_by_displayed_action(self):
        e=self.e
        for ell in [-2.,-.24,0,.13,3.]:
            self.assertEqual(energy(e,0,ell),0.)
            f=source_and_hessian(e,0,ell)
            self.assertEqual(f["strain_force"],0.)
            self.assertEqual(f["H_ell_ell"],0.)
            self.assertEqual(f["H_u_ell"],0.)
        with self.assertRaises(ValueError):
            static_positive_B(e,0)

    def test_W_over_B_has_nonzero_right_limit_but_not_at_exact_B_zero(self):
        e=self.e
        l=right_limit_B_zero(e)
        self.assertGreater(l["Lambda_limit"],0.)
        for B in [1e-2,1e-3,1e-5,1e-7]:
            c=coefficient(e,B)
            self.assertAlmostEqual(c["Lambda"],
                  (1-B)*l["Lambda_limit"],delta=1e-12)
        self.assertTrue(l["B0_lagrangian_vanishes"])
        self.assertFalse(l["initial_strain_determined_at_B_zero"])

    def test_right_limit_of_stationary_strain_not_necessarily_flat(self):
        e=self.e
        near=right_limit_B_zero(e)
        ell=near["conditional_stationary_strain_right_limit"]
        self.assertNotAlmostEqual(ell,matched_flat_initial_shape(e),places=4)
        self.assertNotEqual(
             near["conditional_stationary_curvature_right_limit"],0.)
        self.assertNotEqual(gamma_curvature_mode(e,ell)["q"],0.)

    def test_E_and_Hessian_vanish_linearly_near_B_zero(self):
        e=self.e
        ell=.3
        a=1e-6;b=2e-6
        self.assertAlmostEqual(energy(e,b,ell)/energy(e,a,ell),2.,
                              delta=3e-6)
        H1=source_and_hessian(e,a,ell)["H_ell_ell"]
        H2=source_and_hessian(e,b,ell)["H_ell_ell"]
        self.assertAlmostEqual(H2/H1,2.,delta=3e-6)

    def test_finite_W_sources_Jordan_metric_into_strain_during_handoff(self):
        e=self.e
        for B in [.02,.1,.5,.8,.99]:
            h=source_and_hessian(e,B,0.)
            self.assertNotEqual(h["selector_source"],0.)
            self.assertEqual(h["physical_matter_ell_source_direct"],0.)
            s=static_positive_B(e,B)
            self.assertLess(s["ell"],0.)
            self.assertAlmostEqual(s["strain_force"],0.,delta=1e-13)

    def test_B1_handoff_locking_restores_relative_strain_zero(self):
        e=self.e
        a=static_positive_B(e,1.)
        self.assertEqual(a["ell"],0.)
        self.assertEqual(a["field_energy"],0.)
        self.assertAlmostEqual(a["q"],e.u)
        self.assertNotEqual(a["gamma_R3"],0.)
        self.assertEqual(coefficient(e,1.)["W"],0.)

    def test_metric_to_relative_strain_mixed_hessian_zero_at_both_ends(self):
        e=self.e
        self.assertEqual(source_and_hessian(e,0,0)["H_u_ell"],0.)
        self.assertGreater(source_and_hessian(e,.5,0)["H_u_ell"],0.)
        self.assertEqual(source_and_hessian(e,1.,0)["H_u_ell"],0.)
        self.assertAlmostEqual(
            source_and_hessian(e,.3,0)["H_u_ell"],
            source_and_hessian(e,.7,0)["H_u_ell"])

    def test_energy_derivative_in_B_matches_action_finite_difference(self):
        e=self.e;ell=-.12;u=.14
        for B in [.001,.07,.3,.5,.9,.99]:
            h=1e-6
            num=(energy(e,B+h,ell,u=u)-
                 energy(e,B-h,ell,u=u))/(2*h)
            self.assertAlmostEqual(num,activation_energy_derivative(
                            e,B,ell,u=u),delta=1e-9)

    def test_flat_initial_shape_has_positive_activation_energy_derivative(self):
        e=self.e
        ell=matched_flat_initial_shape(e)
        for B in [0,.2,.5,.8,1.]:
            d=activation_energy_derivative(e,B,ell)
            self.assertGreater(d,0.)
            self.assertAlmostEqual(d,e.fL2/3*e.omega0_sq*ell*ell/4,
                                   delta=1e-13)

    def test_kinetic_momentum_zero_at_B0(self):
        e=self.e
        self.assertEqual(constraint_rank_check(e,0.)["canonical_momentum_for_fixed_velocity"],0.)
        self.assertGreater(constraint_rank_check(e,1e-5)["canonical_momentum_for_fixed_velocity"],0.)

    def test_zero_well_u_has_zero_strain_and_zero_curvature(self):
        e=FlatHandoff(u=0.)
        for B in [.01,.2,.5,.9,1.]:
            s=static_positive_B(e,B)
            self.assertEqual(s["ell"],0.)
            self.assertEqual(s["gamma_R3"],0.)

    def test_zero_k_mode_does_not_derive_local_flatness_condition(self):
        with self.assertRaises(ValueError):
            matched_flat_initial_shape(FlatHandoff(k=0.))
        e=FlatHandoff(k=0.)
        self.assertEqual(coefficient(e,.5)["C"],0.)
        self.assertEqual(static_positive_B(e,.5)["ell"],0.)

    def test_input_validation(self):
        for args in [{"eta":-1},{"fL2":0},{"a_s":0},
                     {"omega0_sq":0},{"k":-1}]:
            with self.assertRaises(ValueError):
                FlatHandoff(**args)
        for B in [-.1,1.1]:
            with self.assertRaises(ValueError):
                coefficient(self.e,B)

    def test_diagnostic_is_conditional_not_CDM_result(self):
        d=diagnostic()
        self.assertTrue(d["physical_initial_state_from_Bzero_undetermined"])
        self.assertTrue(d["source_W_zero_both_Bzero_and_Bone"])
        self.assertTrue(d["material_gravitates_through_Jordan_u_not_inserted_direct_L_source"])
        self.assertTrue(d["independent_extra_CDM_density_not_predicted"])
        self.assertEqual(d["formal_Bzero"]["kinetic_hessian_ell"],0.)
        self.assertIsNone(d["activation_rows"][0]["ell_equilibrium"])

if __name__=="__main__":
    unittest.main()
