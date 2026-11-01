"""Source-supported primordial well: exact Eq106.4, conservation and tail tests."""
import unittest
import numpy as np
from sdmc.pg_primordial_well_retarded import (
  PrimordialWell,material_fraction,material_ledger,
  coefficients,static_minimum,driven_acceleration,evolve,
  followup_free_tail,amplitude_curvature_source,response_memory_demo)

class PrimordialWellRetarded(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.e=PrimordialWell()
        cls.s=evolve(cls.e)

    def test_Mu_conserved_exactly_as_material_redistributes(self):
        e=self.e
        for x in np.linspace(-.6,.9,40):
            m=material_ledger(e,x)
            self.assertAlmostEqual(m["M_total"],e.M_u,places=14)
            self.assertAlmostEqual(m["M_remaining"]+m["M_inferred"],
                                   e.M_u,places=14)
            self.assertTrue(0<m["f"]<1)
        self.assertGreater(material_ledger(e,-.5)["M_remaining"],
                           material_ledger(e,.7)["M_remaining"])

    def test_mass_ledger_cosine_density_scales_as_a_minus_three_at_fixed_shape(self):
        e=PrimordialWell(transfer_center=1e4)
        a=material_ledger(e,-.2)
        b=material_ledger(e,.2)
        # negligible shift in the material transfer profile itself
        self.assertAlmostEqual(
            b["material_density_cosine"]/a["material_density_cosine"],
            np.exp(-1.2),places=12)

    def test_material_fourier_shape_broadens_and_decreases(self):
        c=material_ledger(self.e,-.5)
        d=material_ledger(self.e,.65)
        self.assertGreater(c["normalized_mode_shape"],
                           d["normalized_mode_shape"])
        self.assertGreater(c["core_shape"],c["distributed_shape"])
        self.assertNotEqual(c["normalized_mode_shape"],
                            c["material_density_cosine"])

    def test_curvature_selector_has_compact_support_and_is_positive(self):
        e=self.e
        for x in [-1.,-.3,.3,.9]:
            self.assertEqual(coefficients(e,x)["W"],0.)
        for x in [-.15,0,.2]:
            self.assertGreater(coefficients(e,x)["W"],0.)
        self.assertAlmostEqual(coefficients(e,0.)["W"],e.W_peak)

    def test_mismatch_u_not_equated_with_material_density(self):
        e=self.e
        c=coefficients(e,0)
        m=material_ledger(e,0)
        self.assertAlmostEqual(c["u"],
             e.mismatch_initial*m["normalized_mode_shape"],places=14)
        self.assertNotAlmostEqual(c["u"],
             m["material_density_cosine"],places=8)

    def test_exact_manuscript_1064_Lambda_relation(self):
        e=self.e
        for x in np.linspace(-.2,.2,11):
            c=coefficients(e,x)
            self.assertAlmostEqual(c["Lambda"],8*c["W"]*e.eta_K*e.k**4/
               (3*(c["K"]*3)*c["a_s"]**4),places=12)
            self.assertAlmostEqual(c["source_rhs"],-3*c["Lambda"]*c["u"])

    def test_stationary_trace_is_actual_energy_minimum(self):
        e=self.e
        for x in [-.18,0,.13]:
            s=static_minimum(e,x)
            c=coefficients(e,x)
            self.assertAlmostEqual(s["minimum_residual"],0.,
                            delta=1e-12*max(1.,s["E_min"]))
            h=1e-4
            def energy(A):
                return c["K"]*c["omega0_sq"]*A*A/4+(
                  c["C"]/2)*(c["u"]+A/3)**2
            d=(energy(c["static_ell"]+h)-
               energy(c["static_ell"]-h))/(2*h)
            self.assertAlmostEqual(d,0.,delta=1e-10)
            self.assertLess(energy(c["static_ell"]),
                            energy(0))

    def test_curvature_mismatch_suppressed_during_selector(self):
        c=coefficients(self.e,0.)
        self.assertGreater(c["Lambda"],0)
        self.assertLess(c["static_ell"],0)
        self.assertLess(c["static_q"],c["u"])
        self.assertAlmostEqual(
            c["static_q"]/c["u"],
            c["omega0_sq"]/c["total_omega_sq"],places=13)

    def test_W_to_zero_statically_eliminates_trace_support(self):
        e=self.e
        for x in [.28,.5,.85]:
            c=coefficients(e,x)
            self.assertEqual(c["C"],0.)
            self.assertEqual(c["Lambda"],0.)
            self.assertEqual(c["static_ell"],0.)
            self.assertEqual(c["source_rhs"],0.)
            self.assertEqual(c["C_geom_effective"],0.)

    def test_W_zero_everywhere_cannot_drive_trace_from_zero(self):
        e=PrimordialWell(W_peak=0.)
        out=evolve(e,n=200)
        self.assertLess(float(np.max(np.abs(out["ell"]))),1e-13)
        self.assertLess(float(np.max(np.abs(out["energy_strain"]))),1e-13)
        self.assertGreater(coefficients(e,.0)["u"],0)

    def test_k_zero_curvature_lock_does_not_drive_homogeneous_trace(self):
        e=PrimordialWell(k=0.)
        out=evolve(e,n=200)
        self.assertEqual(coefficients(e,0.)["C"],0.)
        self.assertLess(float(np.max(np.abs(out["ell"]))),1e-13)

    def test_structural_trace_deformed_during_handoff(self):
        out=self.s
        self.assertGreater(float(np.max(np.abs(out["ell"]))),1e-4)
        self.assertGreater(float(np.max(np.abs(out["energy_curvature"]))),1e-8)

    def test_tail_has_zero_direct_curvature_source_but_nonzero_excitation(self):
        v=followup_free_tail(self.e,self.s)
        self.assertTrue(v["all_W_exactly_zero"])
        self.assertTrue(v["all_source_exactly_zero"])
        self.assertTrue(v["static_well_exactly_zero"])
        self.assertTrue(v["geometric_energy_still_positive"])
        self.assertGreater(v["final_tail_energy"],0.)
        self.assertNotEqual(v["end_tail_ell"],0.)

    def test_real_trace_end_not_permanent_static_well(self):
        e=self.e;out=self.s
        c=coefficients(e,float(out["x"][-1]))
        self.assertEqual(c["static_ell"],0)
        self.assertNotAlmostEqual(float(out["ell"][-1]),0,places=7)
        # Since W=0, the tail is homogeneous evolution rather than
        # a static deformation sourced by primordial material.
        self.assertEqual(c["source_rhs"],0.)

    def test_higher_W_strengthens_static_lock_response(self):
        mild=PrimordialWell(W_peak=.05)
        strong=PrimordialWell(W_peak=.4)
        a=coefficients(mild,0)
        b=coefficients(strong,0)
        self.assertGreater(abs(b["static_ell"]),abs(a["static_ell"]))
        self.assertLess(b["static_q"],a["static_q"])

    def test_initial_free_excitation_is_separate_from_forcing(self):
        excited=evolve(self.e,excited=True,n=200)
        plain=evolve(self.e,excited=False,n=200)
        self.assertGreater(float(np.max(abs(excited["ell"]-plain["ell"]))),.003)
        self.assertTrue(np.allclose(excited["material_total"],
                                    plain["material_total"]))
        self.assertTrue(np.allclose(excited["source_rhs"],
                                    plain["source_rhs"]))

    def test_rhs_matches_action_gradient_at_fixed_x(self):
        e=self.e
        x=.03; c=coefficients(e,x)
        h=1e-5
        def potential(A):
            return (c["K"]*c["omega0_sq"]*A*A/4+
                       c["C"]/2*(c["u"]+A/3)**2)
        A=.13; Ax=.0
        energygradient=(potential(A+h)-potential(A-h))/(2*h)
        expected_acc=-2*energygradient/(c["K"]*c["H"]**2)
        self.assertAlmostEqual(driven_acceleration(e,x,A,Ax),
                 expected_acc,delta=1e-8)

    def test_canonical_reduced_energy_well_not_material_source_counted_twice(self):
        c=coefficients(self.e,0)
        s=static_minimum(self.e,0)
        self.assertAlmostEqual(
             s["E_min"],s["E_bare"]+s["E_curv"],places=14)
        self.assertLess(s["E_min"],c["C"]*c["u"]**2/2)
        self.assertAlmostEqual(material_ledger(self.e,0)["M_total"],
                               self.e.M_u)

    def test_limit_W_off_with_nonzero_metric_mismatch_is_not_a_well(self):
        c=coefficients(self.e,.7)
        self.assertGreater(c["u"],0.)
        self.assertEqual(c["source_rhs"],0.)
        self.assertEqual(static_minimum(self.e,.7)["E_min"],0.)

    def test_traceless_not_in_this_scalar_source_calculation(self):
        res=response_memory_demo()
        self.assertIn("metric constraints",str(res["missing"]))
        self.assertIn("no direct material-to-L source",res["origin"])
        self.assertEqual(res["final"]["W"],0.)
        self.assertGreater(res["max_ell"],0.)
        self.assertGreater(res["final"]["energy"],0.)

    def test_invalid_conditions(self):
        with self.assertRaises(ValueError):
            PrimordialWell(B_floor=0.)
        with self.assertRaises(ValueError):
            PrimordialWell(m2=0.)
        with self.assertRaises(ValueError):
            PrimordialWell(W_peak=-1.)
        with self.assertRaises(ValueError):
            evolve(self.e,x_i=1,x_f=0)

if __name__=="__main__":
    unittest.main()
