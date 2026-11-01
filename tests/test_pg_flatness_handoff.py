"""Flatness-branch and surviving-structural-mode handoff diagnostics."""
import unittest
import numpy as np
from sdmc.pg_flatness_handoff import (
  Handoff,equation,two_phase,formal_flatness,
  mature_FLRW_curvature,occupied_trace_energy,
  amplitude_for_fraction,post_handoff_consistency,results)

class FlatnessHandoff(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mass=Handoff(k=.05,W_peak=0.)
        cls.pulse=Handoff(k=.30,W_peak=.20)
        cls.r_mass=two_phase(cls.mass)
        cls.r_pulse=two_phase(cls.pulse)
        cls.rad=Handoff(k=.30,W_peak=.20,p_v=2.,H_power=2.)
        cls.r_rad=two_phase(cls.rad)

    def test_B_activation_positive_and_monotonic(self):
        e=self.pulse
        v0=e.values(-1.6);v1=e.values(0.);v2=e.values(1.6)
        self.assertGreater(v0["B"],e.B_floor)
        self.assertLess(v0["B"],v1["B"])
        self.assertLess(v1["B"],v2["B"])
        self.assertLess(v2["B"],1.)
        for x in [-.7,0.,.7]:
            h=1e-5
            num=(e.values(x+h)["B"]-e.values(x-h)["B"])/(2*h)
            self.assertAlmostEqual(num,e.values(x)["Bprime"],
                        delta=1e-9*max(1,abs(num)))
            self.assertGreater(e.values(x)["K"],0.)

    def test_W_selector_only_around_handoff_and_off_after(self):
        e=self.pulse
        self.assertAlmostEqual(e.values(0.)["W"],e.W_peak,places=13)
        self.assertLess(e.values(-1.6)["W"],.005*e.W_peak)
        self.assertLess(e.values(1.6)["W"],.005*e.W_peak)
        self.assertGreater(e.values(0.)["W"],0.)

    def test_flat_branch_unique_energetic_minimum_while_selector_on(self):
        for x in [-.6,0.,.6]:
            g=formal_flatness(self.pulse,x)
            self.assertTrue(g["flat_energy_minimum"])
            self.assertAlmostEqual(g["curvature_costs"]["0"],0.)
            self.assertGreater(g["curvature_costs"]["-1"],0.)
            self.assertAlmostEqual(
                  g["curvature_costs"]["-1"],
                  g["curvature_costs"]["1"])
            self.assertTrue(g["homogeneous_κ_time_evolution_not_derived"])

    def test_flatness_selection_is_not_mature_evolution(self):
        e=self.pulse
        flat=[mature_FLRW_curvature(0,e,x)["Omega_K"]
                                              for x in [-1.,0.,1.,3.]]
        self.assertEqual(flat,[0.,0.,0.,0.])
        open0=mature_FLRW_curvature(-1,e,0.)["Omega_K"]
        open1=mature_FLRW_curvature(-1,e,1.)["Omega_K"]
        self.assertAlmostEqual(open1/open0,np.exp(2*(e.H_power-1.)))
        self.assertGreater(open1/open0,1.)

    def test_flat_curvature_stays_zero_for_any_finite_expansion(self):
        for h in [.1,.7,1.5,2.,3.]:
            e=Handoff(H_power=h)
            for x in [-2.,0.,1.,3.]:
                self.assertEqual(mature_FLRW_curvature(0,e,x)["Omega_K"],0.)

    def test_selector_ruler_scaling_explicitly_a_s_minus_four(self):
        x=.3;e=self.pulse
        f=Handoff(p_v=2.,k=.30,W_peak=.20)
        ratio=f.values(x)["curvature_selector_intensity"]/(
              e.values(x)["curvature_selector_intensity"])
        self.assertAlmostEqual(ratio,np.exp(-4*(2-1.5)*x),places=13)

    def test_mass_mode_adiabatic_density_survives_B_rise(self):
        r=self.r_mass
        self.assertGreater(r["B"][-1]/r["B"][0],8.)
        self.assertGreater(r["E_ratio"][-1],0.)
        self.assertAlmostEqual(
            r["E_ratio"][-1],r["WKB_E_ratio"][-1],
            delta=.065*r["WKB_E_ratio"][-1])
        self.assertAlmostEqual(
            r["action_invariant"][-1],1.,delta=.06)
        self.assertLess(r["E_curv_fraction"][-1],1e-12)

    def test_amplitude_damping_does_not_imply_zero_energy(self):
        r=self.r_mass
        self.assertLess(r["amplitude_envelope"][-1],.12)
        self.assertAlmostEqual(
            r["amplitude_envelope"][-1],
            r["WKB_envelope_final"],delta=.065*r["WKB_envelope_final"])
        self.assertGreater(r["E_ratio"][-1],.003)

    def test_selector_pulse_can_shutdown_while_energy_survives(self):
        r=self.r_pulse;e=self.pulse
        self.assertLess(r["W"][-1]/e.W_peak,.075)
        self.assertGreater(r["B"][-1],.95)
        self.assertGreater(r["E_ratio"][-1],.001)
        self.assertLess(r["E_curv_fraction"][-1],.002)
        self.assertAlmostEqual(
             r["E_ratio"][-1],r["WKB_E_ratio"][-1],
             delta=.15*r["WKB_E_ratio"][-1])

    def test_radiation_handoff_retains_mode_energy_during_flat_selection(self):
        r=self.r_rad
        self.assertGreater(r["B"][-1],.95)
        self.assertGreater(r["E_ratio"][-1],0.)
        self.assertLess(r["E_curv_fraction"][-1],.002)
        self.assertAlmostEqual(
            r["E_ratio"][-1],r["WKB_E_ratio"][-1],
            delta=.18*r["WKB_E_ratio"][-1])
        self.assertAlmostEqual(r["action_invariant"][-1],1.,delta=.12)

    def test_same_scalar_only_WKB_no_Friedmann_claim(self):
        e=self.mass
        r=self.r_mass
        for a in [r["E_ratio"],r["WKB_E_ratio"],
                  r["action_invariant"],r["amplitude_envelope"]]:
            self.assertTrue(np.isfinite(a).all())
        self.assertAlmostEqual(r["E_mass_fraction"][-1],.5,delta=.035)
        self.assertLess(r["H_over_frequency"][-1],.01)

    def test_exact_locking_mode_has_zero_excitation_energy(self):
        self.assertEqual(occupied_trace_energy(1.,1e3,0.),0.)
        self.assertGreater(occupied_trace_energy(1.,1e3,.001),0.)

    def test_small_field_can_formal_match_any_fraction_but_not_predict_it(self):
        amp=amplitude_for_fraction(.25,1.,1.,1.,1e3)
        self.assertAlmostEqual(amp,np.sqrt(3)/1e3,places=14)
        energy=occupied_trace_energy(1.,1e3,amp)
        self.assertAlmostEqual(energy,3*.25,places=12)
        self.assertLess(amp,.01)
        self.assertNotEqual(amp,
               amplitude_for_fraction(.05,1.,1.,1.,1e3))

    def test_large_kinetic_coefficient_can_suppress_amplitude_at_fixed_energy(self):
        a=amplitude_for_fraction(.2,1.,1.,1.,1e3)
        b=amplitude_for_fraction(.2,1.,1.,4.,1e3)
        self.assertAlmostEqual(b,a/2)

    def test_occupancy_inversion_explicitly_not_prediction(self):
        v=post_handoff_consistency(self.pulse,x_end=1.6)
        self.assertTrue(v["flat_selector_off"])
        self.assertGreater(v["B_final"],.99)
        self.assertLess(v["required_peak_ell_for_assumed_fraction"],.01)
        self.assertFalse(v["empirical_abundance_fitted_or_predicted"])

    def test_W_zero_no_curvature_energy_even_occupied_mode(self):
        e=self.mass
        for x in [-.5,0.,.9]:
            self.assertEqual(e.values(x)["C"],0.)
            self.assertEqual(e.values(x)["Lambda"],0.)
        self.assertGreater(self.r_mass["E_ratio"][-1],0.)

    def test_no_negative_B_or_mass_and_numerical_invalid_rejected(self):
        for args in [{"B_floor":0.},{"width":0.},{"mL":0.},
                     {"W_peak":-1.},{"cL2":-1.}]:
            with self.assertRaises(ValueError):
                Handoff(**args)
        with self.assertRaises(ValueError):
            two_phase(self.mass,x_i=1,x_f=0)
        with self.assertRaises(ValueError):
            amplitude_for_fraction(-1,1,1,1,1)

    def test_full_diagnostic_separates_selection_from_observations(self):
        r=results()
        self.assertIn("mass_only",r["runs"])
        self.assertIn("transient_selector",r["runs"])
        self.assertIn("radiation_handoff",r["runs"])
        self.assertEqual(len(r["homogeneous_curvature_selector"]),4)
        self.assertIn("actual dynamical selection",str(r["missing"]))
        self.assertIn("no-CDM",str(r["missing"]))

if __name__=="__main__":
    unittest.main()
