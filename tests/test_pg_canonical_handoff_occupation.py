"""Post-Planck conditional quantum structural-mode production in finite K handoff."""
import unittest
import numpy as np
from sdmc.pg_canonical_handoff_occupation import (
    CanonicalHandoff,omega_out,K_history,
    canonical_frequency_squared,canonical_mode,
    fast_slow_comparison)

class StructuralCanonicalHandoffOccupation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fast=canonical_mode(CanonicalHandoff(tau=.4))
        cls.medium=canonical_mode(CanonicalHandoff(tau=.8))
        cls.slow=canonical_mode(CanonicalHandoff(tau=1.4))
        cls.unchanged=canonical_mode(CanonicalHandoff(
                            B_initial=.2,B_final=.2,tau=.4))

    def test_kinetic_ramp_remains_positive_but_varies(self):
        e=CanonicalHandoff()
        xs=np.linspace(-5,5,101)
        Ks=np.array([K_history(e,t)["K"] for t in xs])
        self.assertGreater(Ks.min(),0.)
        self.assertAlmostEqual(Ks[0],e.B_initial,delta=.0001)
        self.assertAlmostEqual(Ks[-1],e.B_final,delta=.0001)
        self.assertTrue(np.all(np.diff(Ks)>0))

    def test_positive_finite_effective_frequency_for_selected_fast_ramp(self):
        self.assertGreater(self.fast["min_canonical_omega_squared"],0.)

    def test_quantized_Wronskian_normalization_preserved(self):
        for r in [self.fast,self.medium,self.slow,self.unchanged]:
            self.assertLess(r["max_Wronskian_residual"],2e-8)
            self.assertAlmostEqual(r["alpha_beta_norm"],1.,delta=2e-8)

    def test_fast_kinetic_ramp_creates_nonzero_canonical_mode_occupation(self):
        self.assertGreater(self.fast["n"],.10)
        self.assertLess(self.fast["n"],.20)

    def test_slow_ramp_is_adiabatic_and_has_far_smaller_occupation(self):
        self.assertGreater(self.fast["n"],self.medium["n"])
        self.assertGreater(self.medium["n"],self.slow["n"])
        self.assertLess(self.slow["n"],5e-5)
        self.assertLess(self.slow["n"]/self.fast["n"],1e-3)

    def test_unchanged_kinetic_coefficient_produces_no_quanta(self):
        self.assertLess(self.unchanged["n"],1e-12)
        self.assertAlmostEqual(self.unchanged["min_canonical_omega_squared"],
                               omega_out(CanonicalHandoff())**2,places=12)

    def test_energy_above_vacuum_matches_frequency_times_occupation(self):
        for r in [self.fast,self.medium,self.slow,self.unchanged]:
            self.assertAlmostEqual(r["final_excess_energy"],
                   r["occupation_energy"],delta=2e-8)
            self.assertLess(abs(r["final_excess_energy_error"]),2e-8)

    def test_common_coefficient_rescaling_does_not_change_mode_occupation(self):
        a=canonical_mode(CanonicalHandoff(B_initial=.2,B_final=1.))
        b=canonical_mode(CanonicalHandoff(B_initial=.4,B_final=2.))
        self.assertAlmostEqual(a["n"],b["n"],delta=3e-9)

    def test_larger_mode_frequency_reduces_production_for_same_ramp(self):
        a=canonical_mode(CanonicalHandoff(k=.3,c2=.1))
        b=canonical_mode(CanonicalHandoff(k=3,c2=.1))
        self.assertLess(b["n"],a["n"])

    def test_frequency_formula_from_exact_rescaling_identity(self):
        e=CanonicalHandoff()
        for t in [-.5,-.13,.0,.21,.7]:
            h=K_history(e,t)
            K=h["K"];gamma=h["rate"];gammadot=h["ratedot"]
            q=np.cos(2*t);qd=-2*np.sin(2*t);qdd=-4*q
            A=q/np.sqrt(K)
            Adot=(qd-gamma*q/2)/np.sqrt(K)
            Addot=(qdd-gamma*qd+(gamma**2/4-gammadot/2)*q)/np.sqrt(K)
            lhs=Addot+gamma*Adot+omega_out(e)**2*A
            rhs=(qdd+canonical_frequency_squared(e,t)*q)/np.sqrt(K)
            self.assertAlmostEqual(lhs,rhs,places=12)

    def test_n_is_mode_occupation_not_cdm_density(self):
        d=fast_slow_comparison()
        self.assertFalse(d["physically_fixed_spectrum"])
        self.assertEqual(d["source_curvature_selector_W"],0)
        self.assertTrue(d["classically_exact_zero_field_remains_zero"])
        self.assertIn("renormalized rho_g",str(d["missing"]))

    def test_initial_and_final_ramp_rate_negligible(self):
        e=CanonicalHandoff()
        r=self.fast
        self.assertLess(abs(r["initial_rate"]),5e-5)
        self.assertLess(abs(r["final_rate"]),5e-5)

    def test_equal_B_produces_only_vacuum_fluctuations(self):
        e=CanonicalHandoff(B_initial=.3,B_final=.3)
        d=canonical_mode(e)
        self.assertAlmostEqual(d["final_excess_energy"],0,delta=2e-8)

    def test_canonical_transition_not_valid_at_exact_B_zero(self):
        for kw in [{"B_initial":0},{"B_final":0},{"tau":0},
                   {"m_L":0},{"c2":-1},{"k":-.1},{"t_window":0}]:
            with self.assertRaises(ValueError):
                CanonicalHandoff(**kw)

    def test_limit_handoff_reversed_can_still_mix_vacuum_modes(self):
        r=canonical_mode(CanonicalHandoff(B_initial=1,B_final=.2))
        self.assertGreater(r["n"],0)
        self.assertAlmostEqual(r["alpha_beta_norm"],1.,delta=2e-8)

    def test_input_grid_validation(self):
        with self.assertRaises(ValueError):
            canonical_mode(CanonicalHandoff(),samples=2)

if __name__=="__main__":
    unittest.main()
