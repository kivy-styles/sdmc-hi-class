#!/usr/bin/env python3
"""Unit checks of analytic Late300 transition relations, not a likelihood test."""
import json
import math
import unittest
from pathlib import Path
import numpy as np
from run_unified_transition_audit import audit

PARAM_PATH=Path(__file__).resolve().parent/'late300_candidate.json'

class UnifiedTransitionChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=json.loads(PARAM_PATH.read_text())['parameters']
        cls.p=p
        cls.a=audit(p)
        cls.e=cls.a['exact_late300_transition']
        cls.t=cls.a['two_exponential_candidate']
        cls.n=cls.a['NoSlip_offshell_test']

    def test_exact_kinetic_gravitational_identity(self):
        self.assertLess(self.e['max_abs_identity_error'],1e-12)
        self.assertLess(self.e['max_abs_derivative_identity_error'],1e-12)
        self.assertAlmostEqual(self.e['C_D_over_lnF'],16.713608919761782,places=11)

    def test_two_distinct_transition_epochs(self):
        self.assertAlmostEqual(self.e['e_fold_separation'],1.2496044314376333,places=11)
        self.assertGreater(self.e['z_early_tracker_handoff'],self.e['z_gravitational_center'])
        self.assertAlmostEqual(self.e['U_at_early_handoff'],0.02245285802552488,places=12)

    def test_canonical_tracker_width_mismatch(self):
        self.assertAlmostEqual(self.t['predicted_width_using_exact_matter_tracker_speed'],0.3610760069393739,places=12)
        self.assertGreater(self.t['field_speed_fractional_excess_over_canonical_matter_tracker'],0.09)
        self.assertLess(self.t['field_speed_fractional_excess_over_canonical_matter_tracker'],0.091)

    def test_simple_two_slope_tracker_can_not_explain_half_shutdown(self):
        self.assertGreater(self.t['canonical_matter_scaling_fraction_from_effective_slope_at_zt'],
                           self.t['canonical_matter_scaling_fraction_from_steep'])
        self.assertLess(self.t['Late300_actual_tracker_fraction_at_zt_with_radiation'],
                        .51*self.t['canonical_matter_scaling_fraction_from_steep'])
        self.assertAlmostEqual(self.t['canonical_matter_fraction_relative_change'],0.0427809447534544,places=11)
        self.assertLess(self.t['Gamma_minus_1_at_zt_for_two_exponential'],.02)

    def test_kinetic_composite_required_by_noslip(self):
        for row in self.e['samples_by_redshift']:
            reconstructed=row['Kcov_required']-6*row['alphaM']+6*row['alphaM']**2
            self.assertAlmostEqual(reconstructed,row['D'],places=12)
            self.assertGreater(row['alphaK'],0)
        self.assertAlmostEqual(self.e['Kcov_at_midpoint'],0.314317128,delta=1e-8)
        self.assertLess(self.e['braid_squared_over_midpoint_D_correction'],.01)

    def test_linear_G3_deviation(self):
        d=self.n['same_phi_velocity_deviations']
        self.assertAlmostEqual(d['1.0']['NoSlip_residual_at_center'],0,places=14)
        self.assertAlmostEqual(d['1.1']['NoSlip_residual_over_alphaM'],-.21,places=12)
        self.assertAlmostEqual(d['0.9']['NoSlip_residual_over_alphaM'],.19,places=12)

    def test_offshell_logarithmic_relation(self):
        # Take M²=1, F_phi=.03, Xstar=.1. For exact log G3, X G3_X=-F_phi/2.
        Fphi=.03;Xstar=.1
        for X in (1e-4,.01,.2,1.,20.):
            # Exact logarithm: G3_X=-(1/2) Fphi / X.
            self.assertAlmostEqual(X*(-Fphi/(2*X)),-Fphi/2,places=13)
            # Regularized logarithm G3=-Fphi/2 ln(1+X/Xstar) has a residual.
            xG3X=X*(-Fphi/(2*(X+Xstar)))
            error=xG3X+Fphi/2
            self.assertAlmostEqual(error/(Fphi/2),Xstar/(X+Xstar),places=13)
        self.assertGreater(self.n['regularized_ratio_test']['1.0']['NoSlip_residual_over_alphaM'],.49)

    def test_future_IR_transition_remains_required(self):
        self.assertLess(self.e['D_finite_formula_asymptote'],0.4)
        self.assertNotAlmostEqual(self.e['D_finite_formula_asymptote'],2.,places=2)
        self.assertEqual(self.a['scientific_verdict']['fundamental_delta_k'],'UNDETERMINED')

if __name__=='__main__':unittest.main(verbosity=2)