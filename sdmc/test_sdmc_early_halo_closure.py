#!/usr/bin/env python3
"""Fast reproducibility checks for run_sdmc_early_halo_closure.py."""
import unittest
import math
from run_sdmc_early_halo_closure import age_gyr, masses_and_screening, jacobian_3d, smbhs

class SDMCClosureChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geometry=jacobian_3d(24)
        cls.screen=masses_and_screening()
        cls.bh=smbhs()

    def test_age_matches_late300(self):
        self.assertLess(abs(age_gyr(0)-13.59688902966), 2e-5)
        self.assertGreater(age_gyr(6),age_gyr(20))

    def test_exact_jacobian_invariants(self):
        self.assertLess(self.geometry['max_det_invariant_error'],1e-12)
        self.assertLess(self.geometry['mass_conservation_cellwise_max_error'],1e-12)
        self.assertLess(self.geometry['max_AMGM_violation'],1e-12)

    def test_anisotropic_at_fixed_trace(self):
        case=self.geometry['fixed_trace_case']
        self.assertLess(case['J_aniso'],case['J_iso'])
        self.assertLess(case['first_caustic_D_aniso'],case['first_caustic_D_iso'])

    def test_screening_window_selectivity(self):
        band=self.screen['cutoff_band_kc_Mpc_inverse_90pct_1e11_and_10pct_1e9']
        self.assertLess(band[0],3)
        self.assertGreater(band[1],3)
        xs=self.screen['screening_rows']
        self.assertLess(xs[1]['illustrative_filter_kc_3'],.10)
        self.assertGreater(xs[3]['illustrative_filter_kc_3'],.90)
        self.assertLess(xs[3]['filter_z20_if_fixed_physical_screening'],.90)

    def test_bh_formation_shift(self):
        ages=self.bh['ages_Gyr']
        dt=(ages['20']-ages['23.26'])*1e3
        self.assertAlmostEqual(dt,34.3840745,places=3)
        base=next(x for x in self.bh['growth_rows'] if x['z_observed']==7 and x['z_seed_formation']==20
                  and x['radiative_efficiency']==.1 and x['Eddington_ratio']==1 and x['duty_cycle']==1)
        advanced=next(x for x in self.bh['growth_rows'] if x['z_observed']==7 and x['z_seed_formation']==23.26
                  and x['radiative_efficiency']==.1 and x['Eddington_ratio']==1 and x['duty_cycle']==1)
        ratio=advanced['mass_growth_factor']/base['mass_growth_factor']
        self.assertAlmostEqual(ratio,math.exp((dt/1000)/(base['e_fold_time_Myr']/1000)),places=7)

if __name__=='__main__':unittest.main()