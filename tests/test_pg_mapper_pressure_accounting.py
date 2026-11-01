"""SDMC supporting-space mapper-fixed stress and Ward identity tests."""
import unittest
import numpy as np
from sdmc.pg_adiabatic_dilution import (
    WellEra,coefficients,mode_averaged_terms,
    evolve_two_phase_ensemble)
from sdmc.pg_mapper_pressure_accounting import (
    two_restricted_stresses,on_shell_mapper_exchange,
    wkb_mapper_stress,combine_modes_adiabatically,demo)

class MapperPressureClosure(unittest.TestCase):
    def test_two_spatial_variations_difference_is_exact_curvature_work(self):
        for Kin,M,G,C in [(0,0,0,1.),(.4,.2,.3,.6),(1.,0,0,0.),
                          (.4,1.7,.06,.001)]:
            r=two_restricted_stresses(Kin,M,G,C)
            self.assertAlmostEqual(r["pressure_difference"],
                                   4*C/3,places=12)
            self.assertAlmostEqual(r["prediction"],4*C/3,places=12)

    def test_exact_mapper_ward_identity_for_arbitrary_trace_states(self):
        for e in [
              WellEra(m2=.2,G0=.3,Lambda0=.8,k_run=.3,w_run=.5),
              WellEra(m2=0,G0=0,Lambda0=1,p_v=1.5),
              WellEra(m2=0,G0=0,Lambda0=1,p_v=2,h=2),
              WellEra(m2=1,G0=0,Lambda0=0)]:
            for x in [0.,.2,.5,1.]:
                for A,Ax in [(.5,.2),(-.3,1.5),(.1,-.3)]:
                    r=on_shell_mapper_exchange(e,x,A,Ax)
                    self.assertLess(abs(r["action_ward_residual"]),
                                   1e-10*max(1,r["E"]))
                    self.assertAlmostEqual(
                       r["J_difference_over_H"],
                       r["J_difference_expected"],
                       delta=1e-12*max(1,r["E"]))

    def test_WKB_gradient_has_radiationlike_pressure(self):
        r=wkb_mapper_stress(WellEra(
              m2=0,G0=1,Lambda0=0),0)
        self.assertAlmostEqual(r["P_over_E"],1/3)
        self.assertAlmostEqual(r["J_over_H_E"],0.)
        self.assertAlmostEqual(r["dilution_log_slope"],-4.)

    def test_WKB_mass_has_conserved_dustlike_reduced_stress(self):
        r=wkb_mapper_stress(WellEra(
              m2=1,G0=0,Lambda0=0),0)
        self.assertAlmostEqual(r["P_over_E"],0.)
        self.assertAlmostEqual(r["J_over_H_E"],0.)
        self.assertAlmostEqual(r["dilution_log_slope"],-3.)

    def test_matter_tracker_curvature_fixed_mapper_is_not_dust(self):
        r=wkb_mapper_stress(WellEra(
              m2=0,G0=0,Lambda0=1,p_v=1.5),0)
        self.assertAlmostEqual(r["P_over_E"],2/3)
        self.assertAlmostEqual(r["J_over_H_E"],-1.)
        self.assertAlmostEqual(r["dilution_log_slope"],-6.)

    def test_radiation_tracker_curvature_fixed_mapper_is_not_dust(self):
        r=wkb_mapper_stress(WellEra(
              m2=0,G0=0,Lambda0=1,p_v=2,h=2),0)
        self.assertAlmostEqual(r["P_over_E"],2/3)
        self.assertAlmostEqual(r["J_over_H_E"],-2.)
        self.assertAlmostEqual(r["dilution_log_slope"],-7.)

    def test_compensating_curvature_coefficient_gives_a_minus_three_but_not_dust(self):
        e=WellEra(m2=0,G0=0,Lambda0=1,p_v=1.5,w_run=6.)
        r=wkb_mapper_stress(e,0)
        self.assertAlmostEqual(r["dilution_log_slope"],-3.)
        self.assertAlmostEqual(r["P_over_E"],2/3)
        self.assertAlmostEqual(r["J_over_H_E"],+2.)
        self.assertFalse(r["physical_global_dust"])

    def test_general_mixed_coefficient_WKB_Ward_identity(self):
        for e in [
          WellEra(m2=.1,G0=.7,Lambda0=.6,p_v=2,k_run=.2,w_run=4.),
          WellEra(m2=.4,G0=.2,Lambda0=.7,p_v=1.5,k_run=-.2,w_run=1.),
          WellEra(m2=1,G0=0,Lambda0=1,p_v=1.5,k_run=.3,w_run=5.)
        ]:
            for x in [0,.3,.8]:
                r=wkb_mapper_stress(e,x)
                self.assertAlmostEqual(r["identity_residual"],0.,places=12)

    def test_dynamic_two_phase_average_approaches_mapper_WKB_pressure(self):
        e=WellEra(m2=0,G0=0,Lambda0=1,p_v=1.5,H0=.003)
        out=evolve_two_phase_ensemble(e,x_end=.9,num=160)
        x=out["x"][-1]
        mean_old_P=out["P_over_E"][-1]
        mean_curv_f=0.5 # adiabatic curvature domination
        mean_new_P=mean_old_P+(4/3)*mean_curv_f
        self.assertAlmostEqual(mean_new_P,
                  wkb_mapper_stress(e,x)["P_over_E"],delta=.035)

    def test_positive_mixture_cannot_be_exact_dust_if_curvature_fraction_nonzero(self):
        mass=WellEra(m2=1,G0=0,Lambda0=0)
        curvature=WellEra(m2=0,G0=0,Lambda0=1,p_v=1.5)
        mixed=combine_modes_adiabatically([mass,curvature],
                                         [.6,.4],x=.2)
        self.assertGreater(mixed["curvature_fraction"],0.)
        self.assertGreater(mixed["w_effective"],0.)
        self.assertTrue(mixed["indeterminate_abundance"])

    def test_free_spectrum_weights_change_geometric_abundance(self):
        mass=WellEra(m2=1,G0=0,Lambda0=0)
        curv=WellEra(m2=0,G0=0,Lambda0=1,p_v=1.5)
        a=combine_modes_adiabatically([mass,curv],[.95,.05],x=0.)
        b=combine_modes_adiabatically([mass,curv],[.05,.95],x=0.)
        self.assertGreater(b["w_effective"],a["w_effective"])
        self.assertNotEqual(b["J_over_H_E"],a["J_over_H_E"])

    def test_status_does_not_claim_physical_global_pg(self):
        r=demo()
        self.assertIn("mapper field equation",str(r["missing"]))
        self.assertIn("curvature_matter",r["wkb_cases"])
        self.assertAlmostEqual(
            r["wkb_cases"]["curvature_frequency_compensated"]
                             ["dilution_log_slope"],-3.)

if __name__=="__main__":
    unittest.main()
