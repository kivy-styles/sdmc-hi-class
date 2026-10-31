"""Tests of structural-well virial pressure vs Jordan FLRW dilution.

These concern the existing reduced supporting-space trace (not a new
CDM source). Local curvature dilution is not identified with global
Friedmann density and requires a true on-shell averaged stress closure.
"""
import unittest
from math import exp
import numpy as np
from sdmc.pg_adiabatic_dilution import (
    WellEra,coefficients,well_x_acceleration,mode_averaged_terms,
    on_shell_analytic_energy_x_derivative,exact_reduced_ward_identity,
    adiabatic_formula,wkb_equation_of_state,wkb_exchange_over_H_E,
    wkb_log_dilution_slope,curvature_dust_running_condition,
    evolve_two_phase_ensemble,demo)

class StructuralAdiabaticDilution(unittest.TestCase):
    def test_on_shell_jordan_ward_identity_arbitrary_mode_state(self):
        for e in [
            WellEra(m2=.3,G0=.8,Lambda0=.4,k_run=.3,w_run=.5),
            WellEra(m2=0.,G0=1.,Lambda0=0.),
            WellEra(m2=0.,G0=0.,Lambda0=1.,p_v=2.,h=2.),
            WellEra(m2=2.,G0=.2,Lambda0=.5,k_run=-.4,w_run=3.)
        ]:
            for x in [0.,.25,.7,1.1]:
                for A,Ax in [(.1,.8),(.9,-2.),(-.4,.17)]:
                    r=exact_reduced_ward_identity(e,x,A,Ax)
                    self.assertLess(abs(r["residual"]),
                        1e-11*max(1.,r["terms"]["E"]))

    def test_direct_finite_difference_of_energy_matches_action_rate(self):
        e=WellEra(m2=.4,G0=.3,Lambda0=.7,
                   k_run=.2,w_run=.4,H0=.3)
        x=.5;A=.41;Ax=-.33
        Axx=well_x_acceleration(e,x,A,Ax)
        h=1e-5
        def local_energy(sign):
            xp=x+sign*h
            ap=A+sign*h*Ax+.5*(h**2)*Axx
            axp=Ax+sign*h*Axx
            return mode_averaged_terms(e,xp,ap,axp)["E"]
        finite=(local_energy(1)-local_energy(-1))/(2*h)
        exact=on_shell_analytic_energy_x_derivative(e,x,A,Ax)
        self.assertAlmostEqual(finite,exact,delta=1e-7*max(1.,abs(exact)))

    def test_wkb_pressure_and_exchange_match_adiabatic_log_slope(self):
        e=WellEra(m2=.2,G0=.4,Lambda0=.8,k_run=.15,w_run=4.,
                  h=1.5,p_v=1.5)
        for x in [.0,.4,.9]:
            p=wkb_equation_of_state(e,x)
            j=wkb_exchange_over_H_E(e,x)
            slope=wkb_log_dilution_slope(e,x)
            self.assertAlmostEqual(slope,-3*(1+p)+j,places=12)

    def test_mass_dominated_action_mode_is_dustlike_only_in_WKB(self):
        e=WellEra(m2=1.,G0=0.,Lambda0=0.,H0=.003)
        d=evolve_two_phase_ensemble(e,x_end=1.,num=150)
        self.assertAlmostEqual(d["E_ratio"][-1],exp(-3),delta=.015*exp(-3))
        self.assertLess(abs(d["P_over_E"][-1]),.02)
        self.assertEqual(d["J_over_H_E"][-1],0.)
        self.assertLess(max(d["ward_max_error"]),1e-10)

    def test_gradient_dominated_structural_mode_is_radiationlike(self):
        e=WellEra(m2=0.,G0=1.,Lambda0=0.,H0=.003)
        d=evolve_two_phase_ensemble(e,x_end=1.,num=150)
        self.assertAlmostEqual(d["E_ratio"][-1],exp(-4),delta=.02*exp(-4))
        self.assertAlmostEqual(d["P_over_E"][-1],1/3,delta=.025)
        self.assertAlmostEqual(d["J_over_H_E"][-1],0.)
        self.assertLess(max(d["ward_max_error"]),1e-10)

    def test_curvature_mode_in_matter_era_dilutes_as_a_minus_six(self):
        e=WellEra(m2=0.,G0=0.,Lambda0=1.,H0=.003,p_v=1.5,h=1.5)
        d=evolve_two_phase_ensemble(e,x_end=1.,num=150)
        self.assertAlmostEqual(d["E_ratio"][-1],exp(-6),delta=.07*exp(-6))
        self.assertLess(abs(d["P_over_E"][-1]),.05)
        self.assertAlmostEqual(d["J_over_H_E"][-1],-3.,delta=.15)
        self.assertLess(max(d["ward_max_error"]),1e-10)

    def test_curvature_mode_in_radiation_era_dilutes_as_a_minus_seven(self):
        e=WellEra(m2=0.,G0=0.,Lambda0=1.,H0=.001,
                  p_v=2.,h=2.)
        d=evolve_two_phase_ensemble(e,x_end=1.,num=150)
        self.assertAlmostEqual(d["E_ratio"][-1],exp(-7),delta=.07*exp(-7))
        self.assertLess(abs(d["P_over_E"][-1]),.04)
        self.assertAlmostEqual(d["J_over_H_E"][-1],-4.,delta=.2)

    def test_structural_curvature_dust_requires_running_coupling(self):
        e=WellEra(m2=0.,G0=0.,Lambda0=1.,
                  w_run=6.,H0=.003,p_v=1.5)
        d=evolve_two_phase_ensemble(e,x_end=1.,num=150)
        self.assertAlmostEqual(curvature_dust_running_condition(e),6.)
        self.assertAlmostEqual(d["E_ratio"][-1],exp(-3),delta=.02*exp(-3))
        self.assertAlmostEqual(d["J_over_H_E"][-1],0.)
        self.assertLess(abs(d["P_over_E"][-1]),.02)

    def test_kinetic_running_changes_coupling_condition(self):
        e=WellEra(m2=0.,G0=0.,Lambda0=1.,
                  k_run=2.,w_run=8.,H0=.003,p_v=1.5)
        self.assertAlmostEqual(curvature_dust_running_condition(e),8.)
        self.assertAlmostEqual(coefficients(e,.8)["omega"],1.)
        d=evolve_two_phase_ensemble(e,x_end=.8,num=150)
        self.assertAlmostEqual(d["E_ratio"][-1],exp(-2.4),
                               delta=.03*exp(-2.4))
        self.assertLess(abs(d["J_over_H_E"][-1]),.04)

    def test_general_mixed_mode_wkb_is_not_uniform_dust(self):
        e=WellEra(m2=.5,G0=1.,Lambda0=.2)
        d=evolve_two_phase_ensemble(e,x_end=.7,num=110)
        self.assertGreater(d["P_over_E"][-1],0.)
        self.assertLess(d["P_over_E"][-1],1/3)
        self.assertAlmostEqual(d["E_ratio"][-1],
                               d["WKB_E_ratio"][-1],
                               delta=.035*d["WKB_E_ratio"][-1])

    def test_p_zero_does_not_imply_j_zero_for_curvature(self):
        e=WellEra(m2=0.,G0=0.,Lambda0=1.)
        d=evolve_two_phase_ensemble(e,x_end=.6,num=130)
        self.assertLess(abs(d["P_over_E"][-1]),.03)
        self.assertLess(d["J_over_H_E"][-1],-2.8)

    def test_bare_rho_v_not_equal_reduced_curvature_energy(self):
        e=WellEra(m2=0.,G0=0.,Lambda0=1.,p_v=2.,h=2.)
        bare_rho_v_exponent=2*e.p_v
        geom_wkb_exponent=3+2*e.p_v
        self.assertEqual(bare_rho_v_exponent,4)
        self.assertEqual(geom_wkb_exponent,7)

    def test_wkb_has_admissibility_conditions(self):
        e=WellEra(m2=0.,G0=0.,Lambda0=1.,H0=2.)
        self.assertLess(coefficients(e,0.)["omega_over_H"],1.)
        # A low omega/H invalidates assuming the WKB result is accurate.
        self.assertEqual(wkb_log_dilution_slope(e,0.),-6.)

    def test_invalid_parameters_and_grid(self):
        with self.assertRaises(ValueError):
            WellEra(K0=-1.)
        with self.assertRaises(ValueError):
            WellEra(m2=0.,G0=0.,Lambda0=0.)
        with self.assertRaises(ValueError):
            evolve_two_phase_ensemble(WellEra(),x_end=-1.)

    def test_demo_separates_reduced_source_from_global_pg(self):
        r=demo()
        self.assertEqual(len(r["results"]),5)
        self.assertIn("global FLRW pressureless",str(r["not_derived"]))
        for a in r["results"]:
            self.assertGreater(a["final_omega_over_H"],15.)
            self.assertLess(a["exact_ward_max_error"],1e-10)

if __name__=="__main__":
    unittest.main()
