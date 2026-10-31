"""Action-derived Jordan energy exchange for the restricted SDMC trace sector."""
import unittest
from math import exp,sqrt
import numpy as np
from sdmc.pg_energy_exchange import (
    Epoch, mode_parameters, mode_kinematics,
    reduced_energy_and_pressure, energy_rate_on_shell,
    jordan_energy_exchange_residual, jordan_exchange_analytic,
    required_driver_rate_for_zero_exchange, exchange_demo)

class PGEnergyExchange(unittest.TestCase):
    def setUp(self):
        self.e=Epoch(a_j=.92,a_s=1.34,H_j=.045,H_s=.0675,
               K=1.6,c2=.7,mass2=.35,weta=.08,
               Kdot=.035,c2dot=-.011,mass2dot=.019,wetadot=.004)
        self.k=2.3

    def test_jordan_structural_clocks_are_not_identical(self):
        self.assertAlmostEqual(self.e.p_v,1.5)
        self.assertNotEqual(self.e.H_s,self.e.H_j)

    def test_energy_pressure_and_eom_consistent(self):
        e=self.e
        ell,vel,u,udot=.31,-.13,.19,.047
        t=reduced_energy_and_pressure(ell,vel,u,self.k,e)
        self.assertGreater(t["E"],0)
        self.assertAlmostEqual(
            t["E"],sum(t[n] for n in ("kinetic","mass","gradient","curvature")))
        chi=mode_kinematics(ell,vel,u,udot,self.k,e)
        self.assertAlmostEqual(chi["q"],ell/3+u)
        self.assertAlmostEqual(chi["q_dot"],vel/3+udot)

    def test_energy_rate_matches_finite_difference_with_varying_K_coeff(self):
        e=self.e; k=self.k; ell,vel,u,udot=.31,-.13,.19,.047
        acc=mode_kinematics(ell,vel,u,udot,k,e)["ell_ddot"]
        def epo(t):
            return Epoch(a_j=e.a_j*exp(e.H_j*t),
                a_s=e.a_s*exp(e.H_s*t),
                H_j=e.H_j,H_s=e.H_s,K=e.K+e.Kdot*t,
                c2=e.c2+e.c2dot*t,
                mass2=e.mass2+e.mass2dot*t,
                weta=e.weta+e.wetadot*t,
                Kdot=e.Kdot,c2dot=e.c2dot,
                mass2dot=e.mass2dot,wetadot=e.wetadot)
        def calc(t):
            return reduced_energy_and_pressure(
                ell+vel*t+acc*t*t/2, vel+acc*t,
                u+udot*t,k,epo(t))["E"]
        h=2e-5
        numeric=(calc(h)-calc(-h))/(2*h)
        analytic=energy_rate_on_shell(ell,vel,u,udot,k,e)
        self.assertAlmostEqual(numeric,analytic,
            delta=1e-7*abs(analytic)+1e-8)

    def test_exact_energy_exchange_residual_matches_reduced_ward_identity(self):
        for ell,vel,u,udot in [
                (.31,-.13,.19,.047),(0.,0.,.3,0.),
                (1.,0.,0.,0.),(0.,1.,0.,0.),
                (.8,.5,.2,-.1)]:
            residual=jordan_energy_exchange_residual(
                ell,vel,u,udot,self.k,self.e)
            derived=jordan_exchange_analytic(
                ell,vel,u,udot,self.k,self.e)
            self.assertAlmostEqual(
                residual,derived,
                delta=1e-10*max(abs(residual),abs(derived),1))

    def test_fixed_coeff_no_curvature_trace_energy_separately_conserved(self):
        e=Epoch(1.,1.,.05,.075,1.,.7,.2,0.)
        for ell,vel in [(1.,0.),(0.,1.),(.6,.5)]:
            self.assertAlmostEqual(jordan_energy_exchange_residual(
                ell,vel,0.,0.,2.,e),0.,places=13)

    def test_curvature_trace_sector_exchanges_energy_when_S_expands(self):
        e=Epoch(1.,1.,.05,.075,1.,.7,.2,.04)
        ell=0.;vel=0.;u=.3;k=2.
        q=reduced_energy_and_pressure(ell,vel,u,k,e)
        residual=jordan_energy_exchange_residual(
            ell,vel,u,0.,k,e)
        self.assertAlmostEqual(
            residual,-4*e.H_s*q["curvature"],places=12)
        self.assertAlmostEqual(
            residual/(e.H_j*q["E"]),-4*e.p_v,places=12)

    def test_zero_restricted_pressure_does_not_imply_Jred_zero(self):
        e=Epoch(1.,1.,.05,.075,1.,.7,.2,.04)
        ell=.3;u=.2;k=2.
        q=reduced_energy_and_pressure(ell,0.,u,k,e)
        vel=sqrt(2*(q["mass"]+q["gradient"]/3+q["curvature"])/e.K)
        pressure=reduced_energy_and_pressure(ell,vel,u,k,e)["P_restricted"]
        self.assertAlmostEqual(pressure,0.,places=12)
        residual=jordan_energy_exchange_residual(
            ell,vel,u,0.,k,e)
        self.assertNotAlmostEqual(residual,0.,places=8)

    def test_coefficient_running_contributes_exchange_independent_of_H(self):
        e=Epoch(1.,1.,0.,0.,1.5,.8,.6,0.,
                Kdot=.1,c2dot=.08,mass2dot=.04)
        ell=.6;vel=-.2
        a=jordan_exchange_analytic(ell,vel,0.,0.,2.,e)
        b=jordan_energy_exchange_residual(ell,vel,0.,0.,2.,e)
        self.assertAlmostEqual(a,b,places=12)
        self.assertNotAlmostEqual(a,0.,places=8)

    def test_driver_balance_is_condition_not_dynamical_derivation(self):
        e=self.e;ell=.5;vel=-.3;u=.2;k=3.
        udot=required_driver_rate_for_zero_exchange(
            ell,vel,u,k,e)
        residual=jordan_energy_exchange_residual(
            ell,vel,u,udot,k,e)
        self.assertAlmostEqual(residual,0.,places=10)
        with self.assertRaises(ValueError):
            required_driver_rate_for_zero_exchange(0,0,0,0,e)

    def test_demo_has_distinct_structural_exchange(self):
        r=exchange_demo()
        self.assertAlmostEqual(r["structural_tracker_p_v"],1.5)
        self.assertAlmostEqual(r["unforced_exchange_per_H_E"],-6.)
        self.assertLess(abs(r["balanced_exchange_residual"]),1e-12)
        self.assertLess(r["curvature_only_restricted_pressure"],0)

if __name__=="__main__":
    unittest.main()
