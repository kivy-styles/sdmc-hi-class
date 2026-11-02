"""SDMC material inventory / structural Fourier response test suite.

These validate only the conditional Newtonian matter source and the
reduced structural action, not full Planck-era physics or CDM replacement.
"""
import unittest
import numpy as np
from dataclasses import replace
from sdmc.pg_redistribution_spectrum import (
    SelectorParameters,response,source_energy,off_shell_energy,
    density_mode,effective_count,iid_mode_moments,
    correlated_positions,poisson_source,spectrum_response,
    activation_exponents,diagnostic)

class RedistributedStructuralSpectrum(unittest.TestCase):
    def setUp(self):
        self.p=SelectorParameters()
        self.rng=np.random.default_rng(27413)

    def test_inventory_is_exactly_conserved_under_partition(self):
        for rest in (0.,.2,2.,15.):
            masses=np.array([.5,3.,1.5,2.])
            Mu=rest+sum(masses)
            self.assertAlmostEqual(Mu,rest+masses.sum())
            a=density_mode([0.,1.,2.,3.],k=2,masses=masses,
                           core_mass=rest)
            b=density_mode([0.,1.,2.,3.],k=2,
                           masses=masses/Mu,core_mass=rest/Mu)
            self.assertAlmostEqual(abs(a-b),0.,places=13)

    def test_exact_random_mode_moments_for_equal_masses(self):
        for n in (10,100,1000,10000):
            r=iid_mode_moments(np.ones(n),0.)
            self.assertAlmostEqual(r["power"],1/n,places=14)
            self.assertAlmostEqual(r["effective_count"],n,places=11)

    def test_unequal_mass_effective_count_not_planck_mass_count(self):
        weights=np.array([1.,1.,1.,7.])
        r=iid_mode_moments(weights,0.)
        self.assertAlmostEqual(r["power"],52/100,places=12)
        self.assertAlmostEqual(r["effective_count"],100/52,places=12)

    def test_correlated_expected_mode_and_shot_noise(self):
        eps=.12
        r=iid_mode_moments(np.ones(4096),eps/2)
        self.assertAlmostEqual(r["mean_real"],.06)
        self.assertAlmostEqual(r["power"],.06**2+(1-.06**2)/4096)

    def test_random_simulation_recovers_power_without_mass_creation(self):
        n=4096
        amps=np.array([density_mode(correlated_positions(n,0.,self.rng),
                                      k=1) for _ in range(128)])
        theory=1/n
        sample=float(np.mean(abs(amps)**2))
        self.assertAlmostEqual(sample,theory,delta=.35*theory)
        self.assertLess(abs(amps.mean()),.003)

    def test_nonuniform_positions_create_coherent_mode(self):
        n=4096;eps=.12
        amps=np.array([density_mode(correlated_positions(n,eps,self.rng),
                                      k=1) for _ in range(128)])
        theory=iid_mode_moments(np.ones(n),eps/2)
        self.assertAlmostEqual(float(amps.mean().real),
                theory["mean_real"],delta=.003)
        self.assertAlmostEqual(float(np.mean(abs(amps)**2)),
                theory["power"],delta=.0004)

    def test_core_partition_changes_density_mode_not_total_mass(self):
        loc=[.1,1.,2.,3.]
        pure=density_mode(loc,k=1,masses=[1.,1.,1.,1.])
        core=density_mode(loc,k=1,masses=[.75]*4,core_mass=1.)
        expected=.75*pure+.25
        self.assertAlmostEqual(abs(core-expected),0.,places=12)

    def test_positive_selector_has_bounded_transmission(self):
        for k in (.05,.5,1,2,4,8):
            r=response(k,self.p)
            self.assertLessEqual(r["transmission"],1.)
            self.assertGreater(r["transmission"],0.)
            self.assertAlmostEqual(r["transmission"]+
                r["relative_fraction"],1.,places=14)

    def test_quasistatic_ell_solves_energy_stationarity(self):
        for k in (.5,1,2,4):
            u=.01
            s=source_energy(k,u,self.p)
            h=1e-6
            grad=(off_shell_energy(k,u,s["ell"]+h,self.p)
                 -off_shell_energy(k,u,s["ell"]-h,self.p))/(2*h)
            self.assertLess(abs(grad),2e-10)
            self.assertLessEqual(s["total"],
                  off_shell_energy(k,u,0.,self.p)+1e-15)

    def test_source_on_shell_curvature_and_relative_energy_sum(self):
        for k in (.25,1.,2.,4.):
            for u in (.1,.02,complex(.01,.03)):
                s=source_energy(k,u,self.p)
                r=response(k,self.p)
                self.assertAlmostEqual(
                    s["total"],s["relative"]+s["curvature"],
                    places=13)
                self.assertAlmostEqual(s["expected"],
                    r["C"]*r["transmission"]*abs(u)**2,places=13)
                self.assertAlmostEqual(s["relative"],
                    r["C"]*r["transmission"]*(1-r["transmission"])
                        *abs(u)**2,delta=1e-12)
                self.assertAlmostEqual(s["on_shell_residual"],0.,
                                       delta=1e-12)

    def test_poisson_k4_cancellation_at_fixed_delta(self):
        d=.02+1j*.01
        A=.005
        for k in (1.,2.,4.,8.):
            v=spectrum_response([d],[k],self.p,A)["modes"][0]
            predicted=4*self.p.W*self.p.eta/self.p.aS**4*A*A*\
                v["transmission"]*abs(d)**2
            self.assertAlmostEqual(v["energy"],predicted,places=13)
            self.assertAlmostEqual(
                v["energy"],v["predicted_cancelled_k4_energy"],places=13)

    def test_after_selector_off_forced_energy_and_ell_disappear(self):
        off=replace(self.p,W=0)
        for k in (1.,2.,4.):
            s=source_energy(k,.02,off)
            self.assertEqual(s["ell"],0.)
            self.assertEqual(s["total"],0.)
            self.assertAlmostEqual(s["q"],.02)
        spectra=spectrum_response([.1,.05],[1.,2.],off,.01)
        self.assertEqual(spectra["E_sum"],0.)

    def test_selector_not_independent_material_matter_density(self):
        s=spectrum_response([.05,.02],[1.,2.],self.p,.01)
        self.assertGreater(s["E_sum"],0)
        self.assertFalse(s["observational_gravity_proven"])

    def test_required_spectrum_growth_is_stronger_than_dust_growth(self):
        for p,powers in [(2,7),(1.5,5),(1,3)]:
            r=activation_exponents(p)
            self.assertAlmostEqual(r["required_power_slope_if_W_T_constant"],
                                   powers)
            self.assertAlmostEqual(r["required_delta_amplitude_slope"],
                                   powers/2)

    def test_actual_gamma_slope_includes_chi_transfer(self):
        # fixed W and frozen delta, conserved Jordan nonrelativistic rho
        # A_m~a_J^-1; rho_v~a_s^-2, and a_s=a_J^p.
        p_v=1.5
        h=1e-4
        def gamma(a):
            params=replace(self.p,aJ=a,aS=a**p_v)
            E=spectrum_response([.03],[2.],params,A_m=.02/a)["E_sum"]
            return E*a**(2*p_v)
        numerical=(np.log(gamma(np.exp(h)))-
                np.log(gamma(np.exp(-h))))/(2*h)
        params_p=replace(self.p,aJ=np.exp(h),aS=np.exp(h*p_v))
        params_m=replace(self.p,aJ=np.exp(-h),aS=np.exp(-h*p_v))
        dlnT=(np.log(response(2.,params_p)["transmission"])-
              np.log(response(2.,params_m)["transmission"]))/(2*h)
        expected=-2*p_v-2+dlnT
        self.assertAlmostEqual(numerical,expected,delta=2e-8)

    def test_global_K_zero_mode_not_a_poisson_source(self):
        with self.assertRaises(ValueError):
            response(0.,self.p)
        with self.assertRaises(ValueError):
            poisson_source(0.,1.,1.)

    def test_illegal_coefficients_rejected(self):
        with self.assertRaises(ValueError):
            SelectorParameters(B=0)
        with self.assertRaises(ValueError):
            SelectorParameters(W=-1)
        with self.assertRaises(ValueError):
            density_mode([.1,.2],1,masses=[1,-1])

    def test_demo_does_not_claim_physical_cdm_replacement(self):
        d=diagnostic()
        self.assertGreater(d["spectrum"]["E_sum"],0.)
        self.assertEqual(d["selector_off"]["E_sum"],0.)
        self.assertAlmostEqual(d["required_activation"]["radiation"]
                 ["required_power_slope_if_W_T_constant"],7)
        self.assertIn("full Jordan F/G2/G3",str(d["omitted"]))

if __name__=="__main__":
    unittest.main()
