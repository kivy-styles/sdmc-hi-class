"""Tests for SDMC pg action perturbations and conditional growth closure."""
import unittest
import numpy as np
from sdmc.pg_linear_perturbations import (
    Background, archived_reference, dust_conformal_derivatives,
    gamma_for_dust, uniform_nested_sphere_mass_ratio,
    quadratic_energy_variation, trace_static_response,
    solve_subhorizon_growth, psi_subhorizon
)

class PGLinearEquations(unittest.TestCase):
    def test_flat_locked_source_is_zero_at_linear_order(self):
        for dl,dld,dq in [(1.,0.,0.),(0.,3.,0.),(0.,0.,-2.),(.5,.4,.2)]:
            ans=quadratic_energy_variation(.5,1.,1.,2.,1.,.4,.8,
                ell0=0.,ell0_dot=0.,q0=0.,
                dell=dl,dell_dot=dld,dq=dq)
            self.assertEqual(ans,0.)

    def test_background_strain_can_source_first_order_energy(self):
        ans=quadratic_energy_variation(.4,1.,1.,2.,1.,.5,.4,
                ell0=.2,ell0_dot=.1,q0=.3,
                dell=.1,dell_dot=.2,dq=-.1)
        self.assertNotEqual(ans,0.)

    def test_static_curvature_stiffness_screens_q(self):
        weak=trace_static_response(.5,1.,1.,2.,1.,.3,.1,.2)
        stiff=trace_static_response(.5,1.,1.,2.,1.,.3,10.,.2)
        self.assertGreater(weak['q_over_original'],stiff['q_over_original'])
        for a in [weak,stiff]:
            self.assertAlmostEqual(a['q'],-.2*a['q_over_original'],places=13)

    def test_exact_dust_conformal_conservation(self):
        dd,dt=dust_conformal_derivatives(.2,.4,.5,-.1,.03,-.2)
        self.assertAlmostEqual(dd,.19)
        self.assertAlmostEqual(dt,.02-.032)

    def test_gamma_structural_radiation_and_matter_scaling(self):
        for a in [.001,.01,.2,1.]:
            self.assertAlmostEqual(gamma_for_dust(a,a**1.5,.25),.25,places=12)
            self.assertAlmostEqual(gamma_for_dust(a,a**2,.25),.25*a,places=12)

    def test_nested_sphere_is_special_not_general_volume_law(self):
        self.assertAlmostEqual(uniform_nested_sphere_mass_ratio(.94644),6.568,delta=.02)
        with self.assertRaises(ValueError):
            uniform_nested_sphere_mass_ratio(1.01)

    def test_eds_growth_growing_mode(self):
        bg=Background(.7,.05,.95,0.,0.)
        d=solve_subhorizon_growth(bg,z_initial=20.,n=100)
        self.assertLess(np.max(np.abs(d['delta_b']-d['a'])),2e-7)
        self.assertLess(np.max(np.abs(d['delta_pg']-d['a'])),2e-7)
        self.assertLess(np.max(np.abs(d['u_b']+d['a'])),2e-7)

    def test_effective_dust_is_exact_linear_cdm_control(self):
        d=solve_subhorizon_growth(archived_reference())
        self.assertTrue(np.allclose(d['delta_b'],d['delta_pg'],rtol=1e-12,atol=1e-12))
        self.assertGreater(d['delta_pg'][-1],.5)

    def test_smooth_geometric_source_fails_to_cluster(self):
        bg=archived_reference()
        cold=solve_subhorizon_growth(bg)
        smooth=solve_subhorizon_growth(bg,smooth_pg=True)
        self.assertGreater(cold['delta_b'][-1]/smooth['delta_b'][-1],8.)
        self.assertTrue(np.all(smooth['delta_pg']==0))

    def test_small_effective_pressure_suppresses_growth(self):
        bg=archived_reference()
        cold=solve_subhorizon_growth(bg,k_hmpc=.1)
        press=solve_subhorizon_growth(bg,k_hmpc=.1,cs2_pg=1e-4)
        self.assertLess(press['delta_pg'][-1],cold['delta_pg'][-1]*.1)

    def test_subhorizon_potential_sign(self):
        self.assertLess(psi_subhorizon(archived_reference(),1.,.1,.8,.8),0)

if __name__=="__main__":
    unittest.main()
