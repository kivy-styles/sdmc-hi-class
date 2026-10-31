"""Action-variation source regression for an excited SDMC structural well."""
import unittest
import numpy as np
from sdmc.pg_quadratic_source import TraceCoefficients
from sdmc.pg_inhomogeneous_adm_source import (
    fields,local_reduced_sources,lapse_shift_action,
    source_curvature_variation,density_linearization,
    restricted_spatial_pressure,canonical_gradient_traceless_stress,
    Friedmann_source_bookkeeping,spectral_cosine_amplitude,
    structural_source_demo)
from sdmc.pg_quadratic_source_3d import quadratic_energy_3d

class InhomogeneousADMSource(unittest.TestCase):
    def setUp(self):
        self.n=24
        self.L=2*np.pi
        a=np.arange(self.n)*self.L/self.n
        self.x,self.y,self.z=np.meshgrid(a,a,a,indexing="ij")
        self.e=.18*(np.cos(3*self.x)+.7*np.sin(4*self.y)+
                      .4*np.cos(5*self.z))
        self.v=.07*np.sin(3*self.x)+.05*np.cos(4*self.y)
        self.c=TraceCoefficients(bf2=3.,c_l2=.8,m_l2=.2,
                                  w_eta=.005,a_j=1.1,a_s=1.4,
                                  length=self.L)

    def test_lapse_uniform_variation_yields_negative_mean_reduced_energy(self):
        eps=1e-5
        p=np.ones_like(self.e)
        plus=lapse_shift_action(self.e,self.v,self.c,N=1+eps*p)
        minus=lapse_shift_action(self.e,self.v,self.c,N=1-eps*p)
        d=(plus-minus)/(2*eps)
        mean=float(local_reduced_sources(self.e,self.v,self.c)["rho"].mean())
        self.assertAlmostEqual(d,-mean,delta=2e-9*max(1,mean))

    def test_lapse_long_mode_variation_yields_density_fourier_source(self):
        b=np.cos(self.x)+.3*np.cos(2*self.y)
        eps=1e-5
        d=(lapse_shift_action(self.e,self.v,self.c,N=1+eps*b)-
           lapse_shift_action(self.e,self.v,self.c,N=1-eps*b))/(2*eps)
        expected=-np.mean(local_reduced_sources(self.e,self.v,self.c)["rho"]*b)
        self.assertAlmostEqual(d,float(expected),delta=2e-9*max(1,abs(expected)))

    def test_shift_variation_yields_canonical_momentum(self):
        amp=.03
        b=np.sin(6*self.x)+.2*np.sin(2*self.y)
        eps=1e-5
        plus=np.zeros((3,)+(self.n,)*3)
        plus[0]=eps*b
        l1=lapse_shift_action(self.e,self.v,self.c,shift=plus)
        l2=lapse_shift_action(self.e,self.v,self.c,shift=-plus)
        num=(l1-l2)/(2*eps)
        j=local_reduced_sources(self.e,self.v,self.c)["momentum"][0]
        self.assertAlmostEqual(num,-float(np.mean(j*b)),
                               delta=2e-9*max(1,abs(num)))

    def test_curvature_metric_variation_equals_bilaplacian_source(self):
        eps=1e-5
        z=np.cos(3*self.x)+.2*np.sin(4*self.y)
        num=(lapse_shift_action(self.e,self.v,self.c,zeta=eps*z)-
             lapse_shift_action(self.e,self.v,self.c,zeta=-eps*z))/(2*eps)
        source=source_curvature_variation(self.e,self.c)
        self.assertAlmostEqual(num,float(np.mean(source*z)),
                               delta=3e-8*max(1,abs(num)))

    def test_periodic_curvature_metric_source_has_zero_mean(self):
        s=source_curvature_variation(self.e,self.c)
        self.assertAlmostEqual(float(s.mean()),0.,places=10)
        self.assertGreater(float(np.sqrt(np.mean(s**2))),1e-5)

    def test_reduced_source_agrees_with_earlier_3d_energy(self):
        old=quadratic_energy_3d(self.e,self.c,ell_dot=self.v)
        new=local_reduced_sources(self.e,self.v,self.c)
        for term in ("kinetic","mass","gradient","curvature"):
            self.assertTrue(np.allclose(old[term],new[term],rtol=1e-12,atol=1e-12))
        self.assertTrue(np.allclose(new["rho"],sum(old.values())))

    def test_lapse_shift_source_is_quadratic_in_well_amplitude(self):
        n=local_reduced_sources(self.e,self.v,self.c)
        l=local_reduced_sources(3*self.e,3*self.v,self.c)
        self.assertTrue(np.allclose(l["rho"],9*n["rho"]))
        self.assertTrue(np.allclose(l["momentum"],9*n["momentum"]))

    def test_curvature_metric_source_is_linear_in_well_amplitude(self):
        source=source_curvature_variation(self.e,self.c)
        src_new=source_curvature_variation(3*self.e,self.c)
        self.assertTrue(np.allclose(src_new,3*source,atol=1e-11))

    def test_first_density_variation_matches_direct_finite_difference(self):
        de=.09*np.cos(2*self.x)+.05*np.sin(3*self.z)
        dv=.04*np.cos(3*self.x)+.03*np.sin(self.y)
        eps=1e-4
        plus=local_reduced_sources(self.e+eps*de,self.v+eps*dv,self.c)["rho"]
        minus=local_reduced_sources(self.e-eps*de,self.v-eps*dv,self.c)["rho"]
        estimate=(plus-minus)/(2*eps)
        exact=density_linearization(self.e,self.v,de,dv,self.c)
        self.assertLess(float(np.max(np.abs(estimate-exact))),2e-9)

    def test_zero_structural_background_has_no_first_density_variation(self):
        zero=np.zeros_like(self.e)
        first=density_linearization(zero,zero,self.e,self.v,self.c)
        self.assertAlmostEqual(float(np.max(np.abs(first))),0.,places=14)
        second=local_reduced_sources(self.e,self.v,self.c)["rho"]
        self.assertGreater(float(second.mean()),0.)

    def test_static_well_momentum_zero_but_density_positive_and_pressure_negative(self):
        z=np.zeros_like(self.e)
        val=local_reduced_sources(self.e,z,self.c)
        p=restricted_spatial_pressure(self.e,z,self.c)
        self.assertAlmostEqual(float(np.max(np.abs(val["momentum"]))),0.)
        self.assertGreater(float(val["rho"].mean()),0.)
        self.assertLess(float(p.mean()),0.)

    def test_standing_wave_can_have_zero_net_but_nonzero_local_momentum(self):
        e=np.cos(4*self.x)
        vel=.5*np.cos(4*self.x)
        val=local_reduced_sources(e,vel,self.c)
        self.assertLess(abs(float(val["momentum"][0].mean())),1e-11)
        self.assertGreater(float(np.mean(val["momentum"][0]**2)),.01)

    def test_traveling_wave_can_have_nonzero_mean_momentum(self):
        e=np.cos(4*self.x)
        vel=4*np.sin(4*self.x)
        val=local_reduced_sources(e,vel,self.c)
        self.assertAlmostEqual(float(val["momentum"][0].mean()),
                               -8*self.c.K,places=10)

    def test_gradient_stress_is_traceless_and_locally_anisotropic(self):
        pi=canonical_gradient_traceless_stress(self.e,self.c)
        tr=pi[0,0]+pi[1,1]+pi[2,2]
        self.assertLess(float(np.max(np.abs(tr))),1e-12)
        self.assertGreater(float(np.mean(pi[0,0]**2)),1e-5)
        self.assertTrue(np.allclose(pi[0,1],pi[1,0]))

    def test_equal_axis_modes_have_zero_averaged_gradient_anisotropy(self):
        e=(np.cos(4*self.x)+np.cos(4*self.y)+np.cos(4*self.z))/np.sqrt(3)
        pi=canonical_gradient_traceless_stress(e,self.c)
        self.assertLess(float(np.linalg.norm(pi.mean(axis=(2,3,4)))),1e-11)
        self.assertGreater(float(np.linalg.norm(pi)),1e-5)

    def test_perturbed_well_energy_has_a_long_mode(self):
        eps=.06;k=5
        short=(np.cos(k*self.x)+np.cos(k*self.y)+np.cos(k*self.z))/np.sqrt(3)
        ell=(1+eps*np.cos(self.x))*short
        val=local_reduced_sources(ell,np.zeros_like(ell),self.c)["rho"]
        amp=spectral_cosine_amplitude(val,self.L,k_long=1)
        self.assertGreater(amp,0.)
        self.assertAlmostEqual(amp/float(val.mean()),2*eps,delta=.02)

    def test_global_Friedmann_background_requires_a_new_constraint(self):
        val=local_reduced_sources(self.e,self.v,self.c)["rho"]
        rho_g=float(val.mean())
        b=Friedmann_source_bookkeeping(
            F=1.,H_j=1.,rho_b=2.,rho_r=.5,rho_X=.5,
            reduced_rho_g=rho_g)
        self.assertAlmostEqual(b["residual"],-rho_g,places=12)
        self.assertFalse(b["physical_solution_verified"])
        self.assertEqual(b["existing_source"],3.)

    def test_restricted_isotropic_pressure_from_action_scale_variation(self):
        a=self.c.a_j
        eps=1e-5
        def scale_action(s):
            c=TraceCoefficients(
                bf2=self.c.bf2,c_l2=self.c.c_l2,
                m_l2=self.c.m_l2,w_eta=self.c.w_eta,
                a_j=s,a_s=self.c.a_s,length=self.L)
            return s**3*lapse_shift_action(self.e,self.v,c)
        num=(scale_action(a+eps)-scale_action(a-eps))/(6*a*a*eps)
        expected=float(np.mean(restricted_spatial_pressure(
                              self.e,self.v,self.c)))
        self.assertAlmostEqual(num,expected,delta=2e-8*max(1,abs(expected)))

    def test_demo_reports_incomplete_averaged_stress(self):
        r=structural_source_demo()
        self.assertGreater(r["mean_reduced_density"],0)
        self.assertGreater(r["momentum_rms"],0)
        self.assertGreater(r["long_density_contrast_amplitude"],0)
        self.assertAlmostEqual(r["mean_trace"],0.,places=12)
        self.assertFalse(r["mean_sources_nonzero_but_scalar_only_background_unmodified"])
        self.assertIn("complete cubic Jordan metric stress",str(r["missing"]))

    def test_invalid_lapse_and_shift(self):
        with self.assertRaises(ValueError):
            lapse_shift_action(self.e,self.v,self.c,N=-np.ones_like(self.e))
        with self.assertRaises(ValueError):
            lapse_shift_action(self.e,self.v,self.c,shift=np.ones((3,12)))
        with self.assertRaises(ValueError):
            Friedmann_source_bookkeeping(-1,1,1,1,1,1)

if __name__=="__main__":
    unittest.main()
