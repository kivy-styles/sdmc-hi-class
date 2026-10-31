"""Nonlinear/averaging checks for SDMC's reduced structural-well energy.

These tests verify algebra and numerical coarse-graining; they DO NOT
establish that the resulting Jordan averaged energy behaves as CDM.
"""
import unittest
from math import exp,pi
import numpy as np
from sdmc.pg_linear_perturbations import trace_dynamic_derivatives
from sdmc.pg_quadratic_source import (
    TraceCoefficients, derivative, gaussian_coarse_grain,
    quadratic_energy_terms, quadratic_energy_density,
    first_variation_energy, averaged_energy, coarse_grained_contrast,
    energy_mode, mode_energy_rate, harmonic_demo,
    restricted_isotropic_pressure, restricted_isotropic_pressure_mean,
    integrated_static_kernel, integrated_static_energy
)

class GeometricQuadraticSource(unittest.TestCase):
    def setUp(self):
        self.n=1024
        self.L=2*pi
        self.x=np.arange(self.n)*self.L/self.n
        self.c=TraceCoefficients(bf2=3.,c_l2=1.3,m_l2=.9,
              w_eta=.25,a_j=.87,a_s=1.25,length=self.L)

    def test_spectral_derivatives(self):
        q=np.cos(19*self.x)
        self.assertLess(np.max(np.abs(derivative(q,self.L,1)+
                                 19*np.sin(19*self.x))),1e-10)
        self.assertLess(np.max(np.abs(derivative(q,self.L,2)+
                                 19**2*q)),1e-8)
        self.assertLess(np.max(np.abs(derivative(q,self.L,4)-
                                 19**4*q))/19**4,1e-8)

    def test_quadratic_energy_positive_and_zero_limit(self):
        q=np.cos(11*self.x)
        terms=quadratic_energy_terms(q,self.c,ell_dot=.2*q)
        for v in terms.values():
            self.assertTrue(np.all(v>=0))
        self.assertGreater(sum(v.mean() for v in terms.values()),0)
        z=np.zeros(self.n)
        self.assertEqual(averaged_energy(z,self.c),0.)
        with self.assertRaises(ValueError):
            coarse_grained_contrast(z,self.c,width=.15)

    def test_first_order_source_zero_on_locked_flat_branch(self):
        z=np.zeros(self.n)
        pert=np.sin(9*self.x)
        val=first_variation_energy(z,pert,self.c,
                       dell_dot=pert,dm=.1*pert,dphi=.2*pert)
        self.assertEqual(float(np.max(np.abs(val))),0.)

    def test_variation_matches_direct_finite_difference(self):
        ell=.3*np.cos(11*self.x)+.2*np.sin(13*self.x)
        vel=.04*np.sin(17*self.x)
        mm=.02*np.sin(3*self.x)
        phi=.03*np.cos(2*self.x)
        de=.05*np.cos(5*self.x)
        dv=.02*np.sin(7*self.x)
        dm=.01*np.cos(2*self.x)
        dp=.04*np.sin(3*self.x)
        h=1e-3
        up=quadratic_energy_density(ell+h*de,self.c,
                        vel+h*dv,mm+h*dm,phi+h*dp)
        dn=quadratic_energy_density(ell-h*de,self.c,
                        vel-h*dv,mm-h*dm,phi-h*dp)
        exact=first_variation_energy(ell,de,self.c,vel,dv,
                                     mm,dm,phi,dp)
        scale=max(float(np.max(np.abs(exact))),1.)
        self.assertLess(float(np.max(np.abs((up-dn)/(2*h)-exact)))/scale,
                        2e-8)

    def test_gaussian_averaging_preserves_total_mean(self):
        q=np.cos(7*self.x)+.2*np.sin(2*self.x)
        avg=gaussian_coarse_grain(q,self.L,.2)
        self.assertAlmostEqual(float(avg.mean()),float(q.mean()),places=14)
        self.assertAlmostEqual(float(np.max(np.abs(
            gaussian_coarse_grain(q,self.L,0)-q))),0.)

    def test_short_well_envelope_generates_long_source(self):
        ks,kl,eps=32,2,.1
        ell=(1+eps*np.cos(kl*self.x))*np.cos(ks*self.x)
        coeff=TraceCoefficients(bf2=3.,c_l2=1.,m_l2=.05,
                   w_eta=.01,a_j=1.,a_s=1.,length=self.L)
        width=.16
        contrast=coarse_grained_contrast(ell,coeff,width)
        long_coeff=2*np.mean(contrast*np.cos(kl*self.x))
        expected=2*eps*np.exp(-.5*(kl*width)**2)
        self.assertLess(abs(long_coeff-expected)/expected,.03)
        self.assertAlmostEqual(float(np.mean(contrast)),0.,places=12)

    def test_long_mode_disappears_without_amplitude_modulation(self):
        ks,kl=32,2
        ell=np.cos(ks*self.x)
        coeff=TraceCoefficients(length=self.L)
        d=coarse_grained_contrast(ell,coeff,width=.16)
        long_coeff=2*np.mean(d*np.cos(kl*self.x))
        self.assertAlmostEqual(float(long_coeff),0.,places=12)

    def test_dilution_of_individual_terms_not_dust_by_default(self):
        ell=np.cos(8*self.x)
        p=1.5
        c1=TraceCoefficients(bf2=2,c_l2=.4,m_l2=.8,
                             w_eta=.1,a_j=1,a_s=1,length=self.L)
        c2=TraceCoefficients(bf2=2,c_l2=.4,m_l2=.8,
                             w_eta=.1,a_j=2,a_s=2**p,length=self.L)
        t1=quadratic_energy_terms(ell,c1)
        t2=quadratic_energy_terms(ell,c2)
        self.assertAlmostEqual(t2['mass'].mean()/t1['mass'].mean(),1.)
        self.assertAlmostEqual(t2['gradient'].mean()/t1['gradient'].mean(),2**-2)
        self.assertAlmostEqual(t2['curvature'].mean()/t1['curvature'].mean(),2**(-4*p))

    def test_mode_energy_balance_matches_dynamical_trace_eom(self):
        c=TraceCoefficients(bf2=3.,c_l2=.8,m_l2=.4,
                       w_eta=.02,a_j=.91,a_s=1.12,length=self.L)
        ell=.35; vel=-.11; u=.17; udot=.03
        h=.04; hs=.06; k=1.5
        _,acc=trace_dynamic_derivatives(ell,vel,h,0.,-u,k,
                       c.a_j,c.a_s,c.bf2,c.c_l2,c.m_l2,c.w_eta)
        predicted=mode_energy_rate(ell,vel,u,udot,k,c,h,hs)
        step=2e-6
        def en(t):
            ct=TraceCoefficients(bf2=c.bf2,c_l2=c.c_l2,
                 m_l2=c.m_l2,w_eta=c.w_eta,
                 a_j=c.a_j*exp(h*t),a_s=c.a_s*exp(hs*t),length=c.length)
            return energy_mode(ell+vel*t+.5*acc*t*t,vel+acc*t,
                              u+udot*t,k,ct)
        measured=(en(step)-en(-step))/(2*step)
        self.assertAlmostEqual(measured,predicted,delta=abs(predicted)*1e-7+1e-8)

    def test_unforced_energy_is_nonincreasing_in_expanding_epoch(self):
        c=self.c
        for ell,vel in [(1,0),(0,1),(.7,.6)]:
            ans=mode_energy_rate(ell,vel,0.,0.,k=3.,coeff=c,
                                 H_j=.1,H_s=.15)
            self.assertLessEqual(ans,0.)

    def test_sourced_well_can_exchange_energy(self):
        c=self.c
        src=mode_energy_rate(.4,0.,.2,20.,2.,c,.1,.15)
        free=mode_energy_rate(.4,0.,.2,0.,2.,c,.1,.15)
        self.assertGreater(src,free)


    def test_static_well_has_negative_restricted_isotropic_pressure(self):
        ell=np.cos(8*self.x)
        p=restricted_isotropic_pressure_mean(ell,self.c)
        rho=averaged_energy(ell,self.c)
        self.assertGreater(rho,0)
        self.assertLess(p,0)
        # Fixed-S, reduced isotropic metric variation ONLY.

    def test_restricted_pressure_equals_scale_derivative_of_action(self):
        ell=.4*np.cos(5*self.x)
        vel=.2*np.sin(3*self.x)
        a=self.c.a_j
        def averaged_action(a_value):
            c=TraceCoefficients(bf2=self.c.bf2,c_l2=self.c.c_l2,
              m_l2=self.c.m_l2,w_eta=self.c.w_eta,
              a_j=a_value,a_s=self.c.a_s,length=self.L)
            t=quadratic_energy_terms(ell,c,ell_dot=vel)
            lag=float(np.mean(t['kinetic']-t['mass']-
                              t['gradient']-t['curvature']))
            return a_value**3*lag
        h=1e-5
        scale_deriv=(averaged_action(a+h)-averaged_action(a-h))/(6*a*a*h)
        pressure=restricted_isotropic_pressure_mean(ell,self.c,ell_dot=vel)
        self.assertAlmostEqual(scale_deriv,pressure,
                               delta=1e-6*abs(pressure)+1e-9)

    def test_dust_pressure_requires_energy_balance_in_reduced_action(self):
        ell=np.cos(7*self.x)
        baseline=quadratic_energy_terms(ell,self.c)
        required=(baseline['mass'].mean()+
                  baseline['gradient'].mean()/3+
                  baseline['curvature'].mean())
        vel=np.full(self.n,np.sqrt(2*required/self.c.K))
        p=restricted_isotropic_pressure_mean(ell,self.c,ell_dot=vel)
        self.assertAlmostEqual(p,0.,delta=1e-10*required)

    def test_integrated_trace_kernel_matches_direct_minimization(self):
        c=self.c
        for k in [.01,.5,3.,8.]:
            u=.3
            w=c.m_l2+c.c_l2*k*k/c.a_j**2
            ck=c.C*k**4
            lam=2*ck/(9*c.K)
            ell=-3*lam/(w+lam)*u
            total=.5*c.K*w*ell**2+ck*(u+ell/3)**2
            reduced=integrated_static_energy(u,k,c)
            self.assertAlmostEqual(total,reduced,
                                   delta=abs(reduced)*1e-12+1e-12)
            self.assertGreaterEqual(integrated_static_kernel(k,c),0.)

    def test_curvature_metric_variation_is_not_a_density_prediction(self):
        c=self.c;k=3.;u=.4;eps=1e-5
        # u=m-Phi; variation w.r.t Phi gives -2*Keff*u.
        num=(integrated_static_energy(u-eps,k,c)-
             integrated_static_energy(u+eps,k,c))/(2*eps)
        self.assertAlmostEqual(num,-2*integrated_static_kernel(k,c)*u,
                               delta=1e-8*abs(num)+1e-9)
        self.assertEqual(integrated_static_kernel(0.,c),0.)

    def test_demo_is_finite_and_nonzero(self):
        d=harmonic_demo()
        self.assertGreater(d['mean_energy_unmodulated_arbitrary_units'],0)
        self.assertGreater(d['long_density_contrast_cosine_amplitude'],0)
        self.assertLess(abs(d['long_density_contrast_cosine_amplitude']
                            -d['leading_separation_approximation']),.01)

if __name__=="__main__":
    unittest.main()
