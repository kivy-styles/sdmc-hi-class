"""Nonlinear conformal supporting-space curvature and cubic metric variation."""
import unittest
import numpy as np
from sdmc.pg_quadratic_source import TraceCoefficients
from sdmc.pg_quadratic_source_3d import geometry_derivatives
from sdmc.pg_nonlinear_conformal_stress import (
    exact_conformal_riemann2,exact_curvature_action,
    quadratic_curvature_action,through_cubic_curvature_action,
    cubic_spatial_metric_source,quadratic_spatial_metric_source,
    conformal_demo)

class NonlinearConformalStress(unittest.TestCase):
    def setUp(self):
        self.n=24
        self.L=2*np.pi
        a=np.arange(self.n)*self.L/self.n
        self.x,self.y,self.z=np.meshgrid(a,a,a,indexing="ij")
        self.ell=(.13*np.cos(3*self.x)+.09*np.sin(4*self.y)+
                  .035*np.cos(5*self.z))
        self.zero=np.zeros_like(self.ell)
        self.c=TraceCoefficients(bf2=3.,c_l2=1.,
              m_l2=.2,w_eta=.01,a_j=1.,a_s=1.2,
              length=self.L)

    def test_exact_3d_riemann_squared_matches_independent_ricci_identity(self):
        q=self.ell/3
        grad,H,L=geometry_derivatives(q,self.L)
        w=sum(x*x for x in grad)
        Rij=np.asarray([[ -H[i][j]-(L+w if i==j else 0.)
                          +grad[i]*grad[j] for j in range(3)]
                          for i in range(3)])
        tr=sum(Rij[i,i] for i in range(3))
        ric_sq=sum(Rij[i,j]**2 for i in range(3) for j in range(3))
        direct=np.exp(-4*q)/self.c.a_s**4*(4*ric_sq-tr**2)
        exact=exact_conformal_riemann2(q,self.c)
        self.assertTrue(np.allclose(exact,direct,rtol=1e-10,atol=1e-10))
        self.assertGreater(float(exact.mean()),0.)

    def test_quadratic_curvature_recovers_previous_structural_energy(self):
        q=.001*self.ell
        exact=-exact_curvature_action(q,self.zero,self.c)
        quad=-quadratic_curvature_action(q,self.zero,self.c)
        self.assertAlmostEqual(exact,quad,delta=.01*quad)

    def test_through_cubic_action_is_fourth_order_accurate(self):
        a=self.ell
        def residual(f):
            x=f*a
            return abs(exact_curvature_action(x,self.zero,self.c)-
                     through_cubic_curvature_action(x,self.zero,self.c))
        a1=residual(.15)
        a2=residual(.075)
        self.assertGreater(a1,a2)
        self.assertAlmostEqual(a1/a2,16.,delta=2.)

    def test_quadratic_spatial_curvature_metric_source_finite_difference(self):
        z=np.cos(3*self.x)+.3*np.sin(4*self.y)
        h=1e-5
        numerical=(quadratic_curvature_action(self.ell,h*z,self.c)-
                   quadratic_curvature_action(self.ell,-h*z,self.c))/(2*h)
        predicted=float(np.mean(
              quadratic_spatial_metric_source(self.ell,self.c)*z))
        self.assertAlmostEqual(numerical,predicted,
                               delta=1e-8*max(1,abs(predicted)))

    def test_cubic_spatial_metric_source_matches_action_variation(self):
        z=np.cos(2*self.x)+.15*np.sin(3*self.y)
        h=1e-5
        def cubic_part(s):
            return (through_cubic_curvature_action(self.ell,s*z,self.c)-
                   quadratic_curvature_action(self.ell,s*z,self.c))
        numerical=(cubic_part(h)-cubic_part(-h))/(2*h)
        predicted=float(np.mean(
              cubic_spatial_metric_source(self.ell,self.c)*z))
        self.assertAlmostEqual(numerical,predicted,
                               delta=2e-8*max(1,abs(predicted)))

    def test_mean_cubic_metric_source_equals_quadratic_curvature_energy(self):
        r=cubic_spatial_metric_source(self.ell,self.c)
        quad=-quadratic_curvature_action(self.ell,self.zero,self.c)
        self.assertAlmostEqual(float(r.mean()),quad,
                               delta=3e-11*max(1,abs(quad)))

    def test_second_order_source_scales_quadratically_with_well(self):
        r=cubic_spatial_metric_source(self.ell,self.c)
        scaled=cubic_spatial_metric_source(2*self.ell,self.c)
        self.assertTrue(np.allclose(scaled,4*r,rtol=1e-11,atol=1e-10))

    def test_linear_curvature_source_has_zero_mean_on_periodic_domain(self):
        r=quadratic_spatial_metric_source(self.ell,self.c)
        self.assertAlmostEqual(float(r.mean()),0.,places=12)
        self.assertGreater(float(np.sqrt(np.mean(r*r))),1e-6)

    def test_full_nonlinear_uniform_zeta_variation_has_exact_volume_scaling(self):
        # q=e+zeta: a constant ζ leaves its derivatives unchanged,
        # but sqrt(h) Riem² scales as exp(3ζ-4ζ)=exp(-ζ).
        h=1e-4
        a=exact_curvature_action(self.ell,self.zero,self.c)
        plus=exact_curvature_action(self.ell,h*np.ones_like(self.ell),self.c)
        minus=exact_curvature_action(self.ell,-h*np.ones_like(self.ell),self.c)
        self.assertAlmostEqual(plus,np.exp(-h)*a,delta=3e-12)
        self.assertAlmostEqual(minus,np.exp(h)*a,delta=3e-12)
        deriv=(plus-minus)/(2*h)
        self.assertAlmostEqual(deriv,-a,delta=1e-7*max(1,abs(a)))

    def test_cubic_terms_explain_first_nonlinear_metric_source(self):
        z=np.cos(3*self.x)
        h=1e-5
        scale=.02
        ell=scale*self.ell
        exact=(exact_curvature_action(ell,h*z,self.c)-
               exact_curvature_action(ell,-h*z,self.c))/(2*h)
        model=float(np.mean(
           (quadratic_spatial_metric_source(ell,self.c)+
            cubic_spatial_metric_source(ell,self.c))*z))
        self.assertAlmostEqual(exact,model,delta=2e-7*max(1,abs(model)))

    def test_curvature_source_is_zero_for_constant_structural_trace(self):
        q=np.ones_like(self.ell)*.15
        r=exact_conformal_riemann2(q,self.c)
        self.assertLess(float(np.max(np.abs(r))),1e-10)
        self.assertLess(abs(exact_curvature_action(3*q,self.zero,self.c)),1e-10)

    def test_demo_does_not_claim_physical_full_stress(self):
        r=conformal_demo()
        self.assertGreater(r["mean_quadratic_curvature_energy"],0.)
        self.assertGreater(r["exact_riemann2_mean"],0.)
        self.assertLess(abs(r["cubic_vs_quadratic_mean_identity_residual"]),1e-10)
        self.assertIn("full nonlinear trace/ADM kinetic",str(r["missing"]))

if __name__=="__main__":
    unittest.main()
