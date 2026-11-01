"""Excited traceless supporting-space ADM lapse/shift source validation."""
import unittest
import numpy as np
from sdmc.pg_full_spatial_stress import SupportCoefficients
from sdmc.pg_tensor_adm_constraints import (
    fields,local_3d_energy,full_trace_tensor_shift_momentum,
    action_3d_lapse_shift,demo)

class TensorADMSource(unittest.TestCase):
    def setUp(self):
        self.n=20; self.L=2*np.pi
        a=np.arange(self.n)*self.L/self.n
        self.x,self.y,self.z=np.meshgrid(a,a,a,indexing="ij")
        self.ell=.06*(np.cos(2*self.x)+.6*np.sin(3*self.y))
        self.vel=.025*np.sin(2*self.x)+.01*np.cos(3*self.y)
        self.Q=np.zeros((3,3,self.n,self.n,self.n))
        self.dQ=self.Q.copy()
        self.Q[0,0]=.03*np.cos(2*self.x)
        self.Q[1,1]=-self.Q[0,0]
        self.dQ[0,1]=.025*np.sin(3*self.y)
        self.dQ[1,0]=self.dQ[0,1]
        self.c=SupportCoefficients(
             bf2=3.,c2=.7,m2=.2,weta=.012,
             a_j=1.2,a_s=1.4,length=self.L)

    def test_lapse_uniform_derivative_equals_negative_mean_energy(self):
        h=1e-5
        ones=np.ones_like(self.ell)
        plus=action_3d_lapse_shift(
            self.ell,self.vel,self.Q,self.dQ,self.c,lapse=1+h*ones)
        minus=action_3d_lapse_shift(
            self.ell,self.vel,self.Q,self.dQ,self.c,lapse=1-h*ones)
        rho=local_3d_energy(
             self.ell,self.vel,self.Q,self.dQ,self.c)["rho"]
        self.assertAlmostEqual((plus-minus)/(2*h),
                   -float(rho.mean()),delta=1e-7*max(1,abs(rho.mean())))

    def test_lapse_long_projection_matches_local_density(self):
        basis=np.cos(self.x)+.2*np.sin(2*self.y)
        h=1e-5
        plus=action_3d_lapse_shift(
            self.ell,self.vel,self.Q,self.dQ,self.c,lapse=1+h*basis)
        minus=action_3d_lapse_shift(
            self.ell,self.vel,self.Q,self.dQ,self.c,lapse=1-h*basis)
        rho=local_3d_energy(
             self.ell,self.vel,self.Q,self.dQ,self.c)["rho"]
        expected=-float(np.mean(rho*basis))
        self.assertAlmostEqual((plus-minus)/(2*h),expected,
                delta=1e-7*max(1,abs(expected)))

    def test_full_shift_source_includes_traceless_commutator(self):
        m=full_trace_tensor_shift_momentum(
              self.ell,self.vel,self.Q,self.dQ,self.c)
        rms=float(np.sqrt(np.mean(m["commutator_divergence"]**2)))
        self.assertGreater(rms,1e-4)
        self.assertGreater(float(np.sqrt(np.mean(
                         m["momentum"]**2))),1e-4)
        self.assertLess(np.max(abs(m["tensor_commutator"]+
                       m["tensor_commutator"].swapaxes(0,1))),1e-12)

    def test_shift_action_variation_matches_analytic_commutator_momentum(self):
        basis=np.zeros((3,self.n,self.n,self.n))
        basis[0]=np.cos(2*self.x)*np.cos(3*self.y)
        basis[1]=.3*np.sin(self.x)
        h=1e-5
        plus=action_3d_lapse_shift(
            self.ell,self.vel,self.Q,self.dQ,self.c,shift=h*basis)
        minus=action_3d_lapse_shift(
            self.ell,self.vel,self.Q,self.dQ,self.c,shift=-h*basis)
        mom=full_trace_tensor_shift_momentum(
              self.ell,self.vel,self.Q,self.dQ,self.c)["momentum"]
        expected=-float(np.mean(np.sum(mom*basis,axis=0)))
        self.assertGreater(abs(expected),1e-5)
        self.assertAlmostEqual((plus-minus)/(2*h),expected,
                              delta=1e-7*max(1,abs(expected)))

    def test_pure_trace_momentum_reduces_to_existing_scalar_approximation(self):
        zero=np.zeros_like(self.Q)
        m=full_trace_tensor_shift_momentum(
              self.ell,self.vel,zero,zero,self.c)
        grad=np.array([np.gradient(self.ell,axis=i)
                       for i in range(3)])
        self.assertAlmostEqual(float(np.max(
                  np.abs(m["commutator_divergence"]))),0.,places=12)
        # independent spectral derivative of analytic cos/sin initial field
        ex=self.c.bf2/3*self.vel*(-.12*np.sin(2*self.x))
        ey=self.c.bf2/3*self.vel*(.18*np.cos(3*self.y))
        self.assertTrue(np.allclose(m["momentum"][0],ex,atol=1e-11))
        self.assertTrue(np.allclose(m["momentum"][1],ey,atol=1e-11))

    def test_noncommuting_tensor_modes_can_have_zero_mean_momentum(self):
        m=full_trace_tensor_shift_momentum(
             self.ell,self.vel,self.Q,self.dQ,self.c)
        self.assertLess(float(np.linalg.norm(
                          m["momentum"].mean(axis=(1,2,3)))),1e-9)
        self.assertGreater(float(np.linalg.norm(m["momentum"])),1e-5)

    def test_empty_tensor_state_has_no_ADM_sources(self):
        z=np.zeros_like(self.ell)
        Q=np.zeros_like(self.Q)
        e=local_3d_energy(z,z,Q,Q,self.c)["rho"]
        m=full_trace_tensor_shift_momentum(z,z,Q,Q,self.c)
        self.assertAlmostEqual(float(np.max(np.abs(e))),0.,places=12)
        self.assertAlmostEqual(float(np.max(np.abs(m["momentum"]))),0.,
                               places=12)

    def test_tensor_energy_scaling_quadratic_at_small_amplitude(self):
        r=local_3d_energy(self.ell,self.vel,self.Q,self.dQ,self.c)
        scaled=local_3d_energy(.1*self.ell,.1*self.vel,
                       .1*self.Q,.1*self.dQ,self.c)
        self.assertAlmostEqual(float(scaled["rho"].mean()/
                       r["rho"].mean()),.01,delta=.0003)

    def test_shift_invalid_and_lapse_nonpositive_rejected(self):
        with self.assertRaises(ValueError):
            action_3d_lapse_shift(self.ell,self.vel,
                   self.Q,self.dQ,self.c,
                   lapse=-np.ones_like(self.ell))
        with self.assertRaises(ValueError):
            action_3d_lapse_shift(self.ell,self.vel,
                   self.Q,self.dQ,self.c,
                   shift=np.ones((3,4)))

    def test_tensor_shift_status_is_restricted_parent_action(self):
        d=demo()
        self.assertGreater(d["mean_energy"],0)
        self.assertGreater(d["momentum_rms"],0)
        self.assertGreater(d["commutator_divergence_rms"],0)
        self.assertTrue(d["source_requires_full_mapper_constraint"])

if __name__=="__main__":
    unittest.main()
