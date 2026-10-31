"""Three-dimensional scalar structural geometry and coarse-grained source tests."""
import unittest
import numpy as np
from sdmc.pg_quadratic_source import (
    TraceCoefficients,quadratic_energy_density)
from sdmc.pg_quadratic_source_3d import (
    geometry_derivatives,curvature_energy_3d,quadratic_energy_3d,
    mean_quadratic_source_3d,gaussian_coarse_grain_3d,
    long_mode_contrast_3d,gradient_anisotropy_proxy_3d,three_d_demo)

class Geometric3D(unittest.TestCase):
    def setUp(self):
        self.n=32;self.L=2*np.pi
        ax=np.arange(self.n)*self.L/self.n
        self.x,self.y,self.z=np.meshgrid(ax,ax,ax,indexing="ij")
        self.c=TraceCoefficients(bf2=2,c_l2=.7,m_l2=.3,
                               w_eta=.05,a_j=1.3,a_s=1.1,length=self.L)

    def test_exact_3d_conformal_curvature_single_mode(self):
        kx,ky=3,4
        A=.6
        q=A*np.cos(kx*self.x+ky*self.y)
        obs=float(curvature_energy_3d(q,self.c).mean())
        wanted=2*self.c.w_eta/self.c.a_s**4*A**2*(kx*kx+ky*ky)**2
        self.assertAlmostEqual(obs,wanted,delta=abs(wanted)*1e-10)

    def test_parseval_hessian_mean_equal_laplacian_mean(self):
        q=np.cos(3*self.x+4*self.y)+.7*np.sin(5*self.z+2*self.y)
        _,H,lap=geometry_derivatives(q,self.L)
        hess_mean=sum(np.mean(H[i][j]**2)
                      for i in range(3) for j in range(3))
        self.assertAlmostEqual(hess_mean,float(np.mean(lap**2)),
                               delta=1e-9*float(np.mean(lap**2)))

    def test_3d_source_reduces_to_1d_embedded_profile(self):
        a=np.arange(self.n)*self.L/self.n
        u=.4*np.cos(5*a)
        field=np.broadcast_to(u[:,None,None],(self.n,)*3).copy()
        mean_3=mean_quadratic_source_3d(field,self.c)
        mean_1=float(quadratic_energy_density(u,self.c).mean())
        self.assertAlmostEqual(mean_3,mean_1,delta=1e-9*abs(mean_1))

    def test_isotropic_gradient_correlation_for_three_axes(self):
        k=6
        ell=(np.cos(k*self.x)+np.cos(k*self.y)+
             np.cos(k*self.z))/np.sqrt(3)
        T=gradient_anisotropy_proxy_3d(ell,self.L)
        self.assertLess(float(np.linalg.norm(T)),1e-10)

    def test_single_axis_has_nonzero_gradient_anisotropy(self):
        ell=np.cos(5*self.x)
        self.assertGreater(float(np.linalg.norm(
            gradient_anisotropy_proxy_3d(ell,self.L))),1)

    def test_long_density_mode_from_3d_short_geometry(self):
        k=6; eps=.08;width=.3
        short=(np.cos(k*self.x)+np.cos(k*self.y)+
               np.cos(k*self.z))/np.sqrt(3)
        ell=(1+eps*np.cos(self.x))*short
        obs,contrast=long_mode_contrast_3d(ell,self.c,k_long=1,width=width)
        expected=2*eps*np.exp(-.5*width**2)
        self.assertLess(abs(obs-expected),.01)
        self.assertAlmostEqual(float(contrast.mean()),0.,places=11)

    def test_three_d_demo_is_finite_and_consistent(self):
        result=three_d_demo()
        self.assertAlmostEqual(result["short_field_mean"],0.,places=12)
        self.assertLess(abs(result["long_density_contrast_cosine_amplitude"]-
                     result["separated_scale_approximation"]),.01)
        self.assertLess(result["unmodulated_gradient_anisotropy_norm"],1e-9)
        self.assertGreater(result["mean_action_energy_arbitrary_units"],0.)

if __name__=="__main__":
    unittest.main()
