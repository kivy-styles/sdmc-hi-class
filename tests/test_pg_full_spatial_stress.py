"""Independent 3D geometry, trace/traceless and off-shell Jordan stress checks."""
import unittest
import numpy as np
from sdmc.pg_full_spatial_stress import (
    SupportCoefficients,exact_riemann_squared_3d,supporting_gamma,
    action_spatial_parts,spatial_metric_projection,
    pressure_fixed_mapper,stress_projection_diagnostic,sample)
from sdmc.pg_nonlinear_conformal_stress import exact_conformal_riemann2
from sdmc.pg_quadratic_source import TraceCoefficients

class FullSpatialWellStress(unittest.TestCase):
    def setUp(self):
        self.n=24
        self.length=2*np.pi
        a=np.arange(self.n)*self.length/self.n
        self.x,self.y,self.z=np.meshgrid(a,a,a,indexing="ij")
        self.ell=.055*(np.cos(2*self.x)+.55*np.sin(3*self.y)
                          +.25*np.cos(4*self.z))
        self.vel=.02*np.sin(2*self.x)
        self.Q=np.zeros((3,3,self.n,self.n,self.n))
        self.Q[0,0]=.024*np.cos(3*self.x)
        self.Q[1,1]=-.024*np.cos(3*self.x)
        self.dQ=np.zeros_like(self.Q)
        self.c=SupportCoefficients(bf2=3.,c2=.7,m2=.2,
                 weta=.02,a_j=1.2,a_s=1.4,length=self.length)

    def test_conformal_no_traceless_Riemann_recovers_analytic_curvature(self):
        Q=np.zeros_like(self.Q)
        gamma=supporting_gamma(self.ell,Q,self.c)
        numerical=exact_riemann_squared_3d(gamma,self.length)["riemann_squared"]
        cc=TraceCoefficients(bf2=3,c_l2=.7,m_l2=.2,
              w_eta=.02,a_j=1.2,a_s=1.4,length=self.length)
        analytical=exact_conformal_riemann2(self.ell/3,cc)
        self.assertLess(np.max(np.abs(numerical-analytical))/
                        max(np.max(abs(analytical)),1e-12),2e-7)

    def test_extrinsic_traceless_spatial_state_stays_positive_metric(self):
        gamma=supporting_gamma(self.ell,self.Q,self.c)
        matrices=np.moveaxis(gamma,[0,1],[-2,-1])
        eigen=np.linalg.eigvalsh(matrices)
        self.assertGreater(float(eigen.min()),0.)
        source=exact_riemann_squared_3d(gamma,self.length)
        self.assertGreater(float(source["riemann_squared"].mean()),0.)
        self.assertGreater(float(source["riemann_squared"].min()),-1e-7)

    def test_zero_trace_nonzero_traceless_well_has_positive_energy(self):
        zero=np.zeros_like(self.ell)
        p=action_spatial_parts(zero,zero,self.Q,self.dQ,self.c)
        self.assertGreater(p["curvature"],0.)
        self.assertGreater(p["gradient"],0.)
        self.assertGreater(p["mass"],0.)
        self.assertGreater(p["density"],0.)

    def test_empty_support_has_zero_stress_and_energy(self):
        zero=np.zeros_like(self.ell)
        tensor=np.zeros_like(self.Q)
        parts=action_spatial_parts(zero,zero,tensor,tensor,self.c)
        self.assertAlmostEqual(parts["density"],0.,places=12)
        val=spatial_metric_projection(zero,zero,tensor,tensor,
                     self.c,np.eye(3))
        self.assertAlmostEqual(val,0.,places=12)

    def test_full_isotropic_spatial_stress_matches_direct_action_variation(self):
        p=action_spatial_parts(self.ell,self.vel,self.Q,self.dQ,self.c)
        numerical=spatial_metric_projection(
            self.ell,self.vel,self.Q,self.dQ,self.c,np.eye(3))
        predicted=3*pressure_fixed_mapper(p)
        self.assertAlmostEqual(numerical,predicted,
                               delta=3e-7*max(1,abs(predicted)))

    def test_curvature_isotropic_pressure_positive_one_third_energy(self):
        c=self.c;h=1e-4
        plus=action_spatial_parts(self.ell,self.vel,
                 self.Q,self.dQ,c,strain=h*np.eye(3))
        minus=action_spatial_parts(self.ell,self.vel,
                 self.Q,self.dQ,c,strain=-h*np.eye(3))
        dLcurv=-(plus["curvature"]-minus["curvature"])/(2*h)
        baseline=action_spatial_parts(self.ell,self.vel,
                 self.Q,self.dQ,c)
        self.assertAlmostEqual(dLcurv,baseline["curvature"],
                   delta=2e-7*max(1.,baseline["curvature"]))

    def test_mass_kinetic_gradient_have_expected_isotropic_pressure(self):
        p=action_spatial_parts(self.ell,self.vel,self.Q,self.dQ,self.c)
        h=1e-4
        plus=action_spatial_parts(self.ell,self.vel,
              self.Q,self.dQ,self.c,strain=h*np.eye(3))
        minus=action_spatial_parts(self.ell,self.vel,
              self.Q,self.dQ,self.c,strain=-h*np.eye(3))
        for name,sign in [("kinetic",3.),("mass",3.),("gradient",1.)]:
            d=(plus[name]-minus[name])/(2*h)
            self.assertAlmostEqual(d,sign*p[name],
                delta=2e-7*max(1,abs(d)))

    def test_mean_spatial_tensor_trace_matches_isotropic_projection(self):
        s=stress_projection_diagnostic(self.ell,self.vel,
                                       self.Q,self.dQ,self.c)
        tr=sum(s["stress_projections"][f"diag_{i}"] for i in range(3))
        self.assertAlmostEqual(tr,3*s["pressure_from_metric"],places=12)
        self.assertAlmostEqual(s["pressure_check_error"],0.,
                             delta=3e-7)
        self.assertGreater(s["anisotropic_diagonal_spread"],1e-5)

    def test_traceless_excitation_has_nonzero_homogeneous_shear_stress(self):
        z=np.zeros_like(self.ell)
        stress=stress_projection_diagnostic(z,z,self.Q,self.dQ,self.c)
        self.assertGreater(stress["anisotropic_diagonal_spread"],1e-5)

    def test_off_diagonal_metric_variation_can_detect_shear(self):
        phase=2*self.x+3*self.y
        ell=.05*np.cos(phase)
        vel=.01*np.cos(phase)
        Q=np.zeros_like(self.Q)
        Q[0,1]=.025*np.cos(phase)
        Q[1,0]=Q[0,1]
        Q[0,0]=.015*np.cos(phase)
        Q[1,1]=-Q[0,0]
        xy=np.zeros((3,3))
        xy[0,1]=xy[1,0]=.5
        p=spatial_metric_projection(ell,vel,Q,np.zeros_like(Q),
                self.c,xy)
        self.assertGreater(abs(p),1e-3)

    def test_spatially_varying_shear_metric_variation_nonzero(self):
        phase=2*self.x+3*self.y
        ell=.05*np.cos(phase)
        Q=np.zeros_like(self.Q)
        Q[0,1]=.02*np.cos(phase)
        Q[1,0]=Q[0,1]
        basis=np.zeros_like(Q)
        basis[0,1]=.5*np.cos(2*phase)
        basis[1,0]=basis[0,1]
        source=spatial_metric_projection(ell,np.zeros_like(ell),
                        Q,np.zeros_like(Q),self.c,basis)
        # A genuine inhomogeneous metric projection is supported.
        self.assertTrue(np.isfinite(source))
        self.assertGreater(abs(source),1e-5)

    def test_curvature_riemann_invariant_scaling_with_structural_ruler(self):
        p=action_spatial_parts(self.ell,self.vel,self.Q,self.dQ,self.c)
        cc=SupportCoefficients(bf2=3.,c2=.7,m2=.2,
              weta=.02,a_j=1.2,a_s=2.8,length=self.length)
        twice=action_spatial_parts(self.ell,self.vel,
                             self.Q,self.dQ,cc)
        self.assertAlmostEqual(twice["curvature"],
                        p["curvature"]/16,
                        delta=1e-9*max(1,abs(p["curvature"])))
        for name in ["kinetic","mass","gradient"]:
            self.assertAlmostEqual(twice[name],p[name],places=11)

    def test_background_scaling_cannot_fix_amplitude_from_geometry(self):
        zero=np.zeros_like(self.Q)
        p=action_spatial_parts(self.ell,self.vel,zero,zero,self.c)
        for scale in [.04,.08]:
            q=action_spatial_parts(scale*self.ell,scale*self.vel,
                                  zero,zero,self.c)
            self.assertAlmostEqual(q["density"]/p["density"],
                    scale*scale,delta=.025*scale*scale)

    def test_no_independent_CDMLike_density_claim(self):
        x=sample()
        self.assertIn("complete 4D parent action",str(x["missing"]))
        self.assertGreater(x["energy"],0.)
        self.assertGreater(abs(x["pressure_from_metric"]),1e-6)
        self.assertLess(abs(x["pressure_check_error"]),1e-6)

    def test_invalid_shapes_and_nonpositive_metric_are_rejected(self):
        with self.assertRaises(ValueError):
            supporting_gamma(self.ell,self.Q[:2],self.c)
        with self.assertRaises(ValueError):
            action_spatial_parts(self.ell,self.vel,self.Q,
               self.dQ,self.c,strain=np.ones((5,5)))
        with self.assertRaises(ValueError):
            SupportCoefficients(weta=0.)

if __name__=="__main__":
    unittest.main()
