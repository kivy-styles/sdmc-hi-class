"""Locked-flat quadratic SDMC G3+structural-trace two-scalar action tests."""
import unittest
import numpy as np
from sdmc.pg_constraint_schur import TraceAction,effective_kernel,trace_stiffness
from sdmc.pg_rolling_g3_adm import RollingAction
from sdmc.pg_coupled_rolling_g3_trace import (
  CoupledBranch,kinetic_and_stiffness_matrices,
  quadratic_mode_lagrangian,dispersion,static_trace_elimination,
  twofield_energy,integrate_locked_flat_modes,example)

class CoupledRollingG3Trace(unittest.TestCase):
    def setUp(self):
        self.rolling=RollingAction(F=1.05,F_sigma=.12,H=.8,
                  sigma_dot=.4,k1=1.7,k2=.3,g=-.75,g_sigma=.15)
        self.trace=TraceAction(K=.6,c2=.7,mass2=.2,weta=.012,
                  a_j=1.,a_s=1.)
        self.b=CoupledBranch(self.rolling,self.trace,
                     Hdot=-.05,alpha_B_dot=-.08)

    def test_full_rolling_g3_Qs_and_structural_K_positive(self):
        v=kinetic_and_stiffness_matrices(self.b,1.)
        self.assertAlmostEqual(v["Qs"],.04902300949598315,places=12)
        self.assertAlmostEqual(v["A"][1,1],self.trace.K/2,places=13)
        self.assertEqual(v["A"][0,1],0.)
        self.assertGreater(v["Fs"],0.)
        self.assertLess(v["cs2_scalar"],1.)

    def test_quadratic_action_Hessian_against_finite_differences(self):
        k=1.35;v=kinetic_and_stiffness_matrices(self.b,k)
        h=1e-4
        def lag(x): return quadratic_mode_lagrangian(self.b,k,*x)
        base=np.array([.06,.02,-.03,.01])
        Hess=np.zeros((4,4))
        for i in range(4):
            for j in range(4):
                if i==j:
                    a=base.copy();b=base.copy();a[i]+=h;b[i]-=h
                    Hess[i,j]=(lag(a)-2*lag(base)+lag(b))/h**2
                else:
                    a=base.copy();b=base.copy();c=base.copy();d=base.copy()
                    a[i]+=h;a[j]+=h;b[i]+=h;b[j]-=h
                    c[i]-=h;c[j]+=h;d[i]-=h;d[j]-=h
                    Hess[i,j]=(lag(a)-lag(b)-lag(c)+lag(d))/(4*h*h)
        ref=np.zeros((4,4))
        ref[:2,:2]=-v["M"]
        ref[2:,2:]=v["A"]
        self.assertTrue(np.allclose(Hess,ref,rtol=1e-7,atol=1e-8))

    def test_locked_relative_curvature_is_zeta_plus_ell_third(self):
        k=1.2;v=kinetic_and_stiffness_matrices(self.b,k)
        C=v["C"]
        self.assertAlmostEqual(v["M"][0,1],C/3,places=12)
        self.assertAlmostEqual(v["M"][0,0]-v["Fs"]*k*k,C)
        self.assertAlmostEqual(v["M"][1,1]-
                   self.trace.K*v["omega0_squared"]/2,C/9)

    def test_structural_ell_mode_recovers_retarded_restoring_frequency(self):
        for k in [.02,.5,2.,10.]:
            v=kinetic_and_stiffness_matrices(self.b,k)
            eff=v["M"][1,1]/v["A"][1,1]
            self.assertAlmostEqual(eff,
                   trace_stiffness(k,self.trace)["trace_frequency_squared"],
                   delta=1e-11*max(1,eff))

    def test_joint_static_schur_matches_earlier_structural_kernel(self):
        for k in [0.,.001,.1,1.,5.,40.]:
            r=static_trace_elimination(self.b,k)
            self.assertAlmostEqual(r["schur_identity_residual"],0.,
                delta=1e-10*max(1,abs(r["predicted_stiffness"])))
            self.assertAlmostEqual(
                r["scalar_static_stiffness"],
                self.b.scalar()["Fs"]*k*k+effective_kernel(k,self.trace),
                delta=1e-10*max(1,abs(r["predicted_stiffness"])))

    def test_joint_static_trace_extremum_matches_prior_ell_solution(self):
        for k in [.01,.3,1.,3.]:
            q=static_trace_elimination(self.b,k)
            v=trace_stiffness(k,self.trace)
            expected=-3*v["Lambda"]/v["trace_frequency_squared"]
            self.assertAlmostEqual(q["ell_per_zeta"],expected,places=12)

    def test_no_curvature_gives_two_independent_modes(self):
        t=TraceAction(K=.6,c2=.7,mass2=.2,weta=0.,a_j=1.,a_s=1.)
        b=CoupledBranch(self.rolling,t,Hdot=-.05,alpha_B_dot=-.08)
        for k in [.1,1.,3.]:
            d=dispersion(b,k)
            frequencies=sorted([b.scalar()["Fs"]*k*k/b.scalar()["Qs"],
                                t.mass2+t.c2*k*k/t.a_j**2])
            self.assertTrue(np.allclose(d["eigenvalues_omega_squared"],
                                        frequencies,rtol=1e-11,atol=1e-11))

    def test_positive_energy_branch_has_positive_frozen_eigenvalues(self):
        for k in [.01,.1,1.,3.,10.]:
            d=dispersion(self.b,k)
            self.assertTrue(d["stable_kinetic"])
            self.assertTrue(d["stable_frozen_modes"])
            self.assertGreater(min(d["eigenvalues_omega_squared"]),0.)

    def test_dispersion_trace_determinant(self):
        for k in [.1,1.,4.]:
            d=dispersion(self.b,k)
            v=kinetic_and_stiffness_matrices(self.b,k)
            normalized=np.linalg.det(v["M"])/np.linalg.det(v["A"])
            self.assertAlmostEqual(
                float(np.prod(d["eigenvalues_omega_squared"])),
                normalized,delta=1e-10*max(1,abs(normalized)))

    def test_both_modes_energy_conserved_when_background_frozen(self):
        t=np.linspace(0,20,180)
        r=integrate_locked_flat_modes(
               self.b,k=1.,times=t,q0=[.2,.03],qdot0=[0.,0.])
        self.assertLess(r["fractional_energy_range"],1e-8)
        self.assertGreater(abs(r["ell"][-1]-r["ell"][0]),1e-4)

    def test_gradient_instability_if_isolated_Fs_negative_and_no_curvature(self):
        t=TraceAction(K=.6,c2=.7,mass2=.2,weta=0.,a_j=1.,a_s=1.)
        b=CoupledBranch(self.rolling,t,Hdot=-.1,alpha_B_dot=-1.)
        d=dispersion(b,1.)
        self.assertTrue(d["stable_kinetic"])
        self.assertFalse(d["stable_frozen_modes"])
        self.assertLess(min(d["eigenvalues_omega_squared"]),0.)

    def test_ghost_branch_not_hidden_by_structural_kinetic(self):
        r=RollingAction(F=1.,F_sigma=0.,H=.8,sigma_dot=.4,
           k1=-10.,k2=0.,g=0.,g_sigma=0.)
        b=CoupledBranch(r,self.trace,Hdot=-.08,alpha_B_dot=0.)
        d=dispersion(b,1.)
        self.assertFalse(d["stable_kinetic"])
        self.assertIsNone(d["eigenvalues_omega_squared"])
        with self.assertRaises(ValueError):
            integrate_locked_flat_modes(b,1.,[0.,1.],[.1,0.],[0.,0.])

    def test_structural_and_jordan_scales_cannot_be_confused(self):
        t=TraceAction(K=.6,c2=.7,mass2=.2,weta=.012,a_j=1.2,a_s=3.)
        with self.assertRaises(ValueError):
            CoupledBranch(self.rolling,t,Hdot=-.05,alpha_B_dot=-.08)

    def test_locked_flat_stability_does_not_claim_rho_g(self):
        v=example()
        self.assertLess(v["energy_fractional_range"],1e-8)
        self.assertEqual(len(v["dispersion"]),6)
        self.assertIn("physical rho_g",str(v["missing"]))
        self.assertTrue(all(d["stable"] for d in v["dispersion"]))

if __name__=="__main__":unittest.main()
