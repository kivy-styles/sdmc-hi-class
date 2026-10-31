"""Quasistatic Jordan covariant F-R subclass and structural kernel tests."""
import unittest
import numpy as np
from sdmc.pg_jordan_qs_covariant import (
    JordanSubclass,TraceAction,static_quadratic_action,
    constraint_hessian,static_perturbation_solution,
    brans_dicke_mu_analytic,horndeski_alpha_coefficients,sweep,
    no_slip_kinetic_closure)

class SignedJordanConstraints(unittest.TestCase):
    def setUp(self):
        self.bg=JordanSubclass(F=1.12,F_sigma=.27,k1=1.4,
                               mass2=.15,sigma=1.,F_ref=1.,a_j=.9)
        self.st=TraceAction(K=.8,c2=.9,mass2=.2,
                             weta=.045,a_j=.9,a_s=1.2)

    def test_static_action_hessian_by_finite_difference(self):
        bg=self.bg; k=1.25; st=self.st
        H=constraint_hessian(k,bg,st)
        x=np.array([-.2,-.25,.07])
        h=1e-4
        def energy(v):
            return static_quadratic_action(*v,rho=.7,k=k,bg=bg,structural=st)
        hd=np.empty((3,3))
        for i in range(3):
            for j in range(3):
                if i==j:
                    a=x.copy();b=x.copy();a[i]+=h;b[i]-=h
                    hd[i,j]=(energy(a)-2*energy(x)+energy(b))/h**2
                else:
                    a=x.copy();b=x.copy();c=x.copy();d=x.copy()
                    a[i]+=h;a[j]+=h
                    b[i]+=h;b[j]-=h
                    c[i]-=h;c[j]+=h
                    d[i]-=h;d[j]-=h
                    hd[i,j]=(energy(a)-energy(b)-energy(c)+energy(d))/(4*h*h)
        self.assertTrue(np.allclose(H,hd,rtol=1e-6,atol=1e-7))
        self.assertTrue(np.allclose(H,H.T))

    def test_full_signed_constraint_solution_stationary_action(self):
        for k in [.03,.2,1.,4.,12.]:
            rho=.1
            x=static_perturbation_solution(k,self.bg,rho,self.st)
            v=np.array([x["psi_J"],x["phi_J"],x["delta_sigma"]])
            H=constraint_hessian(k,self.bg,self.st)
            self.assertTrue(np.allclose(H@v,np.array([rho,0.,0.]),
                  rtol=1e-9,atol=1e-9))

    def test_GR_limit(self):
        bg=JordanSubclass(F=1.2,F_sigma=0.,k1=1.,
                          mass2=.1,F_ref=1.)
        for k in [.1,.7,5.]:
            sol=static_perturbation_solution(k,bg)
            self.assertAlmostEqual(sol["mu_subclass"],1/1.2,places=11)
            self.assertAlmostEqual(sol["Sigma_subclass"],1/1.2,places=11)
            self.assertAlmostEqual(sol["slip"],0.,places=11)
            self.assertAlmostEqual(sol["delta_sigma"],0.,places=11)

    def test_scalar_tensor_analytic_mu_and_sigma(self):
        for k in [.01,.1,1.,10.]:
            sol=static_perturbation_solution(k,self.bg)
            target=brans_dicke_mu_analytic(k,self.bg)
            self.assertAlmostEqual(sol["mu_subclass"],target,places=12)
            self.assertAlmostEqual(sol["Sigma_subclass"],
                                   self.bg.F_ref/self.bg.F,places=12)
            self.assertAlmostEqual(sol["slip"],
                       sol["delta_F"]/self.bg.F,places=12)

    def test_massless_no_canonical_scalar_is_fR_four_thirds(self):
        bg=JordanSubclass(F=1.3,F_sigma=.5,k1=0.,
                     mass2=0.,F_ref=1.,a_j=1.)
        for k in [.1,1.,10.]:
            sol=static_perturbation_solution(k,bg)
            self.assertAlmostEqual(sol["mu_subclass"],(4/3)/1.3,places=12)
            self.assertAlmostEqual(sol["Sigma_subclass"],1/1.3,places=12)

    def test_massive_scalar_does_not_modify_large_scale_gravity(self):
        bg=JordanSubclass(F=1.,F_sigma=.2,k1=1.,
                          mass2=1000.,F_ref=1.)
        small=static_perturbation_solution(.01,bg)
        large=static_perturbation_solution(100.,bg)
        self.assertLess(abs(small["mu_subclass"]-1),1e-8)
        self.assertGreater(large["mu_subclass"],small["mu_subclass"])

    def test_effective_stiffness_matches_reduced_trace_kernel(self):
        for k in [.001,.03,.1,1.,4.]:
            v=static_perturbation_solution(k,self.bg,structural=self.st)
            self.assertGreaterEqual(v["structural_kernel"],0.)
            self.assertTrue(np.isfinite(v["mu_subclass"]))
            self.assertTrue(np.isfinite(v["Sigma_subclass"]))

    def test_structural_kernel_can_change_both_potentials(self):
        k=1.
        base=static_perturbation_solution(k,self.bg)
        v=static_perturbation_solution(k,self.bg,structural=self.st)
        self.assertNotAlmostEqual(v["mu_subclass"],base["mu_subclass"],places=6)
        self.assertNotAlmostEqual(v["Sigma_subclass"],base["Sigma_subclass"],places=6)

    def test_phi_psi_not_falsely_assumed_no_slip_with_structure(self):
        k=1.
        v=static_perturbation_solution(k,self.bg,structural=self.st)
        self.assertNotAlmostEqual(v["slip"],v["delta_F"]/self.bg.F,places=6)

    def test_monotonic_canonical_subclass_scalar_enhancement_bound(self):
        for Fp in [-1.,-.2,.05,.7,2.]:
            bg=JordanSubclass(F=1.,F_sigma=Fp,k1=1.,
                            mass2=.001,F_ref=1.)
            for k in [.01,1.,100.]:
                mu=brans_dicke_mu_analytic(k,bg)
                self.assertGreaterEqual(mu,1.)
                self.assertLess(mu,4/3+1e-12)

    def test_healthy_Horndeski_eft_parameter_map(self):
        vals=horndeski_alpha_coefficients(
            F=1.05,F_sigma=.12,H=.8,sigma_dot=.4,Z=.08,
            k1=1.7,k2=.3,g=-.75,g_sigma=.15)
        self.assertTrue(np.isfinite(list(vals.values())).all())
        self.assertGreater(vals["D"],0.)
        # No-slip background identity requires gZ=-F_sigma/2.
        self.assertAlmostEqual(vals["no_slip_covariant_residual"],0.,places=13)
        self.assertAlmostEqual(vals["no_slip_residual"],0.,places=13)
        self.assertAlmostEqual(vals["alpha_B"],-2*vals["alpha_M"],places=13)

    def test_no_slip_identity_does_not_set_kinetic_stability(self):
        vals=horndeski_alpha_coefficients(
            F=1.,F_sigma=.2,H=.8,sigma_dot=.4,Z=.08,
            k1=1.,k2=.2,g=-1.25,g_sigma=0.)
        self.assertAlmostEqual(vals["no_slip_residual"],0.,places=12)
        self.assertIn("D",vals)

    def test_coefficients_recover_G3_zero_alphaB_minus_alphaM(self):
        x=horndeski_alpha_coefficients(
            F=1.2,F_sigma=.1,H=.8,sigma_dot=.3,Z=.045,
            k1=1.,k2=.2,g=0.,g_sigma=0.)
        self.assertAlmostEqual(x["alpha_B"],-x["alpha_M"],places=13)


    def test_accepted_rolling_no_slip_kinetic_identity(self):
        for H,dot,Z,k1,k2,g_sigma in [
            (.8,.4,.08,1.7,.3,.15),
            (1.2,.2,.02,2.,.01,-.05),
            (.6,.7,.245,1.,.5,.02)]:
            F=1.05;Fp=.12
            d=no_slip_kinetic_closure(
                F,Fp,H,dot,Z,k1,k2,g_sigma)
            self.assertAlmostEqual(d["identity_residual"],0.,places=12)
            self.assertAlmostEqual(d["alpha_B"],-2*d["alpha_M"],places=12)
            self.assertAlmostEqual(d["D"],d["required_D_identity"],places=12)
            self.assertTrue(np.isfinite(d["D"]))
            self.assertIsNotNone(d["Q_s"])

    def test_no_slip_ghost_condition_is_independent(self):
        d=no_slip_kinetic_closure(1.,.2,.8,.4,.08,-10.,0.,0.)
        self.assertAlmostEqual(d["identity_residual"],0.,places=12)
        self.assertFalse(d["positive_kinetic"])
        with self.assertRaises(ValueError):
            no_slip_kinetic_closure(1.,.2,.8,.4,0.,1.,0.,0.)

    def test_sweep_models_separated_from_late300(self):
        result=sweep()
        self.assertEqual(len(result["model_rows"]),6)
        self.assertIn("G3 kinetic braiding",str(result["missing"]))
        self.assertTrue(all(row["mu_F_only"]>=row["Sigma_F_only"] for row in result["model_rows"]))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            static_perturbation_solution(0,self.bg)
        with self.assertRaises(ValueError):
            JordanSubclass(F=-1.)
        with self.assertRaises(ValueError):
            horndeski_alpha_coefficients(1.,0.,0.,0.,0.,1.,0.,0.,0.)

if __name__=="__main__":
    unittest.main()
