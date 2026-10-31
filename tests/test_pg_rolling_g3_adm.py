"""Independent finite-difference and EFT benchmarks for rolling G3 ADM."""
import unittest
import numpy as np
from sdmc.pg_rolling_g3_adm import (
    RollingAction,coefficients_at_lapse,lapse_derivatives,alphas,
    adm_kinetic_hessian,kinetic_from_constraints,exact_adm_density_cosine,
    numeric_kinetic_hessian,isolated_scalar_gradient,
    no_slip_gradient_identity,diagnostics)
from sdmc.pg_jordan_qs_covariant import (
    horndeski_alpha_coefficients,no_slip_kinetic_closure)

class RollingADMChecks(unittest.TestCase):
    def setUp(self):
        self.a=RollingAction(F=1.05,F_sigma=.12,H=.8,
             sigma_dot=.4,k1=1.7,k2=.3,g=-.75,g_sigma=.15)
        self.b=RollingAction(F=1.12,F_sigma=-.08,H=.95,
             sigma_dot=.3,k1=2.1,k2=.2,g=.15,g_sigma=-.1,V=.23,a=.91)

    def test_adm_coefficients_and_derivatives_match_finite_difference(self):
        for a in [self.a,self.b]:
            d=lapse_derivatives(a)
            h=1e-4
            a0,b0=coefficients_at_lapse(a,1.)
            plus=coefficients_at_lapse(a,1+h)
            minus=coefficients_at_lapse(a,1-h)
            self.assertAlmostEqual(a0,d["A20"],places=12)
            self.assertAlmostEqual(b0,d["A30"],places=12)
            for key,i in [("A21",0),("A31",1)]:
                num=(plus[i]-minus[i])/(2*h)
                self.assertAlmostEqual(num,d[key],delta=1e-7*max(1,abs(num)))
            for key,i in [("A22",0),("A32",1)]:
                num=(plus[i]-2*(a0,b0)[i]+minus[i])/h**2
                self.assertAlmostEqual(num,d[key],delta=1e-7*max(1,abs(num)))

    def test_full_exact_adm_action_numeric_hessian(self):
        for a in [self.a,self.b]:
            analytic=adm_kinetic_hessian(a)["Hessian"]
            numeric=numeric_kinetic_hessian(a,k=3,n=512,eps=3e-4)
            self.assertTrue(np.allclose(analytic,numeric,rtol=2e-5,atol=2e-6),
                       msg=f"analytic {analytic}, numeric {numeric}")

    def test_hessian_has_shift_as_constraint_not_kinetic_field(self):
        H=adm_kinetic_hessian(self.a)["Hessian"]
        self.assertEqual(H[1,1],0.)
        self.assertAlmostEqual(H[1,2],-self.a.F)
        self.assertAlmostEqual(H[0,2],3*H[0,1],places=12)

    def test_lapse_shift_constraints_match_analytic_Qs(self):
        for a in [self.a,self.b]:
            r=kinetic_from_constraints(a)
            H=adm_kinetic_hessian(a)["Hessian"]
            vel=np.array([r["alpha_per_zdot"],r["B_per_zdot"],1.])
            grad=H@vel
            self.assertAlmostEqual(grad[0],0.,places=12)
            self.assertAlmostEqual(grad[1],0.,places=12)
            self.assertAlmostEqual(float(vel@H@vel),r["Qs_ADM"],places=11)
            self.assertAlmostEqual(r["identity_residual"],0.,
                delta=2e-11*max(1,abs(r["Qs_ADM"])))

    def test_EFT_alpha_coefficients_match_previously_accepted_action(self):
        a=self.a;r=alphas(a)
        earlier=horndeski_alpha_coefficients(
          F=a.F,F_sigma=a.F_sigma,H=a.H,sigma_dot=a.sigma_dot,
          Z=a.Z0,k1=a.k1,k2=a.k2,g=a.g,g_sigma=a.g_sigma)
        for key in ["alpha_M","alpha_B","alpha_K","D"]:
            self.assertAlmostEqual(r[key],earlier[key],places=12)

    def test_alphaB_braiding_from_lapse_derivative(self):
        for a in [self.a,self.b]:
            d=lapse_derivatives(a)
            self.assertAlmostEqual(
                alphas(a)["alpha_B"],-d["A31"]/(a.F*a.H),
                places=13)

    def test_canonical_GR_scalar_Qs(self):
        a=RollingAction(F=1.,F_sigma=0.,H=.8,sigma_dot=.4,
              k1=2.,k2=0.,g=0.,g_sigma=0.)
        r=kinetic_from_constraints(a)
        self.assertAlmostEqual(r["Qs_ADM"],
              a.k1*a.Z0/a.H**2,places=12)
        self.assertAlmostEqual(r["alpha_per_zdot"],1/a.H,places=12)

    def test_positive_kessence_Qs(self):
        a=RollingAction(F=1.,F_sigma=0.,H=.8,sigma_dot=.4,
             k1=1.,k2=.5,g=0.,g_sigma=0.)
        expected=a.Z0*(a.k1+6*a.k2*a.Z0)/a.H**2
        self.assertAlmostEqual(
            kinetic_from_constraints(a)["Qs_ADM"],expected,places=12)

    def test_no_slip_K_kinetic_constraint_matches_covariant_closure(self):
        a=self.a
        self.assertAlmostEqual(alphas(a)["no_slip_residual"],0.,places=12)
        k=no_slip_kinetic_closure(
             F=a.F,F_sigma=a.F_sigma,H=a.H,
             sigma_dot=a.sigma_dot,Z=a.Z0,k1=a.k1,
             k2=a.k2,g_sigma=a.g_sigma)
        r=kinetic_from_constraints(a)
        self.assertAlmostEqual(r["Qs_ADM"],k["Q_s"],places=12)
        self.assertTrue(r["scalar_no_ghost"])

    def test_ghost_diagnostic_rejects_negative_Qs(self):
        a=RollingAction(F=1.,F_sigma=0.,H=.8,sigma_dot=.4,
               k1=-10.,k2=0.,g=0.,g_sigma=0.)
        r=kinetic_from_constraints(a)
        self.assertLess(r["Qs_ADM"],0.)
        self.assertFalse(r["scalar_no_ghost"])

    def test_degenerate_braiding_constraint_is_detected(self):
        a=RollingAction(F=1.,F_sigma=0.,H=1.,sigma_dot=1.,
              k1=1.,k2=0.,g=2.,g_sigma=0.)
        self.assertAlmostEqual(alphas(a)["alpha_B"],2.)
        with self.assertRaises(ValueError):
            kinetic_from_constraints(a)

    def test_isolated_canonical_scalar_sound_speed(self):
        a=RollingAction(F=1.,F_sigma=0.,H=.8,sigma_dot=.4,
            k1=1.,k2=0.,g=0.,g_sigma=0.)
        Hdot=-a.Z0/a.F
        c=isolated_scalar_gradient(a,Hdot)
        self.assertAlmostEqual(c["Qs"],c["Fs"],places=12)
        self.assertAlmostEqual(c["cs2"],1.,places=12)

    def test_no_slip_gradient_two_derivations_agree(self):
        a=self.a
        for Hdot,mdot in [(-.08,.0),(-.1,.02),(-.3,-.04)]:
            r=no_slip_gradient_identity(a,Hdot,mdot)
            self.assertAlmostEqual(r["residual"],0.,places=11)
            self.assertAlmostEqual(r["Fs_NoSlip"],r["Fs_direct"],places=11)

    def test_no_slip_gradient_stability_not_automatic(self):
        a=self.a
        p=no_slip_gradient_identity(a,Hdot=-.1,alpha_M_dot=.01)
        n=no_slip_gradient_identity(a,Hdot=-.1,alpha_M_dot=.5)
        self.assertGreater(p["Fs_NoSlip"],0.)
        self.assertLess(n["Fs_NoSlip"],0.)
        self.assertTrue(p["kinetic_positive"])
        self.assertTrue(n["kinetic_positive"])

    def test_no_slip_gradient_rejects_off_branch(self):
        with self.assertRaises(ValueError):
            no_slip_gradient_identity(self.b,Hdot=-.1,alpha_M_dot=0.)

    def test_background_no_tadpole_is_not_automatically_imposed(self):
        r=diagnostics()
        self.assertLess(r["max_Hessian_difference"],2e-5)
        self.assertIn("full supporting-space ADM Hessian",str(r["missing"]))

if __name__=="__main__":unittest.main()
