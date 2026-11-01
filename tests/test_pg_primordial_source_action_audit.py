"""Accepted SDMC source-channel variation: exact mixed-Hessian and real signed-Jordan checkpoint."""
import unittest
import numpy as np
from sdmc.pg_primordial_source_action_audit import (
  SourceAudit,curvature_stiffness,lambda_trace,
  primordial_support_energy,bare_jordan_matter_action,
  analytic_action_gradient,analytic_hessian,
  finite_difference_hessian,stationary_solution,
  metric_source_grad_ell,curvature_loads_scalar,
  scalar_kinetic_coefficient_mixing,
  independent_well_energy_at_W_zero,
  source_free_static_lump_identity,
  candidate_explicit_direct_source,candidate_source_solution_at_W_zero,
  on_shell_charged_source_audit,report)
from sdmc.pg_jordan_qs_covariant import (
    JordanSubclass,static_perturbation_solution)
from sdmc.pg_constraint_schur import TraceAction

class PrimordialSourceActionAudit(unittest.TestCase):
    def setUp(self):
        self.e=SourceAudit()
        self.off=SourceAudit(w_eta=0.)

    def test_minimally_coupled_Jordan_material_action_has_no_direct_ell_source(self):
        e=self.e
        for w in [0.,.001,.2]:
            src=SourceAudit(w_eta=w)
            self.assertEqual(metric_source_grad_ell(src),0.)
            self.assertEqual(bare_jordan_matter_action(src,.3,-.7),
                             bare_jordan_matter_action(src,.3,+.5))

    def test_curvature_selector_off_blocks_linear_metric_trace_mixing(self):
        H=analytic_hessian(self.off)
        self.assertEqual(H[0,1],0.)
        self.assertEqual(H[1,0],0.)
        self.assertAlmostEqual(H[0,0],self.off.Qmetric)
        self.assertAlmostEqual(H[1,1],
                            self.off.K*self.off.omega0_sq/2)

    def test_curvature_selector_nonzero_generates_exact_W_mixing(self):
        e=self.e;H=analytic_hessian(e)
        C=curvature_stiffness(e)
        self.assertGreater(C,0.)
        self.assertAlmostEqual(H[0,1],C/3)
        self.assertAlmostEqual(H[1,0],C/3)
        self.assertAlmostEqual(H[0,0]-e.Qmetric,C)

    def test_metric_trace_coupling_scales_as_W(self):
        a=SourceAudit(w_eta=.02)
        b=SourceAudit(w_eta=.06)
        self.assertAlmostEqual(analytic_hessian(b)[0,1],
               3*analytic_hessian(a)[0,1],places=13)

    def test_curvature_mixing_is_fourth_order_in_k(self):
        a=SourceAudit(k=.5)
        b=SourceAudit(k=1.)
        self.assertAlmostEqual(analytic_hessian(b)[0,1],
                    16*analytic_hessian(a)[0,1],places=13)
        self.assertEqual(analytic_hessian(SourceAudit(k=0.))[0,1],0.)

    def test_full_quadratic_Hessian_finite_differences(self):
        for e in [self.e,self.off,SourceAudit(k=.05),
                  SourceAudit(K=.3,w_eta=.02)]:
            n=finite_difference_hessian(e)
            exact=analytic_hessian(e)
            self.assertTrue(np.allclose(n,exact,rtol=1e-7,atol=1e-8))

    def test_first_action_gradients_against_independent_finite_difference(self):
        e=self.e;h=1e-5
        p=np.array([.23,-.1])
        g=analytic_action_gradient(e,*p)
        for j in [0,1]:
            pos=p.copy();neg=p.copy();pos[j]+=h;neg[j]-=h
            num=(primordial_support_energy(e,*pos)-
                 primordial_support_energy(e,*neg))/(2*h)
            self.assertAlmostEqual(num,g[j],delta=2e-10)

    def test_stationary_two_field_vs_static_Schur(self):
        for e in [self.e,self.off,SourceAudit(w_eta=.002),
                  SourceAudit(k=.1),SourceAudit(K=.07)]:
            r=stationary_solution(e)
            self.assertAlmostEqual(r["u"],r["schur_u"],places=12)
            self.assertAlmostEqual(r["ell"],r["schur_ell"],places=12)
            self.assertGreater(r["hessian_det"],0.)
            self.assertFalse(r["physical_mu_or_lensing_predicted"])

    def test_W_zero_static_ell_vanishes_despite_nonzero_metric_source(self):
        r=stationary_solution(self.off)
        self.assertEqual(r["ell"],0.)
        self.assertAlmostEqual(r["u"],self.off.j_metric/
                              self.off.Qmetric)
        self.assertEqual(r["extra_metric_stiffness"],0.)
        self.assertNotEqual(r["u"],0.)

    def test_curvature_kernel_suppresses_positive_static_metric_stiffness_toy(self):
        r=stationary_solution(self.e)
        self.assertGreater(r["extra_metric_stiffness"],0.)
        self.assertLess(r["u"],r["baseline_u"])
        # These are not the signed real Jordan constraints.
        self.assertFalse(r["physical_mu_or_lensing_predicted"])

    def test_tiny_W_limit_response_tends_to_baseline(self):
        e=SourceAudit(w_eta=1e-8)
        r=stationary_solution(e)
        self.assertLess(abs(r["u"]-r["baseline_u"]),1e-7)
        self.assertLess(abs(r["ell"]),1e-7)

    def test_nonzero_F_sigma_curvature_loads_scalar_even_when_W_zero(self):
        e=self.off
        d=curvature_loads_scalar(1.1,e.F_sigma,.8)
        self.assertTrue(d["nonzero_if_curvature_and_Fsigma"])
        self.assertFalse(d["depends_on_structural_L"])
        self.assertAlmostEqual(d["sigma_equation_curvature_term"],.08)

    def test_scalar_curvature_variation_independent_finite_difference(self):
        F=1.1;Fp=.18;R=.17;h=1e-5
        p=curvature_loads_scalar(F,Fp,R,h)["action_density"]
        m=curvature_loads_scalar(F,Fp,R,-h)["action_density"]
        num=(p-m)/(2*h)
        self.assertAlmostEqual(num,Fp*R/2,places=11)

    def test_excited_L_can_mix_with_sigma_via_running_K_without_direct_source(self):
        e=self.e
        self.assertEqual(scalar_kinetic_coefficient_mixing(e,0.),0.)
        self.assertAlmostEqual(
            scalar_kinetic_coefficient_mixing(e,.12),
            .5*e.K_sigma*e.omega0_sq*.12)
        # finite difference of coefficient-dependent energy
        ell=.15;h=1e-5
        def E(sig,l):
            return .25*(e.K+e.K_sigma*sig)*e.omega0_sq*l*l
        num=((E(h,ell+h)-E(h,ell-h))-
             (E(-h,ell+h)-E(-h,ell-h)))/(4*h*h)
        self.assertAlmostEqual(num,
             scalar_kinetic_coefficient_mixing(e,ell),places=10)

    def test_preexisting_free_mode_can_gravitate_quadratically(self):
        K=1;mass=2
        self.assertEqual(independent_well_energy_at_W_zero(K,mass,0,0),0.)
        a=independent_well_energy_at_W_zero(K,mass,.1,.05)
        b=independent_well_energy_at_W_zero(K,mass,.2,.1)
        self.assertGreater(a,0.)
        self.assertAlmostEqual(b,4*a)
        self.assertGreater(a,
             independent_well_energy_at_W_zero(K,mass,.1,0))

    def test_nonzero_static_ell_fails_positive_unforced_spatial_equation(self):
        n=24;L=2*np.pi
        a=np.arange(n)*L/n
        x,y,z=np.meshgrid(a,a,a,indexing="ij")
        ell=.04*np.cos(3*x)+.06*np.sin(2*y)+.02*np.cos(4*z)
        r=source_free_static_lump_identity(
              ell,K=1.,mass2=.7,c2=.5,a_j=1.2,length=L)
        self.assertGreater(r["integrated_field_times_EL"],0.)
        self.assertAlmostEqual(r["identity_residual"],0.,delta=1e-13)
        self.assertFalse(r["supports_nonzero_sourcefree_static_lump"])
        self.assertFalse(r["Jordan_gravitational_well_forbidden"])

    def test_static_lump_identity_zero_only_for_zero_deformation(self):
        n=20;f=np.zeros((n,n,n))
        r=source_free_static_lump_identity(
             f,K=1.,mass2=.7,c2=.5,a_j=1.,length=2*np.pi)
        self.assertEqual(r["integrated_field_times_EL"],0.)
        self.assertEqual(r["positive_mass_gradient_norm"],0.)

    def test_mass_only_nonzero_constant_field_has_positive_EL(self):
        n=20;f=np.ones((n,n,n))*.12
        r=source_free_static_lump_identity(
            f,K=1.,mass2=1.,c2=0.,a_j=1.,length=2*np.pi)
        self.assertAlmostEqual(r["integrated_field_times_EL"],.12**2)
        self.assertAlmostEqual(r["identity_residual"],0.,delta=1e-12)

    def test_new_direct_source_operator_NOT_in_accepted_action(self):
        e=self.off
        gnew=.5
        ell=candidate_source_solution_at_W_zero(e,gnew)
        self.assertNotEqual(ell,0.)
        h=1e-5
        u=.1
        p=candidate_explicit_direct_source(e,u,ell+h,gnew)
        m=candidate_explicit_direct_source(e,u,ell-h,gnew)
        self.assertAlmostEqual((p-m)/(2*h),0.,delta=1e-10)
        self.assertAlmostEqual(
            candidate_source_solution_at_W_zero(e,0.),0.)

    def test_no_W_no_k_mode_requires_no_new_material_source(self):
        for e in [SourceAudit(k=0.),self.off]:
            self.assertEqual(curvature_stiffness(e),0.)
            self.assertEqual(lambda_trace(e),0.)
            self.assertEqual(stationary_solution(e)["ell"],0.)

    def test_signed_Jordan_constraints_return_scalar_only_when_W_off(self):
        bg=JordanSubclass(F=1.,F_sigma=.2,k1=1.,
                          mass2=.1,sigma=1.,a_j=1.)
        t=TraceAction(K=1.,c2=.7,mass2=.2,
                         weta=0.,a_j=1.,a_s=1.)
        for k in [.1,.5,1.,3.]:
            a=static_perturbation_solution(k,bg,rho=1.)
            b=static_perturbation_solution(k,bg,rho=1.,structural=t)
            for field in ["mu_subclass","Sigma_subclass",
                          "psi_J","phi_J","delta_sigma"]:
                self.assertAlmostEqual(a[field],b[field],places=12)
            self.assertEqual(b["structural_kernel"],0.)

    def test_signed_Jordan_W_on_response_not_identical_to_scalar_only(self):
        bg=JordanSubclass()
        t=TraceAction(K=1.,c2=.8,mass2=.1,
                         weta=.03,a_j=1.,a_s=1.)
        a=static_perturbation_solution(2,bg,rho=1.)
        b=static_perturbation_solution(2,bg,rho=1.,structural=t)
        self.assertNotAlmostEqual(a["mu_subclass"],b["mu_subclass"],places=5)
        self.assertGreater(b["structural_kernel"],0.)

    def test_baryonic_source_zero_and_W_off_gives_vacuum_in_linear_response(self):
        e=SourceAudit(w_eta=0.,j_metric=0.)
        r=stationary_solution(e)
        self.assertEqual(r["ell"],0.)
        self.assertEqual(r["u"],0.)
        self.assertEqual(r["baseline_u"],0.)

    def test_status_uses_sources_not_assumed_CDM_abundance(self):
        d=report()
        self.assertEqual(d["cases"]["W_exactly_zero"]["linear_metric_trace_cross_Hessian"],0.)
        self.assertEqual(d["cases"]["W_exactly_zero"]["minimal_material_direct_ell_source"],0.)
        self.assertGreater(
            d["cases"]["W_exactly_zero"]["L_energy_at_nonzero_initial_state"],0.)
        self.assertEqual(
            d["conditional_new_operator_example"]["extra_operator_not_in_accepted_action"],
            "-g_new*j_metric*ell")
        self.assertIn("physical complete 4D mapper",str(d["remaining"]))

    def test_input_validation(self):
        for args in [{"K":0},{"omega0_sq":0},{"Qmetric":0},
                     {"w_eta":-1},{"k":-1}]:
            with self.assertRaises(ValueError):
                SourceAudit(**args)
        with self.assertRaises(ValueError):
            candidate_source_solution_at_W_zero(self.e,1.)

if __name__=="__main__":
    unittest.main()
