"""Retarded structural trace, time-memory and Jordan rank-one response tests.

These verify reduced ADM structural dynamics and exact linear algebra.
They do not establish the complete accepted late300 rolling-G3 Hessian.
"""
import unittest
import numpy as np
from sdmc.pg_constraint_schur import TraceAction, effective_kernel
from sdmc.pg_jordan_qs_covariant import JordanSubclass, constraint_hessian
from sdmc.pg_retarded_jordan_response import (
    ResponseEpoch,trace_response,static_adequacy,step_trace_solution,
    analytic_step_trace,solve_rank_one,structural_jordan_rank_one,
    example,local_trace_stability,projected_metric_response)

class RetardedStructuralResponse(unittest.TestCase):
    def setUp(self):
        self.act=TraceAction(K=1.3,c2=.75,mass2=.5,weta=.06,
                             a_j=.94,a_s=1.1)
        self.epoch=ResponseEpoch(self.act,H_j=.1,Kdot_over_K=.035)

    def test_exact_zero_frequency_kernel_recovers_old_static_constraint(self):
        for k in [0.,.01,.2,1.,4.,15.]:
            d=trace_response(k,0.,self.epoch)
            self.assertAlmostEqual(d["K_response"].real,
                                   effective_kernel(k,self.act),
                                   delta=1e-11*max(1,abs(d["K_response"].real)))
            self.assertAlmostEqual(d["K_response"].imag,0.)
            self.assertAlmostEqual(
                d["q_over_u"].real,
                d["coefficients"]["omega0_squared"]/
                    d["coefficients"]["omegaeff_squared"],
                places=13)

    def test_trace_equation_exact_harmonic_response(self):
        k=1.4;freq=.8
        d=trace_response(k,freq,self.epoch)
        c=d["coefficients"]
        l=d["ell_over_u"]
        u=1+0j
        lhs=(c["omegaeff_squared"]-freq**2-1j*c["Gamma"]*freq)*l
        rhs=-3*c["Lambda"]*u
        self.assertAlmostEqual(abs(lhs-rhs),0.,places=12)
        self.assertAlmostEqual(
            abs(d["q_over_u"]-(1+l/3)),0.,places=12)

    def test_quasistatic_error_small_when_frequencies_small(self):
        k=1.7
        w0=np.sqrt(self.epoch.coefficients(k)["omegaeff_squared"])
        q=static_adequacy(k,w0*1e-4,self.epoch)
        self.assertLess(q["scaled_kernel_difference"],1e-4)
        self.assertLess(q["frequency_ratio"],.001)

    def test_near_resonance_static_kernel_unreliable(self):
        k=1.2
        epoch=ResponseEpoch(self.act,H_j=.003,Kdot_over_K=0.)
        omega=np.sqrt(epoch.coefficients(k)["omegaeff_squared"])
        q=static_adequacy(k,omega,epoch)
        self.assertGreater(q["scaled_kernel_difference"],.5)
        self.assertAlmostEqual(q["frequency_ratio"],1.)

    def test_retarded_damping_sign(self):
        for freq in [.1,.5,1.,5.]:
            r=trace_response(2.,freq,self.epoch)
            self.assertLess(r["K_response"].imag,0.)

    def test_high_frequency_structural_trace_cannot_follow_metric(self):
        r=trace_response(1.3,1e6,self.epoch)
        self.assertLess(abs(r["ell_over_u"]),1e-9)
        self.assertLess(abs(r["q_over_u"]-1),1e-9)
        self.assertAlmostEqual(r["K_response"].real,r["coefficients"]["C"],
                               delta=1e-9*max(1,r["coefficients"]["C"]))

    def test_frictionless_pole_is_not_inverted(self):
        epoch=ResponseEpoch(self.act,H_j=0.)
        w2=epoch.coefficients(1.)["omegaeff_squared"]
        with self.assertRaises(ValueError):
            trace_response(1.,np.sqrt(w2),epoch)

    def test_causal_step_matches_closed_form_time_domain(self):
        ep=ResponseEpoch(self.act,H_j=.3)
        k=1.3;t=np.linspace(0,20,130)
        numeric=step_trace_solution(t,k,ep,.4)
        analytic=analytic_step_trace(t,k,ep,.4)
        self.assertLess(float(np.max(np.abs(numeric["ell"]-analytic))),2e-8)
        self.assertAlmostEqual(numeric["ell"][0],0.,places=14)
        self.assertAlmostEqual(numeric["q"][0],.4,places=14)
        self.assertAlmostEqual(numeric["ell"][-1],numeric["ell_static"],
                 delta=1e-4*max(1,abs(numeric["ell_static"])))

    def test_static_adequacy_requires_more_than_hubble_ratio(self):
        ep=ResponseEpoch(self.act,H_j=.02)
        c=ep.coefficients(1.3)
        # Even with Omega_eff / H >= 30, a driver near resonance
        # violates the static elimination assumption.
        w=np.sqrt(c["omegaeff_squared"])
        self.assertGreater(w/ep.H_j,30.)
        self.assertGreater(static_adequacy(1.3,w,ep)[
            "scaled_kernel_difference"],.5)


    def test_positive_damping_gives_stable_trace_roots(self):
        gate=local_trace_stability(1.3,self.epoch)
        self.assertTrue(gate["no_exponential_trace_growth"])
        self.assertLess(gate["max_local_growth_rate"],0.)
        self.assertEqual(len(gate["roots"]),2)

    def test_negative_running_could_destabilize_trace_despite_positive_energy(self):
        H=.05
        ep=ResponseEpoch(self.act,H_j=H,Kdot_over_K=-4*H)
        gate=local_trace_stability(1.3,ep)
        self.assertFalse(gate["no_exponential_trace_growth"])
        self.assertAlmostEqual(gate["max_local_growth_rate"],H/2,places=12)
        self.assertAlmostEqual(gate["Jordan_kinetic_running"],-4.)

    def test_marginal_damping_boundary(self):
        H=.05
        ep=ResponseEpoch(self.act,H_j=H,Kdot_over_K=-3*H)
        gate=local_trace_stability(1.3,ep)
        self.assertTrue(gate["no_exponential_trace_growth"])
        self.assertAlmostEqual(gate["max_local_growth_rate"],0.,places=12)

    def test_k_zero_curvature_decouples(self):
        q=trace_response(0.,.1,self.epoch)
        self.assertEqual(q["K_response"],0.)
        self.assertEqual(q["ell_over_u"],0.)

    def test_rank_one_random_symmetric_indefinite_operators(self):
        rng=np.random.default_rng(2026)
        for dim in [3,5,8]:
            m=rng.normal(size=(dim,dim))
            H=m+m.T+np.diag(np.linspace(-3.,3.,dim))
            H=H+1j*.04*np.eye(dim)
            v=rng.normal(size=dim)
            b=rng.normal(size=dim)
            lam=.21-.15j
            r=solve_rank_one(H,v,b,lam)
            direct=np.linalg.solve(r["updated_operator"],b)
            self.assertTrue(np.allclose(r["response"],direct,
                                   rtol=1e-10,atol=1e-10))

    def test_physical_signed_subclass_static_reduces_to_previous_Hessian(self):
        bg=JordanSubclass(F=1.1,F_sigma=.2,k1=1.1,mass2=.25,
                            sigma=1.2,F_ref=1.,a_j=self.act.a_j)
        k=.9
        H0=constraint_hessian(k,bg,None)
        Hfull=constraint_hessian(k,bg,self.act)
        v=np.array([0.,-1.,1/bg.sigma])
        r=structural_jordan_rank_one(
            k,0.,self.epoch,H0,np.array([1.,0.,0.]),v)
        self.assertTrue(np.allclose(r["updated_operator"],Hfull,
                                   rtol=1e-12,atol=1e-12))
        direct=np.linalg.solve(Hfull,np.array([1.,0.,0.]))
        self.assertTrue(np.allclose(r["response"],direct,
                                   rtol=1e-12,atol=1e-12))

    def test_signed_subclass_finite_frequency_is_algebraically_exact(self):
        bg=JordanSubclass(F=1.1,F_sigma=.2,k1=1.1,mass2=.25,
                            sigma=1.2,F_ref=1.,a_j=self.act.a_j)
        k=1.1;om=.6
        H0=constraint_hessian(k,bg,None)
        v=np.array([0.,-1.,1/bg.sigma])
        r=structural_jordan_rank_one(k,om,self.epoch,
                                     H0,np.array([1.,0.,0.]),v)
        direct=np.linalg.solve(r["updated_operator"],[1.,0.,0.])
        self.assertTrue(np.allclose(r["response"],direct,
                                   rtol=1e-11,atol=1e-11))
        self.assertNotAlmostEqual(
            float(r["structural_response"]["K_response"].imag),0.,places=8)

    def test_modified_constraint_near_pole_is_flagged(self):
        H=np.eye(3)
        v=np.array([1.,0.,0.])
        b=np.array([1.,0.,0.])
        with self.assertRaises(ValueError):
            solve_rank_one(H,v,b,-1.)


    def test_geometric_source_can_induce_slip_from_zero_slip_base(self):
        H0=np.eye(3)
        source=np.array([1.,1.,0.])  # toy baseline Psi=Phi
        v=np.array([0.,-1.,1.])
        out=projected_metric_response(H0,v,source,lam=-.2)
        self.assertAlmostEqual(out["baseline_slip"],0.,places=13)
        self.assertGreater(abs(out["slip"]),.1)
        self.assertAlmostEqual(
            abs(out["induced_slip"]-out["induced_slip_formula"]),0.,places=12)
        self.assertNotAlmostEqual(out["weyl"].real,
                                  out["baseline_weyl"].real,places=9)

    def test_geometric_projection_can_preserve_no_slip_by_alignment(self):
        H0=np.eye(3)
        source=np.array([1.,1.,0.])
        v=np.array([1.,1.,-1.])
        out=projected_metric_response(H0,v,source,lam=.1)
        self.assertAlmostEqual(out["baseline_slip"],0.,places=13)
        self.assertAlmostEqual(out["slip"],0.,places=13)
        self.assertAlmostEqual(out["induced_slip_formula"],0.,places=13)

    def test_dynamic_jerk_splits_metric_potentials_only_conditionally(self):
        H0=np.array([[2.,.4,0.],[.4,1.8,.3],[0.,.3,1.2]])
        b=np.array([.8,1.,0.])
        v=np.array([0.,-1.,1.])
        ep=self.epoch;k=1.2
        dyn=trace_response(k,.6,ep)["K_response"]
        out=projected_metric_response(H0,v,b,lam=-2*dyn)
        full=np.linalg.solve(H0-2*dyn*np.outer(v,v),b)
        self.assertAlmostEqual(abs(out["psi"]-full[0]),0.,places=11)
        self.assertAlmostEqual(abs(out["phi"]-full[1]),0.,places=11)
        self.assertAlmostEqual(abs(out["induced_slip"]-
                                   out["induced_slip_formula"]),0.,places=11)

    def test_invalid_times_and_frequency(self):
        with self.assertRaises(ValueError):
            step_trace_solution([0,1,1],1,self.epoch,.3)
        with self.assertRaises(ValueError):
            trace_response(1,np.inf,self.epoch)
        with self.assertRaises(ValueError):
            ResponseEpoch(self.act,H_j=np.nan)

    def test_example_is_not_cosmological_late300_prediction(self):
        e=example()
        self.assertEqual(len(e["frequency_scan"]),7)
        self.assertIn("full rolling G3",str(e["missing"]))
        self.assertLess(
           abs(e["step_trace_at_end"]-e["static_trace_target"]),.002)

if __name__=="__main__":
    unittest.main()
