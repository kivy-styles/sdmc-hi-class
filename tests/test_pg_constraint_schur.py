"""Structural-trace Schur complement against direct coupled Hessian tests."""
import unittest
from math import isclose
import numpy as np
from sdmc.pg_constraint_schur import (
    TraceAction, trace_stiffness,effective_kernel,
    quasi_static_frequency_ratio,full_static_hessian,
    schur_complement_trace,analytic_effective_hessian,generic_schur,
    toy_jordan_response,direct_full_hessian_solution,diagnostics)
from sdmc.pg_quadratic_source import TraceCoefficients,integrated_static_kernel

class ConstraintSchurTest(unittest.TestCase):
    def setUp(self):
        self.a=TraceAction(K=1.3,c2=.8,mass2=.2,weta=.04,
                            a_j=.89,a_s=1.31)

    def test_schur_equals_analytic_for_scales(self):
        for k in [0.,.01,.1,1.,4.,30.]:
            H=schur_complement_trace(k,self.a,A=2.,D=.8)
            ref=analytic_effective_hessian(k,self.a,A=2.,D=.8)
            self.assertTrue(np.allclose(H,ref,rtol=1e-11,atol=1e-11))

    def test_full_direct_solution_agrees_with_integrating_out(self):
        for k in [.01,.3,2.,8.,100.]:
            sol=direct_full_hessian_solution(k,self.a,A=k*k+1,D=.4,J=2.)
            ref=toy_jordan_response(k,self.a,A=k*k+1,D=.4,J=2.)
            self.assertTrue(np.allclose(sol,[ref["phi"],ref["m"],ref["trace_ell"]],
                                           rtol=2e-11,atol=2e-11))

    def test_static_quadratic_hessian_from_direct_finite_difference(self):
        k=1.4;A=2.;D=.7
        w=trace_stiffness(k,self.a)
        def E(phi,m,ell):
            return (.5*A*phi**2+.5*D*m*m+
              .5*self.a.K*w["omega2"]*ell*ell+
              w["C"]*(m-phi+ell/3)**2)
        vals=[.21,-.11,.07]
        eps=1e-4;Hd=np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                x=np.array(vals)
                if i==j:
                    xp=x.copy();xm=x.copy();xp[i]+=eps;xm[i]-=eps
                    Hd[i,j]=(E(*xp)-2*E(*x)+E(*xm))/eps**2
                else:
                    a=x.copy();b=x.copy();c=x.copy();d=x.copy()
                    a[i]+=eps;a[j]+=eps;b[i]+=eps;b[j]-=eps
                    c[i]-=eps;c[j]+=eps;d[i]-=eps;d[j]-=eps
                    Hd[i,j]=(E(*a)-E(*b)-E(*c)+E(*d))/(4*eps**2)
        H=full_static_hessian(k,self.a,A,D)
        self.assertTrue(np.allclose(Hd,H,rtol=1e-6,atol=1e-7))

    def test_low_k_k4_kernel(self):
        a=self.a
        k=.0001
        ratio=effective_kernel(2*k,a)/effective_kernel(k,a)
        self.assertAlmostEqual(ratio,16.,delta=1e-5)

    def test_large_k_k2_kernel(self):
        a=self.a
        k=1e6
        expected=4.5*a.K*a.c2/a.a_j**2
        measured=effective_kernel(k,a)/k**2
        self.assertAlmostEqual(measured,expected,
                               delta=1e-7*expected)

    def test_old_implementation_kernel_agrees(self):
        a=self.a
        c=TraceCoefficients(bf2=3*a.K,c_l2=a.c2,
                   m_l2=a.mass2,w_eta=a.weta,a_j=a.a_j,a_s=a.a_s)
        for k in [0.,.002,.4,1.,4.,20.]:
            self.assertAlmostEqual(effective_kernel(k,a),
                                   integrated_static_kernel(k,c),
                                   delta=1e-10*max(1,effective_kernel(k,a)))

    def test_zero_curvature_does_not_modify_metric(self):
        a=TraceAction(K=1.,mass2=.2,weta=0.)
        for k in [.1,2.,20.]:
            self.assertEqual(effective_kernel(k,a),0.)
            self.assertAlmostEqual(toy_jordan_response(k,a,1.,.7)["toy_screening_ratio"],1.)

    def test_mapper_free_relaxes_relative_strain(self):
        for k in [.1,1.,20.]:
            out=toy_jordan_response(k,self.a,A=2.,D=0.,J=3.)
            self.assertAlmostEqual(out["toy_screening_ratio"],1.)
            self.assertAlmostEqual(out["m"],out["phi"])

    def test_locked_mapper_adds_stiffness(self):
        k=2.
        s=toy_jordan_response(k,self.a,A=1.5,D=1e12)
        self.assertAlmostEqual(s["delta_jordan_stiffness"],
                               2*effective_kernel(k,self.a),places=9)

    def test_positive_toy_stiffness_cannot_enhance_attraction(self):
        for k in [.01,.1,1.,4.,30.]:
            for D in [0.,.01,.5,2.,1e6]:
                r=toy_jordan_response(k,self.a,A=k*k+.1,D=D)
                self.assertLessEqual(r["toy_screening_ratio"],1.+1e-12)
                self.assertGreater(r["toy_screening_ratio"],0.)

    def test_generic_schur_identity_for_arbitrary_coupled_blocks(self):
        rng=np.random.default_rng(50)
        x=rng.normal(size=(5,5))
        H=x.T@x+np.eye(5)
        ans=generic_schur(H,[0,1,4])
        src=np.array([.2,-.15,.6,0.,0.])
        full=np.linalg.solve(H,src)
        # need source retained indices [0,1,4], no source on 2,3
        src=np.array([.2,-.15,0.,0.,.6])
        full=np.linalg.solve(H,src)
        red=np.linalg.solve(ans,src[[0,1,4]])
        self.assertTrue(np.allclose(full[[0,1,4]],red,
                                    rtol=1e-12,atol=1e-12))

    def test_singular_eliminated_constraint_is_not_pseudoinverted(self):
        H=np.array([[1.,0.],[0.,0.]])
        with self.assertRaises(np.linalg.LinAlgError):
            generic_schur(H,[0])

    def test_quasi_static_validity_is_an_extra_condition(self):
        a=self.a
        self.assertGreater(quasi_static_frequency_ratio(20.,a,.01),
                           quasi_static_frequency_ratio(.01,a,.01))
        self.assertLess(quasi_static_frequency_ratio(0.,a,10.),1.)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):effective_kernel(-1.,self.a)
        with self.assertRaises(ValueError):TraceAction(K=-1.)
        with self.assertRaises(ValueError):toy_jordan_response(1.,self.a,A=0.,D=2.)
        with self.assertRaises(ValueError):generic_schur(np.ones((2,3)),[0])

    def test_demo_requires_independent_jordan_hessian(self):
        d=diagnostics()
        self.assertIn("complete F/G2/G3 Hessian",d["missing"])
        self.assertEqual(len(d["k_sweep"]),5)
        self.assertLess(d["k_sweep"][-1]["toy_screening_ratio"],1.)

if __name__=="__main__":unittest.main()
