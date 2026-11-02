"""SDMC nonflat handoff: action constraint + decelerating transfer checks."""
import unittest
import numpy as np
from scipy.integrate import simpson
from sdmc.flatness_nonflat_handoff import (
 ClassicalTransfer,wx,gamma,transfer_rhs,classical_transfer,
 curved_classical_history,nonflat_restricted_constraint,
 constraint_at_x,archive_clock_check,transition_curvature_residual,
 endpoint_rounding_test,diagnostic)

class NonflatJordanEmergence(unittest.TestCase):
    def setUp(self):
        self.p=ClassicalTransfer()

    def test_archive_lapse_suppression_was_provisional_and_not_K(self):
        r=archive_clock_check()
        self.assertAlmostEqual(r["crossover_growth"],2.25)
        self.assertAlmostEqual(r["thermal_handoff_suppression"],.0381431,
                               delta=3e-7)
        self.assertAlmostEqual(r["total_structural_suppression"],.0858219,
                               delta=3e-7)

    def test_transition_stiff_to_radiation_w_and_positive_thermal_rate(self):
        self.assertGreater(wx(0,self.p),.99)
        self.assertLess(abs(wx(4,self.p)-1/3),.0001)
        self.assertGreater(gamma(1.2,self.p),4.99)
        self.assertGreaterEqual(gamma(0,self.p),0)

    def test_transfer_conserves_total_local_energy_no_created_mass(self):
        for x in [0,.3,.8,1.3,2,3,4]:
            X=1.23;R=.45
            dX,dR=transfer_rhs(x,[X,R],self.p)
            w=wx(x,self.p)
            self.assertAlmostEqual(dX+dR,
                           -3*(1+w)*X-4*R,places=12)

    def test_radiation_is_generated_and_eventually_dominates(self):
        history=classical_transfer()
        self.assertGreater(history[-1]["rho_r"],history[-1]["rho_X"])
        self.assertGreater(history[2]["rho_r"],1e-8)
        for i in history:
            self.assertGreater(i["rho_total"],0)

    def test_deceleration_decreases_while_curvature_fraction_grows(self):
        h=curved_classical_history(K=-1e-6)
        self.assertGreater(h[0]["deceleration"],1.98)
        self.assertLess(abs(h[-1]["deceleration"]-1),.05)
        self.assertTrue(all(t["deceleration"]>0 for t in h))
        self.assertTrue(all(abs(h[i+1]["Omega_k"])>
                         abs(h[i]["Omega_k"]) for i in range(len(h)-1)))
        self.assertGreater(abs(h[-1]["Omega_k"])/abs(h[0]["Omega_k"]),7000)

    def test_both_nonflat_signs_grow_in_abs_when_expanding(self):
        for K in (-1e-6,1e-6):
            h=curved_classical_history(K)
            self.assertGreater(abs(h[-1]["Omega_k"]),abs(h[0]["Omega_k"]))
            self.assertTrue(all(t["H2"]>0 for t in h))
            for state in h:
                self.assertAlmostEqual(state["K_reconstructed"],K,
                        delta=1e-8*abs(K))

    def test_empty_curvature_stays_zero(self):
        h=curved_classical_history(K=0.)
        self.assertTrue(all(item["Omega_k"]==0 for item in h))
        self.assertTrue(all(item["deceleration"]>0 for item in h))

    def test_exact_curvature_log_derivative_matches_integral_2q(self):
        x=np.linspace(0,4,2001)
        h=curved_classical_history(-1e-6,positions=x)
        logratio=np.log(abs(h[-1]["Omega_k"])/abs(h[0]["Omega_k"]))
        integral=2*simpson([r["deceleration"] for r in h],x=x)
        self.assertAlmostEqual(logratio,integral,delta=2e-5)

    def test_selector_action_restrictive_constraint_has_two_roots(self):
        p=nonflat_restricted_constraint(F=1/3,H=1,E0=1,W=1,
                             eta=.01,aJ=1,M=1)
        self.assertEqual(len(p["roots_x"]),2)
        self.assertAlmostEqual(p["roots_x"][0],0,places=11)
        self.assertAlmostEqual(p["roots_x"][1],50/3,places=11)
        self.assertFalse(p["curvature_dynamically_selected"])
        for x in p["roots_x"]:
            self.assertAlmostEqual(constraint_at_x(1/3,1,1,1,.01,1,1,x),0.,
                        places=10)

    def test_nonzero_residual_predicts_conditional_nonzero_x(self):
        params=dict(F=1/3,H=1,E0=1.002,W=1,eta=.01,aJ=1,M=1)
        p=nonflat_restricted_constraint(**params)
        self.assertEqual(len(p["roots_x"]),2)
        self.assertGreater(p["roots_x"][0],0.)
        self.assertLess(p["roots_x"][0],.003)
        for x in p["roots_x"]:
            self.assertAlmostEqual(constraint_at_x(**params,x=x),0.,
                 delta=1e-10)

    def test_uv_off_restores_linear_standard_constraint(self):
        p=nonflat_restricted_constraint(F=1/3,H=1,E0=1.002,
                 W=0,eta=.01,aJ=1,M=1)
        self.assertEqual(len(p["roots_x"]),1)
        self.assertAlmostEqual(p["roots_x"][0],.002,places=12)

    def test_discriminant_excludes_some_prescribed_values(self):
        p=nonflat_restricted_constraint(F=1/3,H=1,E0=10,
                          W=1,eta=1,aJ=1,M=1)
        self.assertFalse(p["real_solution"])

    def test_two_different_K_are_valid_for_same_transfer(self):
        open_model=curved_classical_history(-1e-6)
        closed_model=curved_classical_history(1e-6)
        for a,b in zip(open_model,closed_model):
            self.assertAlmostEqual(a["rho_total"],b["rho_total"],
                delta=1e-13)
            self.assertNotEqual(a["H2"],b["H2"])
            self.assertAlmostEqual(a["K_reconstructed"],-1e-6,
                    delta=1e-15)
            self.assertAlmostEqual(b["K_reconstructed"],1e-6,
                    delta=1e-15)

    def test_excessive_positive_K_turns_around(self):
        with self.assertRaises(ValueError):
            curved_classical_history(K=.01)

    def test_invalid_data_rejected(self):
        with self.assertRaises(ValueError):
            ClassicalTransfer(gamma_width=0.)
        with self.assertRaises(ValueError):
            classical_transfer(positions=[2,1,3])
        with self.assertRaises(ValueError):
            nonflat_restricted_constraint(F=0,H=1,E0=1,
                  W=1,eta=1,aJ=1,M=1)

    def test_endpoint_constraint_from_N_p_chi_f_is_not_a_fitted_Ck(self):
        d=endpoint_rounding_test()
        self.assertAlmostEqual(
             d["rounded"]["omega_k_J"],4.301159775177421e-5,
             delta=5e-14)
        self.assertAlmostEqual(d["exact_flat"]["omega_k_J"],0.,
                               delta=1e-14)
        self.assertFalse(d["rounded"]["independently_derived_curvature"])
        self.assertAlmostEqual(d["flat_chi_required"],
            .0034771495577581015,delta=1e-14)

    def test_curvature_sign_follows_independent_source_residual(self):
        GF=1/(8*np.pi)
        N=3.41350846; p=2.; f=.01
        flat=(3*GF*N*N/(p*p))*f
        for relative_perturbation in (-.002,0,.002):
            r=transition_curvature_residual(
                 N,p,chi=flat*(1+relative_perturbation),
                 f=f,GF=GF)
            self.assertAlmostEqual(r["omega_k_J"],
                        -relative_perturbation,delta=2e-14)

    def test_same_omega_curvature_from_gamma_not_unique_microphysics(self):
        GF=1/(8*np.pi);N=3.41350846;p=2;f=.01
        flat=3*GF*N*N/p**2
        chi=.003
        gamma=flat-chi/f
        a=transition_curvature_residual(N,p,chi,f,GF,gamma_g=gamma)
        b=transition_curvature_residual(N,p,flat*f,f,GF)
        self.assertAlmostEqual(a["omega_k_J"],b["omega_k_J"],delta=1e-14)
        self.assertGreater(gamma,0)
        # Thus metric curvature cannot by itself tell apart
        # active scalar and supporting-source partitions.

    def test_invalid_zero_radiation_fraction_rejected(self):
        with self.assertRaises(ValueError):
            transition_curvature_residual(3,2,.01,0,1)
        with self.assertRaises(ValueError):
            transition_curvature_residual(3,2,.01,1.01,1)

    def test_qualified_no_claim_of_genuine_primordial_selection(self):
        r=diagnostic()
        self.assertIn("not full nonflat SDMC",r["scope"])
        self.assertTrue(any("prethermal" in x for x in r["missing"]))
        self.assertFalse(r["flat_residual_selector_on"][
             "curvature_dynamically_selected"])

if __name__=="__main__":
    unittest.main()
