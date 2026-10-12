#!/usr/bin/env python3
"""Scientific guardrails for paired dynamic candidate search."""
import unittest
import numpy as np
from run_paired_sdmc_lcdm_search import (LCDM_BASE,SDMC_BASE,design,
    flat_omega,build_ini,NAMES,SCALE,STRUCT_AXIS)

class TestPairedSearch(unittest.TestCase):
    def test_matched_evaluation_budget(self):
        for budget in (1,9,17,25):
            a=list(design("lcdm",budget))
            b=list(design("sdmc",budget))
            self.assertEqual(len(a),budget)
            self.assertEqual(len(b),budget)
            self.assertEqual(a[0][0],"baseline")
            self.assertEqual(b[0][0],"baseline")

    def test_both_benchmark_parameters(self):
        self.assertAlmostEqual(LCDM_BASE["H0"],68.56858744695782)
        self.assertAlmostEqual(SDMC_BASE["H0"],69.71482083084993)
        self.assertAlmostEqual(SDMC_BASE["A_F"],.02048146490100771)

    def test_identical_ordinary_sobol_offsets(self):
        a=list(design("lcdm",17))
        b=list(design("sdmc",17))
        for (la,pa),(lb,pb) in zip(a[1:],b[1:]):
            self.assertEqual(la,lb)
            for k in NAMES:
                u=pa[k]-LCDM_BASE[k];v=pb[k]-SDMC_BASE[k]
                self.assertLess(abs(u-v),max(1e-16,abs(u)*1e-12))

    def test_structural_directions_only_in_sdmc(self):
        a=list(design("lcdm",25));b=list(design("sdmc",25))
        self.assertEqual(len(a),25);self.assertEqual(len(b),25)
        self.assertTrue(all("axis_" in label for label,p in a[17:]))
        self.assertTrue(all(any(name in label for name,step in STRUCT_AXIS)
                            for label,p in b[17:]))
        self.assertTrue(all(p["D0"]==SDMC_BASE["D0"] for _,p in b))

    def test_flat_closure_exact_to_primitive_densities(self):
        for p in (LCDM_BASE,SDMC_BASE):
            h=p["H0"]/100.
            omega=flat_omega(p)
            self.assertAlmostEqual(omega+(p["omega_b"]+p["omega_cdm"]+
                                   4.17998772e-5)/h**2,1.,places=12)

    def test_class_evolution_models_are_distinct(self):
        li=build_ini("lcdm",LCDM_BASE,"output/control_")
        si=build_ini("sdmc",SDMC_BASE,"output/target_")
        self.assertIn("Omega_Lambda =",li)
        self.assertNotIn("gravity_model =",li)
        self.assertIn("gravity_model = sdmc_v3_independent_kinetic",si)
        self.assertIn("method_qs_smg = fully_dynamic",si)
        self.assertNotIn("sdmc_v3_covariant_linear_audit",si)

    def test_planck_LCDM_approximation_is_not_exact_H0_inference(self):
        # Planck 2018 flat LCDM geometric-degeneracy approximation:
        # Ωm h^3 = (ωb+ωc)h ≈ 0.09633, with peak-shape/frequency dependence.
        omega_m_physical=.02237+.1200
        H0_approx=100*.09633/omega_m_physical
        self.assertGreater(H0_approx,65.)
        self.assertLess(H0_approx,70.)
        self.assertNotAlmostEqual(H0_approx,LCDM_BASE["H0"],places=2)

    def test_no_double_supernova_counting_in_primary_objective(self):
        # Source-level guardrail: primary score is Planck + BAO + Pantheon+,
        # while Union3/DESY5 are *separate* alternative joint scores.
        import inspect
        from run_paired_sdmc_lcdm_search import evaluate
        s=inspect.getsource(evaluate)
        self.assertIn('planck+bao+sn["pantheonplus"]',s)
        self.assertIn('("union3","desy5")',s)

if __name__=="__main__":
    unittest.main()
