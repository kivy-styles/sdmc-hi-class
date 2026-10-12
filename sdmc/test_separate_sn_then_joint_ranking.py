#!/usr/bin/env python3
import unittest
import tempfile
from pathlib import Path
import pandas as pd
from run_separate_sn_then_joint_ranking import choose,run,DATASETS
class CheckSeparations(unittest.TestCase):
    def setUp(self):
        self.rows=pd.DataFrame([
            dict(family="lcdm",label="baseline",status="OK",H0=68.57,omega_b=.0224,omega_cdm=.1182,A_s=2.05e-9,n_s=.965,tau_reio=.045,
                 chi2_planck=1018.7,chi2_bao=13.4,chi2_pantheonplus=1406.2,chi2_union3=28.8,chi2_desy5=1650.6),
            dict(family="lcdm",label="snleader",status="OK",H0=67.54,omega_b=.0224,omega_cdm=.1182,A_s=2.05e-9,n_s=.965,tau_reio=.045,
                 chi2_planck=1090.,chi2_bao=20.,chi2_pantheonplus=1404.3,chi2_union3=27.0,chi2_desy5=1646.5),
            dict(family="sdmc",label="baseline",status="OK",H0=69.71,omega_b=.0221,omega_cdm=.123,A_s=2.1e-9,n_s=.963,tau_reio=.053,
                 chi2_planck=1021.,chi2_bao=13.16,chi2_pantheonplus=1406.12,chi2_union3=28.79,chi2_desy5=1650.57),
            dict(family="sdmc",label="snleader",status="OK",H0=68.69,omega_b=.0221,omega_cdm=.123,A_s=2.1e-9,n_s=.963,tau_reio=.053,
                 chi2_planck=1097.,chi2_bao=17.,chi2_pantheonplus=1404.21,chi2_union3=26.96,chi2_desy5=1646.45),
        ])
    def test_independent_sn_first(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);data=p/"input.csv";self.rows.to_csv(data,index=False)
            out=run(data,p,"sn")
            self.assertEqual(set(out["datasets"]),set(DATASETS))
            for k in DATASETS:
                self.assertEqual(out["datasets"][k]["results"]["sdmc"]["id"],"snleader")
                self.assertEqual(out["datasets"][k]["results"]["lcdm"]["id"],"snleader")
                self.assertNotIn("joint_",out["datasets"][k]["metric"])
    def test_joint_reorders_candidates(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);file=p/"input.csv";self.rows.to_csv(file,index=False)
            run(file,p,"sn")
            out=run(file,p,"joint")
            for k in DATASETS:
                self.assertEqual(out["datasets"][k]["results"]["sdmc"]["id"],"baseline")
                self.assertEqual(out["datasets"][k]["results"]["lcdm"]["id"],"baseline")
                self.assertIn("joint_",out["datasets"][k]["metric"])
                self.assertTrue((p/f"sn_{k}_rank.csv").exists())
                self.assertTrue((p/f"joint_{k}_rank.csv").exists())
    def test_missing_family_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);f=p/"in.csv";self.rows[self.rows.family=="lcdm"].to_csv(f,index=False)
            with self.assertRaises(ValueError):run(f,p,"sn")
if __name__=="__main__":unittest.main()
