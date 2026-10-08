#!/usr/bin/env python3
"""Repair the historical D0/Dfloor covariant tests with genuine live-input reconstruction.

The previously archived full-Planck/DESI replacement jobs changed the
parameterized target but rebuilt its covariant action with hardcoded original
kinetic coordinates. The workflow replays a validated late300 source recipe
whose linear-G3 reconstruction reads D0 and Dfloor from the generated target
INI. The kinetic D transition window remains distinct from the force F window.

It compares live target, reconstructed free action, perturbation spectra and
optionally a true Planck high-l Plik-lite likelihood. This is NOT a proof of
a microscopic D origin or a full joint cosmological Bayesian evidence.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re, subprocess, sys
from pathlib import Path
import numpy as np
import yaml

SOURCE=Path(".github/workflows/sdmc_late300_derived_bundle_covariant.yml")
OUT=Path("output")
FORCE_AF=0.020520
FORCE_ZC=4.034077502899133
FORCE_WIDTH=0.32925377073762596
D_ZC=3.927876388467848
D_WIDTH=0.33114133956842123
EXPANSION="0.7014079815036496, 18.40625, 16.742289660253434, 0.5, 0.024822032167576252, 0.25, 0.01264704344701022, 1.5"
BASE_INI="parameters_smg = 0.02052, 4.034077502899133, 0.32925377073762596, 0.3419932122356156, 1.0, 0.05181542627513409"
CONFIGS={
 "derived_reference":(.3419932122356156,.05181542627513409),
 "D0_2pi_over_lambda":(.34136151074659893,.05181542627513409),
 "Dfloor_tracker":(.3419932122356156,.047226890271848634),
 "composite_2pi_tracker":(.34136151074659893,.047226890271848634),
 "D0_Xi_over_sqrtF":(.342002216921491,.05181542627513409),
 "historical_D0":(.34231919445927034,.05181542627513409),
}

def run(cmd,name):
 print("=== TRUE_INPUT_REPLAY_STEP "+name+" ===",flush=True)
 subprocess.run(cmd,check=True)

def read_bg(path):
 rows=path.read_text().splitlines()
 hdr=[x for x in rows if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
 mark=list(re.finditer(r"(\d+)\s*:\s*",hdr))
 keys=[hdr[m.end():mark[j+1].start() if j+1<len(mark) else len(hdr)].strip() for j,m in enumerate(mark)]
 arr=np.atleast_2d(np.loadtxt(path))
 return {key:arr[:,j] for j,key in enumerate(keys)}

def get_step(steps,name):
 matches=[x for x in steps if x.get("name")==name and "run" in x]
 if len(matches)!=1:raise RuntimeError("Expected one executable source step "+name)
 return matches[0]["run"]

def run_src(script,name):
 run(["bash","-e","-o","pipefail","-c",script],name)

def Dswitch(z):
 N=-np.log1p(z)
 Nc=-math.log1p(D_ZC)
 return 1./(1.+np.exp(-np.clip((N-Nc)/D_WIDTH,-700,700)))

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--case",required=True,choices=CONFIGS)
 ap.add_argument("--planck",action="store_true")
 args=ap.parse_args()
 case=args.case
 D0,DFLOOR=CONFIGS[case]
 source=SOURCE.read_text()
 if source.count(BASE_INI)!=1:raise RuntimeError("Original late300 target kinetic coordinates absent or duplicated")
 if source.count("Dtarget=DFLOOR+D0*S**POWER")!=1:raise RuntimeError("Known linear-G3 kinetic patch missing")
 srcsha=hashlib.sha256(source.encode()).hexdigest()
 source_steps=yaml.safe_load(source)["jobs"]["linear-covariant-audit"]["steps"]
 OUT.mkdir(exist_ok=True)
 info={"case":case,"D0":D0,"Dfloor":DFLOOR,"power":1,
       "F_window":{"AF":FORCE_AF,"zc":FORCE_ZC,"width":FORCE_WIDTH},
       "D_window":{"zc":D_ZC,"width":D_WIDTH},
       "source_workflow":str(SOURCE),"source_sha256":srcsha,
       "source_reconstruction":"sdmc/run_late300_linear_covariant.py",
       "status":"not yet validated"}
 (OUT/"D_true_inputs_provenance.json").write_text(json.dumps(info,indent=2)+"\n")
 targetline=f"parameters_smg = {FORCE_AF}, {FORCE_ZC}, {FORCE_WIDTH}, {D0}, 1.0, {DFLOOR}"
 run_src(get_step(source_steps,"Apply SDMC parameterized patches"),"apply patches and compile baseline")
 target=get_step(source_steps,"Generate late300 parameterized target")
 if target.count(BASE_INI)!=1:raise RuntimeError("Unique native target kinetic anchor missing")
 run_src(target.replace(BASE_INI,targetline),"generate new target for "+case)
 ini=(OUT/"linear_cov_target.ini").read_text()
 if ini.count(targetline)!=1 or EXPANSION not in ini:raise RuntimeError("Wrong kinetic or expansion INI")
 reconstruct=get_step(source_steps,"Reconstruct and install derived-bundle linear-G3 action")
 if "3.927876388467848" not in reconstruct or "0.33114133956842123" not in reconstruct:
  raise RuntimeError("Separate historical D-window omitted from reconstruction")
 run_src(reconstruct,"regenerate G2,G3,G4 using current target INI")
 summary=dict(line.split("=",1) for line in (OUT/"linear_cov_reconstruction_summary.txt").read_text().splitlines() if "=" in line and "[" not in line)
 if abs(float(summary["D0"])-D0)>1e-13 or abs(float(summary["Dfloor"])-DFLOOR)>1e-13:
  raise RuntimeError("CRITICAL: reconstructed G2 used stale benchmark D inputs")
 if float(summary["Dclosure_max_rel"])>1e-7:raise RuntimeError("Covariant kinetic identity failed")
 run_src(get_step(source_steps,"Free evolve reconstructed linear-G3 action"),"free-evolve NEW covariant action")
 run_src(get_step(source_steps,"Compare free evolution to target"),"compare homogeneous dynamics")
 run_src(get_step(source_steps,"Compare full Boltzmann observables"),"compare TT/TE/EE and P(k)")
 t=read_bg(OUT/"linear_cov_target_00_background.dat")
 f=read_bg(OUT/"linear_cov_free_00_background.dat")
 Dexpected=float(DFLOOR+D0*Dswitch(0.))
 Dt=float(t["kin (D)"][np.argmin(abs(t["z"]))])
 Df=float(f["kin (D)"][np.argmin(abs(f["z"]))])
 if abs(Dexpected-Dt)>1e-8:raise RuntimeError(f"Native kinetic output {Dt} differs from requested {Dexpected}")
 if abs(Df-Dt)>2e-5:raise RuntimeError(f"Free action kinetic {Df} differs from target {Dt}")
 cmp=json.loads((OUT/"linear_cov_replay_comparison.json").read_text())
 boltz=json.loads((OUT/"linear_cov_boltzmann_comparison.json").read_text())
 if cmp["max_abs_dH_H_z100"]>1e-4 or cmp["max_abs_dF_F_z100"]>1e-4:
  raise RuntimeError("Free action and target background do not agree")
 if cmp["min_cs2_z100"]<=0 or cmp["min_D_z100"]<=0 or cmp["max_cs2_z100"]>1:
  raise RuntimeError("Candidate-specific covariant action violates tested scalar gate")
 TT=boltz.get("Cl_TT_max_norm_abs")
 if TT is None or TT>0.005:raise RuntimeError(f"Covariant TT replay check failed {TT}")
 for fn in ("linear_cov_target_00_cl_lensed.dat","linear_cov_free_00_cl_lensed.dat"):
  spec=np.loadtxt(OUT/fn)
  if not np.all(np.isfinite(spec)):raise RuntimeError("Nonfinite "+fn)
 info.update({"status":"PASSED true-input action reproduction",
        "D_today_expected":Dexpected,"D_today_target":Dt,"D_today_covariant":Df,
        "k1_today":float(summary["k1_0"]),"k2X_today":float(summary["k2X_0"]),
        "F_today":float(summary["F0_target"]),
        "D_identity_relative_max":float(summary["Dclosure_max_rel"]),
        "H_replay_relative_max_z100":cmp["max_abs_dH_H_z100"],
        "F_replay_relative_max_z100":cmp["max_abs_dF_F_z100"],
        "TT_replay_peak_normalized_max":TT,
        "min_cs2_z100":cmp["min_cs2_z100"],"max_cs2_z100":cmp["max_cs2_z100"]})
 (OUT/"D_true_inputs_provenance.json").write_text(json.dumps(info,indent=2,sort_keys=True)+"\n")
 print("D_TRUE_INPUTS_ACTION_REPLAY",json.dumps(info,sort_keys=True),flush=True)
 if args.planck:
  code="from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native\nok=TTTEEE_lite_native.install(path='planck_packages',no_progress_bars=True)\nif not ok: raise SystemExit('Plik-lite installation failed')\n"
  run([sys.executable,"-c",code],"install genuine high-l Planck TTTEEE plik-lite data")
  run([sys.executable,"sdmc/evaluate_planck_pliklite.py","--packages","planck_packages",
       "--output",str(OUT/"D_true_inputs_planck_pliklite.csv"),
       f"{case}_TARGET=output/linear_cov_target_00_cl_lensed.dat",
       f"{case}_COVARIANT=output/linear_cov_free_00_cl_lensed.dat"],"Plik-lite on reconstructed covariant spectra")
  with (OUT/"D_true_inputs_planck_pliklite.csv").open() as fcsv:
   rows=list(csv.DictReader(fcsv))
  info["planck_pliklite"]=rows
  (OUT/"D_true_inputs_provenance.json").write_text(json.dumps(info,indent=2,sort_keys=True)+"\n")
  print("D_TRUE_INPUTS_PLANCK",json.dumps({"case":case,"rows":rows},sort_keys=True),flush=True)
 print("D_TRUE_INPUT_REPLAY_COMPLETE",case,D0,DFLOOR,flush=True)

if __name__=="__main__":main()
