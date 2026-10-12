#!/usr/bin/env python3
"""Stage 2: use a frozen reconstructed covariant action and evolve freely.

The action is built separately by run_late300_linear_covariant.py from the
stage-1 winning parameterized SDMC target. This is a test of dynamical
realisability, NOT action-independent prediction or a full model-family fit.
"""
from pathlib import Path
import json, subprocess, traceback
import numpy as np
from run_paired_sdmc_lcdm_search import Likelihoods,read_table,cmb_files,OUT

def covariance_free_ini(original,root):
    source=Path(original).read_text().splitlines()
    result=[]
    replaced=set()
    for row in source:
        strip=row.strip()
        if strip.startswith("gravity_model"):
            result.append("gravity_model = sdmc_v3_covariant_linear_audit")
            replaced.add("gravity")
        elif strip.startswith("parameters_smg"):
            result.append("parameters_smg = 0.0")
            replaced.add("parameters")
        elif strip.startswith("expansion_model") or strip.startswith("expansion_smg"):
            # These are now frozen inside the reconstructed field-equation model.
            pass
        elif strip.startswith("root ="):
            result.append(f"root = {root}");replaced.add("root")
        else:result.append(row)
    if replaced!={"gravity","parameters","root"}:
        raise RuntimeError("Missing covariant input substitutions: "+repr(replaced))
    return "\n".join(result)+"\n"

def report():
    root=OUT/"covariant_promotion"
    root.mkdir(parents=True,exist_ok=True)
    out=OUT/"covariant_replay_result.json"
    summary=json.loads((OUT/"summary.json").read_text())
    top=json.loads((OUT/"promotion_target.json").read_text())
    label=top["label"]
    cache=OUT/"sdmc_screen_cache"
    ini=root/"covariant_free.ini"
    prefix=str(root/"covariant_")
    ini.write_text(covariance_free_ini("output/linear_cov_target.ini",prefix))
    result={"model":"sdmc","candidate":label,
            "note":"Free numerical evolution of reconstructed frozen action; action was selected by target reconstruction.",
            "status":"FAIL","screened_candidate":top}
    try:
        cp=subprocess.run(["./class",str(ini)],capture_output=True,text=True,timeout=420)
        result["returncode"]=cp.returncode
        if cp.returncode!=0:raise RuntimeError(cp.stdout[-1200:]+cp.stderr[-1200:])
        bg=read_table(prefix+"00_background.dat")
        target=read_table(cache/(label+"_00_background.dat"))
        mask=(bg.z.to_numpy(float)<=100)
        zz=bg.z.to_numpy(float)[mask]
        errors={}
        for col in ("H [1/Mpc]","M*^2_smg","kin (D)"):
            if col in bg and col in target:
                v=bg.loc[mask,col].to_numpy(float)
                t=np.interp(zz,target.z.to_numpy(float),target[col].to_numpy(float))
                errors[col]=float(np.max(np.abs(v-t)/np.maximum(np.abs(t),1e-16)))
        result["background_max_relative_errors_z_le_100"]=errors
        def cl_compare(p1,p2):
            a=np.loadtxt(p1);b=np.loadtxt(p2)
            n=min(len(a),len(b));a=a[:n];b=b[:n]
            if not np.allclose(a[:,0],b[:,0]):raise RuntimeError("Cl ell grids differ")
            return {name:float(np.max(np.abs(a[:,j]-b[:,j]))/
                               max(float(np.max(np.abs(a[:,j]))),1e-30))
                    for j,name in [(1,"TT"),(2,"EE"),(3,"TE"),(5,"PP")]}
        result["cl_max_peak_normalized_deviation"]=cl_compare(
            cache/(label+"_00_cl_lensed.dat"),prefix+"00_cl_lensed.dat")
        like=Likelihoods(extra_sn=(len(summary["screen_likelihoods"])>3))
        th=read_table(prefix+"00_thermodynamics.dat")
        cmb,acal=like.planck(prefix+"00_cl_lensed.dat")
        sn=like.score_sn(bg)
        aco=like.acoustics(bg,th)
        bao=like.bao(bg,aco["rd_Mpc"])
        result["scored_covariant"]={"chi2_planck":cmb,"A_planck":acal,
           "chi2_bao":bao,**{f"chi2_{k}":v for k,v in sn.items()},**aco,
           "chi2_primary":cmb+bao+sn["pantheonplus"]}
        result["delta_chi2_primary_vs_parameterized"]=float(
            result["scored_covariant"]["chi2_primary"]-top["chi2_primary"])
        local=summary["best"]["lcdm"]
        result["delta_chi2_covariant_vs_best_screened_lcdm"]=float(
            result["scored_covariant"]["chi2_primary"]-local["chi2_primary"])
        D=bg["kin (D)"].to_numpy(float)
        cs=bg["c_s^2"].to_numpy(float)
        slip=bg["braiding_smg"].to_numpy(float)+2*bg["M2_running_smg"].to_numpy(float)
        result["stability"]={"min_D":float(np.min(D)),"min_cs2":float(np.min(cs)),
              "max_cs2":float(np.max(cs)),
              "max_abs_noslip":float(np.max(np.abs(slip)))}
        result["checks"]={
           "H_F_replay_small_error":all(errors.get(k,float("inf"))<.01 for k in
                 ("H [1/Mpc]","M*^2_smg")),
           "positive_D_and_sound_speed":bool(np.min(D)>0 and np.min(cs)>0),
           "bounded_sound_speed":bool(np.max(cs)<=1.000001),
           "not_new_microscopic_action":True
        }
        result["status"]="OK" if all(result["checks"][k] for k in
                    ("H_F_replay_small_error","positive_D_and_sound_speed","bounded_sound_speed")) else "CHECK_FAILED"
    except Exception as e:
        result["error"]=repr(e)
        result["traceback"]=traceback.format_exc()[-2300:]
    finally:
        out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
        print("COVARIANT_PROMOTION_RESULT",json.dumps({
            key:result.get(key) for key in ("candidate","status",
              "delta_chi2_primary_vs_parameterized",
              "delta_chi2_covariant_vs_best_screened_lcdm","error")},sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    result=report()
    if result["status"]!="OK":raise SystemExit(1)
