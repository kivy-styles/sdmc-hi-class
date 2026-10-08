#!/usr/bin/env python3
import sys, json, math, re
from pathlib import Path
import numpy as np

OUT=Path("output")
OUT.mkdir(exist_ok=True)
CAND=json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]

def prepare():
    xi_inv=math.sqrt(3/(8*math.pi))
    dcanon=16/(CAND["lambda_e"]**2)
    cases={
        "accepted":(CAND["D0"],CAND["D_floor"],1.0),
        "xiD0":(xi_inv,CAND["D_floor"],1.0),
        "canonFloor":(CAND["D0"],dcanon,1.0),
        "structural":(xi_inv,dcanon,1.0),
    }
    rows=[
      f"H0 = {CAND['H0']}",
      f"omega_b = {CAND['omega_b']}",
      f"omega_cdm = {CAND['omega_cdm']}",
      "N_ncdm = 0","N_ur = 3.046","T_cmb = 2.7255","YHe = 0.2453",
      f"A_s = {CAND['A_s']}",f"n_s = {CAND['n_s']}",f"tau_reio = {CAND['tau_reio']}",
      "Omega_Lambda = 0","Omega_fld = 3.1443554e-8",
      "fluid_equation_of_state = SDMC_TRACKER","cs2_fld = 0.003","use_ppf = no",
      "Omega_smg = -1","gravity_model = sdmc_v3_independent_kinetic",
      "expansion_model = sdmc_full",
      f"expansion_smg = {CAND['Omega_x0']}, {CAND['lambda_e']}, {CAND['z_t']}, 0.5, {CAND['A_late']}, 0.25, {CAND['B_late']}, 1.5",
      "pert_initial_conditions_smg = zero","method_qs_smg = fully_dynamic",
      "a_ini_over_a_today_default = 1.e-8","a_ini_test_qs_smg = 1.e-8",
      "pert_ic_ini_z_ref_smg = 1.e7","a_min_stability_test_smg = 1.e-8",
      "output_background_smg = 3","modes = s","output = tCl,pCl,lCl,mPk",
      "lensing = yes","l_max_scalars = 2500","P_k_max_h/Mpc = 1.0","z_pk = 0",
      "write background = yes","write thermodynamics = no","format = class",
      "input_verbose = 0","background_verbose = 0","thermodynamics_verbose = 0",
      "perturbations_verbose = 0","transfer_verbose = 0","primordial_verbose = 0",
      "harmonic_verbose = 0","fourier_verbose = 0","lensing_verbose = 0","output_verbose = 0",
    ]
    base="\n".join(rows)+"\n"
    meta={"Xi_v_inverse":xi_inv,"canonical_radiation_floor":dcanon,"cases":{}}
    for tag,(d0,df,powr) in cases.items():
        s=base+f"parameters_smg = {CAND['A_F']}, {CAND['z_c']}, {CAND['width']}, {d0}, {powr}, {df}\nroot = output/{tag}_\n"
        (OUT/f"{tag}.ini").write_text(s)
        meta["cases"][tag]={"D0":d0,"Dfloor":df,"power":powr}
    (OUT/"structural_kinetic_inputs.json").write_text(json.dumps(meta,indent=2)+"\n")
    print("STRUCTURAL_KINETIC_PREPARED",json.dumps(meta,sort_keys=True))

def read_bg(tag):
    p=OUT/f"{tag}_00_background.dat"
    ls=p.read_text().splitlines()
    hdr=[x for x in ls if x.startswith("#") and re.search(r"(?:^|\s)1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(p)
    if a.ndim==1:a=a[None,:]
    return {n:a[:,i] for i,n in enumerate(names)}

def compare(pa,pb):
    A=np.loadtxt(pa); B=np.loadtxt(pb)
    if A.ndim==1:A=A[None,:]
    if B.ndim==1:B=B[None,:]
    n=min(A.shape[1],B.shape[1])
    x=A[:,0]; xb=B[:,0]
    vals=[]
    for j in range(1,n):
        o=np.argsort(xb)
        bb=np.interp(x,xb[o],B[o,j]); aa=A[:,j]
        peak=max(float(np.max(np.abs(aa))),float(np.max(np.abs(bb))),1e-300)
        vals.append(float(np.max(np.abs(aa-bb))/peak))
    return max(vals) if vals else None

def analyze():
    out={"status":"no-fit structural kinetic hypothesis test","cases":{}}
    for tag in ["accepted","xiD0","canonFloor","structural"]:
        rc=int((OUT/f"{tag}.rc").read_text())
        rec={"returncode":rc}
        if rc==0:
            d=read_bg(tag); z=d["z"]; m=z<=100
            rec["stability"]={
              "min_D_z100":float(np.min(d["kin (D)"][m])),
              "min_cs2_z100":float(np.min(d["c_s^2"][m])),
              "max_cs2_z100":float(np.max(d["c_s^2"][m])),
              "max_abs_noslip_z100":float(np.max(np.abs(d["braiding_smg"][m]+2*d["M2_running_smg"][m])))
            }
            if tag!="accepted":
                rec["observable_residuals"]={
                  "unlensed_cl_peaknorm":compare(OUT/"accepted_00_cl.dat",OUT/f"{tag}_00_cl.dat"),
                  "lensed_cl_peaknorm":compare(OUT/"accepted_00_cl_lensed.dat",OUT/f"{tag}_00_cl_lensed.dat"),
                  "pk_peaknorm":compare(OUT/"accepted_00_pk.dat",OUT/f"{tag}_00_pk.dat")
                }
        else:
            log=(OUT/f"{tag}.log").read_text(errors="replace")
            rec["error_tail"]=[x for x in log.splitlines() if "error" in x.lower() or "instability" in x.lower()][-12:]
        out["cases"][tag]=rec
    (OUT/"structural_kinetic_test.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("STRUCTURAL_KINETIC_TEST",json.dumps(out,sort_keys=True))

if __name__=="__main__":
    if len(sys.argv)!=2 or sys.argv[1] not in ("prepare","analyze"):
        raise SystemExit("usage: run_structural_kinetic_hypothesis.py prepare|analyze")
    globals()[sys.argv[1]]()
