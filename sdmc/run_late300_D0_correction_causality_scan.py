#!/usr/bin/env python3
"""Physical lower bound on the D0 *correction* using c_s^2<=1 and background scan.

Uses the same derived F-window / handoff / D_floor as the passed AF-D0 map.
Runs real hi_class at independent D0 values, samples D and c_s^2, and checks
whether N_s = D*c_s^2 is approximately invariant when D0 is changed.
This is not an independent microscopic prediction of the correction amplitude.
"""
import json
import math
import re
import subprocess
from pathlib import Path
import numpy as np

OUT=Path("output/D0_correction_causality");OUT.mkdir(parents=True,exist_ok=True)
AF=0.020520
ZC=4.034077502899133
WIDTH=0.32925377073762596
DFLOOR=0.05181542627513409
ZT=16.742289660253434
D0_ACCEPTED=0.34231919445927034
D0_VALUES=[0.0005,0.005,0.015,0.018,0.019,0.01915,0.01917,0.01918,0.0192,0.01925,0.0195,0.02,0.025,0.03,0.05,0.075,0.10,0.15,0.20,0.25,0.30,D0_ACCEPTED,0.38]
INI="""H0 = 69.71482083084993
omega_b = 0.022083219194622913
omega_cdm = 0.12299536722293603
N_ncdm = 0
N_ur = 3.046
T_cmb = 2.7255
YHe = 0.2453
A_s = 2.1047019076951746e-09
n_s = 0.9628416801318527
tau_reio = 0.053167105823755265
Omega_Lambda = 0
Omega_fld = 3.1443554e-8
fluid_equation_of_state = SDMC_TRACKER
cs2_fld = 0.003
use_ppf = no
Omega_smg = -1
gravity_model = sdmc_v3_independent_kinetic
expansion_model = sdmc_full
expansion_smg = 0.7014079815036496, 18.40625, 16.742289660253434, 0.5, 0.024822032167576252, 0.25, 0.01264704344701022, 1.5
pert_initial_conditions_smg = zero
method_qs_smg = fully_dynamic
a_ini_over_a_today_default = 1.e-8
a_ini_test_qs_smg = 1.e-8
pert_ic_ini_z_ref_smg = 1.e7
a_min_stability_test_smg = 1.e-8
output_background_smg = 3
modes = s
output =
write background = yes
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
output_verbose = 0
"""

def read_bg(path):
    lines=path.read_text().splitlines()
    hdr=[q for q in lines if q.startswith("#") and re.search(r"(?:^|\s)1\s*:",q)][-1].lstrip("#").strip()
    ma=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    keys=[hdr[m.end():ma[i+1].start() if i+1<len(ma) else len(hdr)].strip() for i,m in enumerate(ma)]
    dat=np.atleast_2d(np.loadtxt(path))
    return {k:dat[:,i] for i,k in enumerate(keys)}

def S_for_z(z):
    N=-np.log1p(z)
    Nc=-np.log1p(ZC)
    u=(N-Nc)/WIDTH
    return 1/(1+np.exp(-np.maximum(-700,np.minimum(700,u))))

def main():
    records=[]
    curves={}
    for d0 in D0_VALUES:
        tag=f"d0_{d0:.9f}".replace(".","p")
        root=OUT/(tag+"_")
        ini=OUT/(tag+".ini")
        ini.write_text(INI+f"parameters_smg = {AF}, {ZC}, {WIDTH}, {d0}, 1.0, {DFLOOR}\nroot = {root}\n")
        p=subprocess.run(["./class",str(ini)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        rec={"D0":d0,"returncode":p.returncode}
        bg=Path(str(root)+"00_background.dat")
        if bg.exists():
            try:
                cols=read_bg(bg)
                z=np.asarray(cols["z"],float); D=np.asarray(cols["kin (D)"],float)
                c=np.asarray(cols["c_s^2"],float)
                mask=(z<=100.)&np.isfinite(c)&np.isfinite(D)
                j=np.argmin(np.abs(z))
                N=-np.log1p(z[mask])
                order=np.argsort(N)
                curves[d0]=(N[order], (D[mask]*c[mask])[order])
                full=(np.isfinite(c)&np.isfinite(D))
                rec.update({
                    "min_cs2_all":float(np.min(c[full])),
                    "max_cs2_all":float(np.max(c[full])),
                    "subluminal_all":bool(np.min(c[full])>=0 and np.max(c[full])<=1),
                    "z_max_cs2_all":float(z[full][np.argmax(c[full])]),
                    "min_cs2_z100":float(np.min(c[mask])),
                    "max_cs2_z100":float(np.max(c[mask])),
                    "min_D_z100":float(np.min(D[mask])),
                    "D_z0":float(D[j]),
                    "cs2_z0":float(c[j]),
                    "z_max_cs2":float(z[mask][np.argmax(c[mask])]),
                    "subluminal_z100":bool(np.min(c[mask])>=0 and np.max(c[mask])<=1),
                    "sample_count_z100":int(np.sum(mask)),
                    "sample_count_all":int(np.sum(full)),
                })
            except Exception as exc:
                rec["parse_error"]=str(exc)
        else:
            rec["CLASS_tail"]=p.stdout[-650:]
        records.append(rec)
        print("D0_CORRECTION_POINT",json.dumps(rec,sort_keys=True),flush=True)

    if D0_ACCEPTED not in curves:
        raise RuntimeError("The accepted benchmark did not produce a background curve")
    Nr,Nsr=curves[D0_ACCEPTED]
    zref=np.expm1(-Nr)
    Ss=S_for_z(zref)
    mask=(Ss>1e-4)&(Ss<1-1e-4)&np.isfinite(Nsr)
    required=(Nsr[mask]-DFLOOR)/Ss[mask]
    necessary=max(0.0,float(np.max(required)))
    for rec in records:
        d0=rec["D0"]
        if d0 not in curves:continue
        N,Ns=curves[d0]
        inter=np.interp(Nr,N,Ns)
        scale=np.maximum(1e-7,np.abs(Nsr))
        rec["sound_numerator_rel_difference_max_vs_accepted"]=float(np.max(np.abs(inter-Nsr)/scale))
    summary={
        "status":"hi_class numerical D0 correction causal bound audit",
        "theory_equation":"D=Dfloor+D0*S; c_s^2=N_s/D. For N_s>0 and S>0, subluminality requires D0 >= max_N ((N_s-Dfloor)/S), if N_s does not depend on D0.",
        "AF":AF,"zc":ZC,"width":WIDTH,"Dfloor":DFLOOR,"z_t":ZT,
        "accepted_D0":D0_ACCEPTED,"records":records,
        "necessary_lower_D0_from_reference_numerator_and_subluminality_z100":necessary,
        "warning":"This is a conditional bound inferred from the already-selected background and Planck-mass trajectory. A bound is not a unique kinetic-correction prediction.",
        "independence_caveat":"The sound-speed numerator must be demonstrated invariant, not simply assumed; see recorded differences across actual CLASS cases."
    }
    (OUT/"causal_scan.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    max_numerator_change=max((q.get("sound_numerator_rel_difference_max_vs_accepted",0.0) for q in records),default=0.0)
    with_cs=sorted([q for q in records if "max_cs2_z100" in q],key=lambda q:q["D0"])
    bracket=next(( (a,b) for a,b in zip(with_cs,with_cs[1:]) if a["max_cs2_z100"]>1.0>=b["max_cs2_z100"]),None)
    numerical_bound=None
    if bracket:
        a,b=bracket
        numerical_bound=a["D0"]+(1-a["max_cs2_z100"])*(b["D0"]-a["D0"])/(b["max_cs2_z100"]-a["max_cs2_z100"])
    summary["numerical_bracket"]=[a["D0"],b["D0"]] if bracket else None
    summary["numerical_interpolated_cs2_eq1"] = numerical_bound
    summary["max_relative_sound_numerator_change_z100_vs_reference"]=max_numerator_change
    (OUT/"causal_scan.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\\n".replace("\\\\n","\\n"))
    print("D0_CAUSALITY_SCAN_SUMMARY",json.dumps({
        "n":len(records),
        "successful":sum(q.get("returncode")==0 for q in records),
        "superluminal":sum(not q.get("subluminal_z100",False) for q in records if "subluminal_z100" in q),
        "necessary_lower_D0":necessary,
        "causal_D0_values":[q["D0"] for q in records if q.get("subluminal_z100")],
        "causal_D0_values_all_z":[q["D0"] for q in records if q.get("subluminal_all")],
        "numerical_threshold_bracket":summary["numerical_bracket"],
        "numerical_interpolated_threshold":numerical_bound,
        "max_numerator_rel_difference_across_scan":max_numerator_change,
    },sort_keys=True))

if __name__=="__main__":
    main()
