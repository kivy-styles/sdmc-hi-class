#!/usr/bin/env python3
"""
Run a background-level basin-of-attraction test for the field-redefined
SDMC structural-clock action.

Chronological closure in the log structural field is

    psi = ln(S/S0) = ln(t/t0),

hence d psi / d ln t = 1 and

    (psi-psi0) - ln(t/t0) = 0.

We perturb the initial scalar position and velocity while leaving the covariant
functions G2,G3,G4 unchanged.  The potential offset (parameters_smg[0]) remains
the hi_class shooting/tuning parameter.  The audit asks whether small
perturbations decay back toward the chronological trajectory.
"""
from pathlib import Path
import json, re, subprocess
import numpy as np
from scipy.interpolate import CubicSpline

CASES = [
    ("nominal", 0.0, 0.0),
    ("psi_p1e6", +1e-6, 0.0),
    ("psi_m1e6", -1e-6, 0.0),
    ("vel_p1e6", 0.0, +1e-6),
    ("vel_m1e6", 0.0, -1e-6),
    ("psi_p1e4", +1e-4, 0.0),
    ("psi_m1e4", -1e-4, 0.0),
    ("vel_p1e4", 0.0, +1e-4),
    ("vel_m1e4", 0.0, -1e-4),
    ("both_p1e3", +1e-3, +1e-3),
    ("both_m1e3", -1e-3, -1e-3),
    ("cross_pm1e3", +1e-3, -1e-3),
    ("cross_mp1e3", -1e-3, +1e-3),
    ("psi_p1e2", +1e-2, 0.0),
    ("psi_m1e2", -1e-2, 0.0),
    ("vel_p1e2", 0.0, +1e-2),
    ("vel_m1e2", 0.0, -1e-2),
]

COMMON = """H0 = 69.45160505326464
omega_b = 0.02208511574370519
omega_cdm = 0.12298509428428875
N_ncdm = 0
N_ur = 3.046
T_cmb = 2.7255
YHe = 0.2453
A_s = 2.1197359302086e-9
n_s = 0.962555355281866
tau_reio = 0.055202901571989066
Omega_Lambda = 0
Omega_fld = 3.1443554e-8
fluid_equation_of_state = SDMC_TRACKER
cs2_fld = 0.003
use_ppf = no
Omega_smg = -1
gravity_model = sdmc_v3_covariant_logstruct_audit
pert_initial_conditions_smg = zero
method_qs_smg = fully_dynamic
a_ini_over_a_today_default = 2.e-8
a_ini_test_qs_smg = 2.e-8
pert_ic_ini_z_ref_smg = 5.e7
a_min_stability_test_smg = 2.e-8
output_background_smg = 3
modes = s
output = mPk
P_k_max_h/Mpc = 0.2
z_pk = 0
write background = yes
format = class
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
output_verbose = 0
"""

def read_bg(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"(?:^|\s)1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    names=[]
    for k,m in enumerate(ms):
        e=ms[k+1].start() if k+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

def scalar_col(d):
    for k in d:
        kl=k.lower()
        if kl.startswith("phi") and "prime" not in kl:
            return k
    raise RuntimeError(f"could not locate scalar field column: {list(d)}")

def prepare(d):
    z=np.asarray(d["z"])
    N=-np.log1p(z)
    t=np.asarray(d["proper time [Gyr]"])
    H=np.asarray(d["H [1/Mpc]"])
    pk=scalar_col(d)
    psi=np.asarray(d[pk])
    o=np.argsort(N)
    N,z,t,H,psi=[x[o] for x in (N,z,t,H,psi)]
    keep=np.r_[True,np.diff(N)>1e-12]
    N,z,t,H,psi=[x[keep] for x in (N,z,t,H,psi)]
    lnt=np.log(t)
    # d psi / d ln t; chronological solution requires q=1.
    q=CubicSpline(lnt,psi)(lnt,1)
    i0=int(np.argmin(np.abs(z)))
    chrono=(psi-psi[i0])-np.log(t/t[i0])
    # u = d(e^psi)/dt = e^psi * q / t. Normalize to present.
    urel=np.exp(psi-psi[i0])*(q/q[i0])*(t[i0]/t)
    return dict(N=N,z=z,t=t,H=H,psi=psi,q=q,chrono=chrono,urel=urel,i0=i0,phi_name=pk)

def interp(src, N):
    o=np.argsort(src["N"])
    out={}
    for k in ["H","psi","q","chrono","urel","z","t"]:
        out[k]=np.interp(N,src["N"][o],src[k][o])
    return out

outdir=Path("output/structural_clock_attractor")
outdir.mkdir(parents=True,exist_ok=True)
runs={}
prepared={}

for name,dpsi,dv in CASES:
    ini=outdir/f"{name}.ini"
    root=f"output/structural_clock_attractor/{name}_"
    ini.write_text(COMMON + f"parameters_smg = 0.0, {dpsi:.17e}, {dv:.17e}\nroot = {root}\n")
    cp=subprocess.run(["./class",str(ini)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (outdir/f"{name}.log").write_text(cp.stdout)
    rec={"delta_psi_ini":dpsi,"delta_v_ini":dv,"returncode":cp.returncode}
    bg=Path(root+"00_background.dat")
    if cp.returncode!=0 or not bg.exists():
        rec["status"]="failed"
        rec["log_tail"]="\n".join(cp.stdout.splitlines()[-30:])
        runs[name]=rec
        continue
    d=prepare(read_bg(bg))
    prepared[name]=d
    i0=d["i0"]
    m100=d["z"]<=100
    m10=d["z"]<=10
    rec.update({
      "status":"ok",
      "phi_column":d["phi_name"],
      "psi0":float(d["psi"][i0]),
      "q0":float(d["q"][i0]),
      "max_abs_q_minus_1_z100":float(np.max(np.abs(d["q"][m100]-1.))),
      "max_abs_chrono_resid_z100":float(np.max(np.abs(d["chrono"][m100]))),
      "max_abs_urel_minus_1_z100":float(np.max(np.abs(d["urel"][m100]-1.))),
      "max_abs_q_minus_1_z10":float(np.max(np.abs(d["q"][m10]-1.))),
      "max_abs_chrono_resid_z10":float(np.max(np.abs(d["chrono"][m10]))),
      "max_abs_urel_minus_1_z10":float(np.max(np.abs(d["urel"][m10]-1.))),
      "age0_Gyr":float(d["t"][i0]),
      "H0_internal_1Mpc":float(d["H"][i0]),
    })
    runs[name]=rec

if "nominal" not in prepared:
    raise SystemExit("nominal structural-clock run failed; cannot perform attractor comparison")

nom=prepared["nominal"]
sample_z=[1e6,1e5,1e4,1e3,100,30,10,3,1,.3,0]
for name,d in prepared.items():
    if name=="nominal":
        continue
    Nlo=max(np.min(nom["N"]),np.min(d["N"]))
    Nhi=min(np.max(nom["N"]),np.max(d["N"]))
    m=(d["N"]>=Nlo)&(d["N"]<=Nhi)
    N=d["N"][m]
    ni=interp(nom,N)
    dphi=d["psi"][m]-ni["psi"]
    dq=d["q"][m]-ni["q"]
    dH=d["H"][m]/ni["H"]-1.
    z=d["z"][m]
    i0=int(np.argmin(np.abs(z)))
    early=int(np.argmax(z))
    amp_ini=float(np.hypot(dphi[early],dq[early]))
    amp_0=float(np.hypot(dphi[i0],dq[i0]))
    rec=runs[name]
    rec.update({
      "delta_psi_vs_nominal_present":float(dphi[i0]),
      "delta_q_vs_nominal_present":float(dq[i0]),
      "deltaH_H_vs_nominal_present":float(dH[i0]),
      "max_abs_deltaH_H_z100":float(np.max(np.abs(dH[z<=100]))),
      "phase_space_amp_earliest_common":amp_ini,
      "phase_space_amp_present":amp_0,
      "phase_space_decay_ratio":float(amp_0/amp_ini) if amp_ini>0 else None,
      "samples_vs_nominal":[]
    })
    for zt in sample_z:
        jj=int(np.argmin(np.abs(z-zt)))
        rec["samples_vs_nominal"].append({
          "z":float(z[jj]),
          "delta_psi":float(dphi[jj]),
          "delta_q":float(dq[jj]),
          "deltaH_H":float(dH[jj]),
        })

summary={
  "purpose":"Test whether the chronological structural-clock solution is a local dynamical attractor under scalar IC perturbations; this does not make the reconstructed action a unique microscopic derivation.",
  "criterion":{
    "chronological_field_relation":"(psi-psi0)-ln(t/t0)=0",
    "chronological_slope":"q=dpsi/dln(t)=1",
    "constant_structural_speed":"u=d(exp(psi))/dt=constant"
  },
  "cases":runs
}
(outdir/"structural_clock_attractor_audit.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")

print("STRUCTURAL_CLOCK_ATTRACTOR_AUDIT")
for name in [x[0] for x in CASES]:
    print(name,json.dumps(runs[name],sort_keys=True))
