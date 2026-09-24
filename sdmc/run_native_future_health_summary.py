#!/usr/bin/env python3
"""
Summarize native future hi_class background runs and compare them with the
independent exact homogeneous structural-action audit.
"""
from pathlib import Path
import json, math, re
import numpy as np
from scipy.interpolate import CubicSpline

OUT=Path("output/native_future_health_summary.json")
REF=Path("output/future_homogeneous_covariant_stitch_audit.json")

TP=5.391247e-44
S0=2.72e61
C_MS=299792458.
MPC_M=3.085677581491367e22
TP_MPC=TP*C_MS/MPC_M

RUNS={
  "matched_chi_F":Path("output/native_future_matched_00_background.dat"),
  "bounded_chi":Path("output/native_future_bounded_00_background.dat"),
}

def read(path):
    lines=path.read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and
         re.search(r"(?:^|\s)1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

ref=json.loads(REF.read_text()) if REF.exists() else {"candidates":{}}
out={"status":"native future hi_class background/stability comparison","runs":{}}

for name,path in RUNS.items():
    if not path.exists():
        out["runs"][name]={"status":"missing","path":str(path)}
        continue
    d=read(path)
    required=["z","H [1/Mpc]","phi_smg","phi'","c_s^2","kin (D)",
              "braiding_smg","M2_running_smg","M*^2_smg"]
    missing=[k for k in required if k not in d]
    if missing:
        out["runs"][name]={"status":"missing_columns","missing":missing}
        continue

    z=np.asarray(d["z"])
    a=1./(1.+z)
    N=np.log(a)
    o=np.argsort(N)
    N=N[o]; a=a[o]
    keep=np.r_[True,np.diff(N)>1e-12]
    N=N[keep]; a=a[keep]
    def A(k): return np.asarray(d[k])[o][keep]

    H=A("H [1/Mpc]")
    psi=A("phi_smg")
    phip=A("phi'")
    F=A("M*^2_smg")
    D=A("kin (D)")
    cs=A("c_s^2")
    bra=A("braiding_smg")
    run=A("M2_running_smg")
    sig=np.exp(psi)
    dotpsi=phip/a
    pstruct=dotpsi/H
    Nstruct=TP_MPC*S0*sig*dotpsi
    q=-1.-CubicSpline(N,np.log(H))(N,1)
    ns=bra+2.*run

    future=N>=-1e-8
    if not np.any(future):
        out["runs"][name]={"status":"no_future_rows"}
        continue
    jf=np.where(future)[0]

    j0=int(np.argmin(np.abs(N)))
    jend=jf[-1]
    jD=jf[np.argmin(D[jf])]
    jcs=jf[np.argmin(cs[jf])]
    jcsmax=jf[np.argmax(cs[jf])]
    jns=jf[np.argmax(np.abs(ns[jf]))]

    # q=0 crossings in the future.
    crossings=[]
    for ia,ib in zip(jf[:-1],jf[1:]):
        if q[ia]==0. or q[ia]*q[ib]<0.:
            frac=-q[ia]/(q[ib]-q[ia])
            nn=N[ia]+frac*(N[ib]-N[ia])
            crossings.append({"ln_a":float(nn),"a":float(math.exp(nn))})

    comparisons=[]
    rcan=ref.get("candidates",{}).get(name,{})
    for rr in rcan.get("samples",[]):
        nt=float(rr.get("ln_a",0.))
        if nt < N[jf[0]]-1e-9 or nt > N[jend]+1e-9:
            continue
        def interp(y):
            return float(CubicSpline(N,y)(nt))
        row={"ln_a":nt}
        pairs=[
          ("H",H,"H"),
          ("p",pstruct,"p"),
          ("q",q,"q"),
          ("N_struct",Nstruct,"N_struct"),
          ("D_horndeski",D,"D_horndeski"),
          ("cs2_horndeski",cs,"cs2_horndeski"),
        ]
        for lab,y,keyref in pairs:
            if keyref not in rr:
                continue
            nat=interp(y); rv=float(rr[keyref])
            row[lab+"_native"]=nat
            row[lab+"_reference"]=rv
            row[lab+"_abs_diff"]=nat-rv
            row[lab+"_rel_diff"]=(nat-rv)/(abs(rv)+1e-300)
        comparisons.append(row)

    out["runs"][name]={
      "status":"ok",
      "path":str(path),
      "future_extent":{"ln_a_max":float(N[jend]),"a_max":float(a[jend])},
      "present":{
        "ln_a":float(N[j0]),"H":float(H[j0]),"psi":float(psi[j0]),
        "sigma":float(sig[j0]),"N_struct":float(Nstruct[j0]),
        "p":float(pstruct[j0]),"q":float(q[j0]),"F":float(F[j0]),
        "D":float(D[j0]),"cs2":float(cs[j0]),
        "alphaB_plus_2alphaM":float(ns[j0]),
      },
      "future_health":{
        "min_D":float(D[jD]),"min_D_ln_a":float(N[jD]),
        "min_cs2":float(cs[jcs]),"min_cs2_ln_a":float(N[jcs]),
        "max_cs2":float(cs[jcsmax]),"max_cs2_ln_a":float(N[jcsmax]),
        "max_abs_alphaB_plus_2alphaM":float(abs(ns[jns])),
        "max_abs_alphaB_plus_2alphaM_ln_a":float(N[jns]),
      },
      "future_kinematics":{
        "q_zero_crossings":crossings,
        "q_max":float(np.max(q[jf])),
        "q_max_ln_a":float(N[jf[np.argmax(q[jf])]]),
      },
      "final":{
        "ln_a":float(N[jend]),"a":float(a[jend]),"sigma":float(sig[jend]),
        "N_struct":float(Nstruct[jend]),"p":float(pstruct[jend]),
        "q":float(q[jend]),"F":float(F[jend]),
        "D":float(D[jend]),"cs2":float(cs[jend]),
        "alphaB_plus_2alphaM":float(ns[jend]),
      },
      "reference_comparison":comparisons,
    }

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NATIVE_FUTURE_HEALTH_SUMMARY")
for name,r in out["runs"].items():
    print(name,json.dumps(r if r.get("status")!="ok" else {
      "future_extent":r["future_extent"],
      "present":r["present"],
      "health":r["future_health"],
      "kinematics":r["future_kinematics"],
      "final":r["final"],
    },sort_keys=True))
