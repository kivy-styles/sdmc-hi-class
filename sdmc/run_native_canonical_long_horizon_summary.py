#!/usr/bin/env python3
"""
Long-horizon native closure summary for the legacy-canonical SDMC endpoint.

This is the stricter follow-up to the ln(a)=5 native gate.  It requires the
native hi_class background to reach ln(a)=10 and to approach the analytically
closed r=0, lambda_l=sqrt(2) endpoint,

    N_inf = sqrt(8*pi/3), p_inf=1, q_inf=0, D_inf=2, c_s^2_inf=1.

For native perturbations, raw synchronous scalar amplitudes are reported but
are not used as a hard physical stability criterion.  The additional hard
late-time check is that both metric potentials remain finite and do not develop
a late exponential runaway between ln(a)=5 and ln(a)=10.
"""
from pathlib import Path
import json,math,re
import numpy as np
from scipy.interpolate import CubicSpline

OUT=Path("output/native_canonical_long_horizon_summary.json")
BG=Path("output/native_future_canonical_long_00_background.dat")
BGRC=Path("output/native_future_canonical_long.rc")
BLOG=Path("output/native_future_canonical_long.log")
PPREFIX="native_future_perturb_canonical_long_"
PRC=Path("output/native_future_perturb_canonical_long.rc")
PLOG=Path("output/native_future_perturb_canonical_long.log")

TP=5.391247e-44
S0=2.72e61
C_MS=299792458.
MPC_M=3.085677581491367e22
TP_MPC=TP*C_MS/MPC_M
XI=math.sqrt(8.*math.pi/3.)

def read_table(path):
    lines=path.read_text(errors="replace").splitlines()
    hdr=None
    kval=None
    for line in lines[:8]:
        m=re.search(r"mode k\s*=\s*([0-9eE+\-.]+)",line)
        if m: kval=float(m.group(1))
    for line in lines:
        if line.startswith("#") and re.search(r"(?:^|\s)1\s*:",line):
            hdr=line.lstrip("#").strip()
    if hdr is None:
        raise RuntimeError(f"no numbered header in {path}")
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    if a.ndim==1: a=a[None,:]
    if a.shape[1]!=len(names):
        raise RuntimeError(f"{path}: {a.shape[1]} cols vs {len(names)} titles")
    return kval,names,a

def rc(path):
    if not path.exists(): return None
    try: return int(path.read_text().strip())
    except Exception: return None

def col(names,a,*opts):
    for q in opts:
        if q in names:
            return np.asarray(a[:,names.index(q)],float),q
    return None,None

out={
  "status":"native canonical long-horizon closure audit",
  "analytic_target":{
    "ln_a":10.0,"N_inf":XI,"p_inf":1.0,"q_inf":0.0,
    "D_inf":2.0,"cs2_inf":1.0
  }
}

# ---------------- Background ----------------
bg={"returncode":rc(BGRC),"path":str(BG)}
if BLOG.exists(): bg["log_tail"]=BLOG.read_text(errors="replace").splitlines()[-30:]
if not BG.exists():
    bg["status"]="missing"
    bg["gate_pass"]=False
else:
    _,names,a=read_table(BG)
    required=["z","H [1/Mpc]","phi_smg","phi'","c_s^2","kin (D)",
              "braiding_smg","M2_running_smg","M*^2_smg"]
    miss=[x for x in required if x not in names]
    if miss:
        bg.update(status="missing_columns",missing=miss,gate_pass=False)
    else:
        def C(name): return np.asarray(a[:,names.index(name)],float)
        z=C("z"); aa=1./(1.+z); N=np.log(aa)
        order=np.argsort(N)
        N=N[order]; aa=aa[order]
        keep=np.r_[True,np.diff(N)>1e-12]
        N=N[keep]; aa=aa[keep]
        def A(name): return C(name)[order][keep]
        H=A("H [1/Mpc]"); psi=A("phi_smg"); phip=A("phi'")
        D=A("kin (D)"); cs=A("c_s^2"); F=A("M*^2_smg")
        ns=A("braiding_smg")+2.*A("M2_running_smg")
        sig=np.exp(psi)
        dotpsi=phip/aa
        p=dotpsi/H
        Nstruct=TP_MPC*S0*sig*dotpsi
        q=-1.-CubicSpline(N,np.log(H))(N,1)
        fut=N>=0.
        jf=np.where(fut)[0]; j=jf[-1]
        finite=bool(np.all(np.isfinite(np.c_[H,psi,phip,D,cs,F,ns,p,Nstruct,q])))
        final={
          "ln_a":float(N[j]),"a":float(aa[j]),"N_struct":float(Nstruct[j]),
          "p":float(p[j]),"q":float(q[j]),"D":float(D[j]),
          "cs2":float(cs[j]),"F":float(F[j]),
          "alphaB_plus_2alphaM":float(ns[j]),"sigma":float(sig[j])
        }
        errs={
          "N_rel":float(Nstruct[j]/XI-1.),
          "p_minus_1":float(p[j]-1.),
          "q":float(q[j]),
          "D_minus_2":float(D[j]-2.),
          "cs2_minus_1":float(cs[j]-1.)
        }
        health={
          "min_D":float(np.min(D[jf])),
          "min_cs2":float(np.min(cs[jf])),
          "max_cs2":float(np.max(cs[jf])),
          "max_abs_alphaB_plus_2alphaM":float(np.max(np.abs(ns[jf])))
        }
        target=(
          abs(errs["N_rel"])<=2e-4 and
          abs(errs["p_minus_1"])<=2e-4 and
          abs(errs["q"])<=2e-4 and
          abs(errs["D_minus_2"])<=1e-2 and
          abs(errs["cs2_minus_1"])<=1e-4
        )
        gate=bool(bg["returncode"]==0 and finite and N[j]>=9.99 and
                  health["min_D"]>0. and health["min_cs2"]>0. and
                  health["max_cs2"]<=1.00001 and target)
        bg.update(status="ok" if gate else "target_or_health_miss",
                  finite=finite,final=final,errors=errs,health=health,
                  target_gate_pass=target,gate_pass=gate)
out["background"]=bg

# ---------------- Perturbations ----------------
pr={"returncode":rc(PRC),"files":[]}
if PLOG.exists(): pr["log_tail"]=PLOG.read_text(errors="replace").splitlines()[-30:]
files=sorted(Path("output").glob(PPREFIX+"*perturbations_k*_s.dat"))
allp=(pr["returncode"]==0 and bool(files))
for path in files:
    try:
        kval,names,a=read_table(path)
        av,_=col(names,a,"a")
        if av is None: raise RuntimeError("missing a")
        Nv=np.log(av)
        finite=bool(np.all(np.isfinite(a)))
        row={"path":str(path),"k_Mpc_inv":kval,
             "a_max":float(np.max(av)),"ln_a_max":float(np.max(Nv)),
             "all_columns_finite":finite}
        metric_gate=True
        for key,opt in [("psi",("psi",)),("phi",("phi",))]:
            y,used=col(names,a,*opt)
            if y is None:
                row[key]={"status":"missing"}
                metric_gate=False
                continue
            late=Nv>=5.
            if np.count_nonzero(late)<3:
                row[key]={"status":"insufficient_late_rows"}
                metric_gate=False
                continue
            yl=np.abs(y[late])
            y0=max(float(yl[0]),1e-300)
            growth=float(np.max(yl)/y0)
            final_ratio=float(yl[-1]/y0)
            # A generous no-runaway condition.  Canonical super-Hubble metric
            # modes should freeze or decay; this only rejects a resolved late
            # explosive mode, not small ringing.
            ok=bool(np.all(np.isfinite(yl)) and growth<=2.0)
            row[key]={
              "column":used,"late_start_ln_a":float(Nv[late][0]),
              "late_growth_factor":growth,
              "late_final_over_start":final_ratio,
              "late_gate_pass":ok
            }
            metric_gate &= ok
        scalar,used=col(names,a,"delta_phi_smg (sync)","V_x_smg (sync)")
        if scalar is not None:
            late=Nv>=5.
            scale=max(abs(float(scalar[late][0])),1e-300) if np.any(late) else 1.
            row["synchronous_scalar"]={
              "column":used,
              "late_max_abs_over_start":
                float(np.max(np.abs(scalar[late]))/scale) if np.any(late) else None,
              "late_final_abs_over_start":
                float(abs(scalar[late][-1])/scale) if np.any(late) else None,
              "note":"reported only; raw synchronous scalar amplitude is gauge-dependent"
            }
        row["gate_pass"]=bool(finite and row["ln_a_max"]>=9.99 and metric_gate)
        allp &= row["gate_pass"]
        pr["files"].append(row)
    except Exception as e:
        pr["files"].append({"path":str(path),"status":"parse_error","error":str(e),"gate_pass":False})
        allp=False
pr["status"]="ok" if allp else "incomplete_or_failed"
pr["gate_pass"]=bool(allp)
out["perturbations"]=pr
out["all_long_horizon_gates_pass"]=bool(bg.get("gate_pass") and pr.get("gate_pass"))

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("NATIVE_CANONICAL_LONG_HORIZON_SUMMARY")
print("BACKGROUND",json.dumps({
  "status":bg.get("status"),"gate":bg.get("gate_pass"),
  "final":bg.get("final"),"errors":bg.get("errors"),"health":bg.get("health")
},sort_keys=True))
for x in pr["files"]:
    print("PERTURB",json.dumps(x,sort_keys=True))
print("ALL_LONG_HORIZON_GATES_PASS",out["all_long_horizon_gates_pass"])
