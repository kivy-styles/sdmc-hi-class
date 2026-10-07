#!/usr/bin/env python3
from pathlib import Path
import json,re
import numpy as np
from scipy.interpolate import CubicSpline

BG=Path("input_action/linear_cov_free_00_background.dat")
if not BG.exists():
    raise RuntimeError(f"missing accepted-action free background: {BG}")

lines=BG.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr))
names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
arr=np.loadtxt(BG)
d={n:arr[:,i] for i,n in enumerate(names)}

def col(name):
    if name not in d:
        raise RuntimeError(f"missing column {name}; available={names}")
    return np.asarray(d[name],float)

z=col("z")
N=-np.log1p(z)
order=np.argsort(N)
z,N=z[order],N[order]
H=col("H [1/Mpc]")[order]
F=col("M*^2_smg")[order]
aM=col("M2_running_smg")[order]
aB=col("braiding_smg")[order]
D=col("kin (D)")[order]
cs2=col("c_s^2")[order]
rho=col("(.)rho_tot")[order]
pre=col("(.)p_tot")[order]

# Geometric kinematics from the freely evolved accepted action.
lnH=CubicSpline(N,np.log(H))
dlnH=lnH(N,1)
q=-1.0-dlnH
wtot=pre/rho
noslip=aB+2.0*aM

# On the accepted trajectory, exact No-Slip reconstruction gives
# X G3_X = -F_phi/2; with phi=N and X=H^2/2, hence the dimensionless
# trajectory braiding contribution H^2 G3_X = -F_phi = -alpha_M F.
gH2=-aM*F

targets={
  "F":1.02388542769,
  "alphaM":0.0,
  "alphaB":0.0,
  "D":2.0,
  "cs2":1.0,
  "q":0.0,
  "wtot":-1.0/3.0,
  "gH2":0.0,
  "noslip":0.0
}
series={"F":F,"alphaM":aM,"alphaB":aB,"D":D,"cs2":cs2,
        "q":q,"wtot":wtot,"gH2":gH2,"noslip":noslip}

sample_z=[5.0,3.0,2.0,1.0,0.5,0.2,0.1,0.0]
rows=[]
for zs in sample_z:
    i=int(np.argmin(np.abs(z-zs)))
    row={"z_requested":zs,"z":float(z[i]),"N":float(N[i])}
    for k,v in series.items():
        row[k]=float(v[i])
        row[k+"_abs_distance_to_mature"]=float(abs(v[i]-targets[k]))
    rows.append(row)

# Endpoint trend from a local cubic fit over the last 8% of the supported N range,
# capped to z<0.5. This is a present-side derivative only, NOT a future extrapolation.
mask=(z<=0.5)
inds=np.where(mask)[0]
if len(inds)<8:
    inds=np.arange(max(0,len(N)-max(8,len(N)//12)),len(N))
else:
    inds=inds[max(0,len(inds)-max(12,len(inds)//3)):]
endpoint={}
for k,v in series.items():
    sp=CubicSpline(N[inds],v[inds])
    slope=float(sp(0.0,1))
    val=float(sp(0.0))
    target=targets[k]
    # Negative derivative of absolute distance means locally approaching target.
    delta=val-target
    distance_derivative=(1.0 if delta>=0 else -1.0)*slope
    endpoint[k]={
      "today":val,
      "mature_target":target,
      "d_dN_today":slope,
      "d_abs_distance_dN_today":float(distance_derivative),
      "local_direction":"toward" if distance_derivative<0 else ("away" if distance_derivative>0 else "flat"),
      "fractional_distance_today":float(abs(delta)/(abs(target)+1e-30)) if target!=0 else None
    }

# Broad supported-domain comparison: distance at z~1 vs today.
iz1=int(np.argmin(np.abs(z-1.0)))
i0=int(np.argmin(np.abs(z)))
progress={}
for k,v in series.items():
    e1=float(abs(v[iz1]-targets[k]))
    e0=float(abs(v[i0]-targets[k]))
    progress[k]={
      "abs_distance_z1":e1,
      "abs_distance_today":e0,
      "distance_ratio_today_over_z1":float(e0/e1) if e1>0 else None,
      "z1_to_today":"closer" if e0<e1 else ("farther" if e0>e1 else "unchanged")
    }

# Mature canonical closure is multi-dimensional. Count only qualitative supported-domain
# tendencies; do not infer future completion from them.
toward=[k for k,v in endpoint.items() if v["local_direction"]=="toward"]
away=[k for k,v in endpoint.items() if v["local_direction"]=="away"]
flat=[k for k,v in endpoint.items() if v["local_direction"]=="flat"]

out={
  "status":"accepted-action supported-domain canonical-approach audit",
  "scope":{
    "action":"unchanged accepted frozen linear-G3 covariant action",
    "domain":"freely evolved reconstructed action over its supported past-to-present table only (N<=0)",
    "future_extrapolation_used":False,
    "reason":"hi_class integrates to log(a/a0)=0 and the reconstructed action table ends at the present endpoint; no arbitrary future spline extrapolation is treated as physics."
  },
  "manuscript_mature_targets":targets,
  "samples":rows,
  "z1_to_today_progress":progress,
  "present_side_endpoint_trends":endpoint,
  "trend_summary":{"toward":toward,"away":away,"flat":flat},
  "interpretation":{
    "positive_result":"Quantities whose supported-domain distance to the mature target is shrinking are already evolving consistently with the manuscript's canonical endpoint before today.",
    "negative_or_flat_result":"A flat or wrong-way present-side trend is evidence that the currently reconstructed late action does not itself exhibit the full canonical transverse flow within its supported domain.",
    "important_limit":"Neither outcome derives the future microscopic transverse flow; the manuscript explicitly leaves that mechanism open."
  }
}
Path("output").mkdir(exist_ok=True)
Path("output/accepted_action_canonical_approach.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("ACCEPTED_ACTION_CANONICAL_APPROACH",json.dumps(out,sort_keys=True))
