#!/usr/bin/env python3
from pathlib import Path
import json,re
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

BG=Path("input_action/linear_cov_free_00_background.dat")
lines=BG.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(BG); d={n:a[:,i] for i,n in enumerate(names)}
z=np.asarray(d["z"]); N=-np.log1p(z); o=np.argsort(N); z=z[o]; N=N[o]
def v(name): return np.asarray(d[name])[o]
H=v("H [1/Mpc]"); F=v("M*^2_smg"); aM=v("M2_running_smg"); aB=v("braiding_smg")
D=v("kin (D)"); cs2=v("c_s^2"); rho=v("(.)rho_tot"); p=v("(.)p_tot")
wt=p/rho; q=-1-CubicSpline(N,np.log(H))(N,1); gH2=-aM*F
series={"F":F,"alphaM":aM,"alphaB":aB,"D":D,"cs2":cs2,"wtot":wt,"q":q,"gH2":gH2}
sp={k:CubicSpline(N,x) for k,x in series.items()}

def nz(Nx): return float(np.exp(-Nx)-1)
def state(Nx):
    return {"N":float(Nx),"a":float(np.exp(Nx)),"z":nz(Nx),
            **{k:float(s(Nx)) for k,s in sp.items()}}

def latest_root(key,target):
    y=sp[key](N)-target; rr=[]
    for i in range(len(N)-1):
        if y[i]*y[i+1]<0:
            rr.append(brentq(lambda x: float(sp[key](x)-target),N[i],N[i+1]))
    return max(rr) if rr else None

events=[]

# Peak Planck-mass running magnitude: marks fastest F transition.
res=minimize_scalar(lambda x:-abs(float(sp["alphaM"](x))),bounds=(max(N.min(),-3.0),0),method="bounded")
events.append({"name":"peak_abs_alphaM","state":state(float(res.x))})

# Minimum scalar sound speed in late universe z<=10.
lo=max(N.min(),-np.log(11.0))
res=minimize_scalar(lambda x:float(sp["cs2"](x)),bounds=(lo,0),method="bounded")
events.append({"name":"minimum_cs2_zle10","state":state(float(res.x))})

# D halfway from its earliest z<=10 value to present is a useful kinetic-response marker.
Nlo=lo; Dlo=float(sp["D"](Nlo)); D0=float(sp["D"](0)); Dmid=.5*(Dlo+D0)
try:
    Nm=brentq(lambda x:float(sp["D"](x)-Dmid),Nlo,0)
    events.append({"name":"D_half_late_response","target_D":Dmid,"state":state(Nm)})
except ValueError:
    pass

# Running/braid threshold crossings on their final decay branch.
for key,thr in [("alphaM",1e-2),("alphaM",1e-3),("gH2",1e-2),("gH2",1e-3)]:
    vals=np.abs(series[key])-thr; rr=[]
    for i in range(len(N)-1):
        if vals[i]*vals[i+1]<0:
            rr.append(brentq(lambda x:abs(float(sp[key](x)))-thr,N[i],N[i+1]))
    if rr:
        events.append({"name":f"{key}_abs_below_{thr:g}_final_crossing","state":state(max(rr))})

Nq=latest_root("q",0.0)
if Nq is not None:
    events.append({"name":"kinematic_coasting_q0","state":state(Nq)})

# Present supported endpoint.
events.append({"name":"today_supported_endpoint","state":state(0.0)})

events=sorted(events,key=lambda e:e["state"]["N"])

targets={"F":1.02388542769,"D":2.0,"cs2":1.0,"alphaM":0.0,"alphaB":0.0,"gH2":0.0,"q":0.0,"wtot":-1/3}
for e in events:
    e["distance_to_mature"]={k:float(abs(e["state"][k]-t)) for k,t in targets.items()}

# Determine chronological ordering and main residual at today.
today=state(0.0)
residuals={
 "F_fractional":abs(today["F"]-targets["F"])/targets["F"],
 "D_fractional":abs(today["D"]-2)/2,
 "cs2_absolute":abs(today["cs2"]-1),
 "alphaM_absolute":abs(today["alphaM"]),
 "alphaB_absolute":abs(today["alphaB"]),
 "gH2_absolute":abs(today["gH2"]),
 "q_absolute":abs(today["q"]),
 "wtot_absolute":abs(today["wtot"]+1/3)
}
out={
 "status":"accepted-action finite-handoff chronology audit",
 "scope":"unchanged accepted frozen action, free evolution over supported N<=0 domain",
 "events":events,
 "today_residuals_to_mature":residuals,
 "chronology":[e["name"] for e in events],
 "interpretation":{
   "sector_ordering":"The audit orders the Planck-mass/braiding transition, kinetic sound-speed response, total-background coasting crossing, and present endpoint without extrapolating beyond today.",
   "canonicalization_test":"If coasting follows the F/G3 transition while D and c_s^2 remain far from their mature values, the finite handoff is demonstrably multi-sector and non-one-dimensional.",
   "overshoot_test":"If q and w_tot cross their mature coasting values and continue to more negative present values, the supported late-time background contains a kinematic overshoot rather than monotonic approach to the mature fixed point.",
   "future_limit":"The mechanism that later reverses that overshoot and completes canonicalization is not derived by this audit."
 }
}
Path("output").mkdir(exist_ok=True)
Path("output/accepted_action_handoff_chronology.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("ACCEPTED_ACTION_HANDOFF_CHRONOLOGY",json.dumps(out,sort_keys=True))
