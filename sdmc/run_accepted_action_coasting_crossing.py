#!/usr/bin/env python3
from pathlib import Path
import json,re
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

BG=Path("input_action/linear_cov_free_00_background.dat")
lines=BG.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(BG)
d={n:a[:,i] for i,n in enumerate(names)}
z=np.asarray(d["z"]); N=-np.log1p(z); o=np.argsort(N); z=z[o]; N=N[o]
def v(name): return np.asarray(d[name])[o]
H=v("H [1/Mpc]"); F=v("M*^2_smg"); aM=v("M2_running_smg"); aB=v("braiding_smg")
D=v("kin (D)"); cs2=v("c_s^2"); rho=v("(.)rho_tot"); p=v("(.)p_tot")
wt=p/rho
q=-1-CubicSpline(N,np.log(H))(N,1)
gH2=-aM*F
slip=aB+2*aM

sp={k:CubicSpline(N,x) for k,x in {
    "F":F,"alphaM":aM,"alphaB":aB,"D":D,"cs2":cs2,"wtot":wt,"q":q,"gH2":gH2,"noslip":slip
}.items()}

def roots_of(key,target):
    y=sp[key](N)-target
    roots=[]
    for i in range(len(N)-1):
        if y[i]==0: roots.append(N[i])
        elif y[i]*y[i+1]<0:
            roots.append(brentq(lambda x: float(sp[key](x)-target),N[i],N[i+1]))
    return roots

qroots=roots_of("q",0.0)
wroots=roots_of("wtot",-1/3)
# focus on latest crossing before today
Nq=max([x for x in qroots if x<=0],default=None)
Nw=max([x for x in wroots if x<=0],default=None)

targets={"F":1.02388542769,"alphaM":0.0,"alphaB":0.0,"D":2.0,"cs2":1.0,
         "q":0.0,"wtot":-1/3,"gH2":0.0,"noslip":0.0}

def state(Nx):
    if Nx is None: return None
    out={"N":float(Nx),"a":float(np.exp(Nx)),"z":float(np.exp(-Nx)-1)}
    for k,s in sp.items():
        out[k]=float(s(Nx))
        out[k+"_abs_distance_to_mature"]=float(abs(s(Nx)-targets[k]))
    return out

sq=state(Nq); sw=state(Nw)
# FLRW consistency between the two crossing locations
deltaN=None if Nq is None or Nw is None else abs(Nq-Nw)

# Measure whether action itself is canonical at the kinematic coasting crossing.
canonical_tolerances={
  "F_frac":0.01,
  "abs_alphaM":1e-3,
  "abs_alphaB":2e-3,
  "D_frac":0.1,
  "cs2_abs":0.1,
  "abs_gH2":1e-3,
}
if sq:
    tests={
      "F":abs(sq["F"]-targets["F"])/targets["F"]<canonical_tolerances["F_frac"],
      "alphaM":abs(sq["alphaM"])<canonical_tolerances["abs_alphaM"],
      "alphaB":abs(sq["alphaB"])<canonical_tolerances["abs_alphaB"],
      "D":abs(sq["D"]-2)/2<canonical_tolerances["D_frac"],
      "cs2":abs(sq["cs2"]-1)<canonical_tolerances["cs2_abs"],
      "gH2":abs(sq["gH2"])<canonical_tolerances["abs_gH2"],
    }
else: tests={}

out={
 "status":"accepted-action kinematic coasting versus action canonicalization audit",
 "scope":"unchanged freely evolved accepted action, supported N<=0 domain only",
 "q_zero_crossings_N":[float(x) for x in qroots],
 "wtot_minus_third_crossings_N":[float(x) for x in wroots],
 "latest_q0_crossing":sq,
 "latest_w_minus_third_crossing":sw,
 "crossing_delta_N":None if deltaN is None else float(deltaN),
 "canonical_tolerances_used_only_as_diagnostic":canonical_tolerances,
 "canonical_tests_at_q0":tests,
 "all_action_sector_tests_pass_at_q0":bool(tests and all(tests.values())),
 "interpretation":{
   "kinematic_statement":"q=0 and w_tot=-1/3 mark total-background coasting.",
   "action_statement":"The mature manuscript branch additionally requires F->F_inf, vanishing running/braiding/G3, D->2 and c_s^2->1.",
   "separation_test":"If coasting occurs while one or more action-sector tests fail strongly, kinematic coasting is not the same event as action canonicalization.",
   "future_limit":"No future extrapolation is used."
 }
}
Path("output").mkdir(exist_ok=True)
Path("output/accepted_action_coasting_crossing.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("ACCEPTED_ACTION_COASTING_CROSSING",json.dumps(out,sort_keys=True))
