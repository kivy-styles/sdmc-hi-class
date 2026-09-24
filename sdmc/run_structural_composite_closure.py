#!/usr/bin/env python3
"""
SDMC structural composite closure audit.

Input: a hi_class background table from the accepted covariant solution.

The covariant action determines the active structural density
  rho_X = 3 Mpl^2 H0^2 [ E^2(F+F') - Omega_m a^-3 - Omega_r a^-4 ].
The original SDMC inverse-square law applies to a bare density
  rho_v = rho_P / S^2,
with rho_X = chi rho_v.

Hence
  p_X = -1/2 d ln rho_X / d ln a = p - A_chi/2.

This script evaluates the minimal-activation closure
  p = max(1,p_X),
  A_chi = 2 max(0,1-p_X),
which is the unique pointwise choice satisfying:
  (i) mature structural floor p>=1,
  (ii) nondecreasing activation A_chi>=0,
  (iii) zero activation whenever the action density can be explained
       by bare structural dilution alone.

This is a model closure hypothesis, not yet a theorem of the Horndeski action.
"""
from pathlib import Path
import argparse, json, math, re
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid

P=argparse.ArgumentParser()
P.add_argument("background")
P.add_argument("--H0",type=float,default=69.45160505326464)
P.add_argument("--S0",type=float,default=2.72e61)
P.add_argument("--rhoP",type=float,default=5.155e96)
P.add_argument("--omega-b",type=float,default=0.02208511574370519)
P.add_argument("--omega-cdm",type=float,default=0.12298509428428875)
P.add_argument("--out",default="structural_composite_closure.json")
a=P.parse_args()

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    x=np.loadtxt(path)
    return {n:x[:,i] for i,n in enumerate(names)}

d=table(a.background)
req=["z","H [1/Mpc]","M*^2_smg","M2_running_smg",
     "(.)rho_g","(.)rho_b","(.)rho_cdm","(.)rho_ur"]
miss=[k for k in req if k not in d]
if miss: raise RuntimeError(f"missing columns {miss}")

z=np.asarray(d["z"])
ln_a=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
F=np.asarray(d["M*^2_smg"])
alphaM=np.asarray(d["M2_running_smg"])
rho_ord=(np.asarray(d["(.)rho_g"])+np.asarray(d["(.)rho_b"]) +
         np.asarray(d["(.)rho_cdm"])+np.asarray(d["(.)rho_ur"]))
rhoX=H*H*F*(1.+alphaM)-rho_ord

o=np.argsort(ln_a)
z,ln_a,H,F,alphaM,rhoX=[x[o] for x in (z,ln_a,H,F,alphaM,rhoX)]
keep=np.r_[True,np.diff(ln_a)>1e-12]
z,ln_a,H,F,alphaM,rhoX=[x[keep] for x in (z,ln_a,H,F,alphaM,rhoX)]
if np.any(rhoX<=0):
    raise RuntimeError("active structural density becomes non-positive")

pX=-0.5*CubicSpline(ln_a,np.log(rhoX))(ln_a,1)
p=np.maximum(1.,pX)
Achi=2.*(p-pX)

# S/S0 from d ln S / d ln a = p, normalized to S(0)=S0.
integ=cumulative_trapezoid(p,ln_a,initial=0.)
i0=int(np.argmin(np.abs(ln_a)))
Srel=np.exp(integ-integ[i0])
ascale=np.exp(ln_a)
Kmrel=Srel/ascale
Kprel=Kmrel

# Lapse N = t_P S H p.
tP=5.391247e-44
MPC_KM=3.085677581491367e19
H0_s=a.H0/MPC_KM
pref=tP*a.S0*H0_s
E=H/H[i0]
Nlapse=pref*Srel*E*p

# Absolute present structural/matter/baryon mapping anchors.
G=6.67430e-11
MPC_M=3.085677581491367e22
H100=100.*1000./MPC_M
rho_c100=3.*H100*H100/(8.*math.pi*G)
rho_m0=(a.omega_b+a.omega_cdm)*rho_c100
rho_b0=a.omega_b*rho_c100
Sm0=(a.rhoP/rho_m0)**(1./3.)
Sb0=(a.rhoP/rho_b0)**(1./3.)
Km0=a.S0/Sm0
Kp0=a.S0/Sb0
fb=a.omega_b/(a.omega_b+a.omega_cdm)

# Present activation normalization from rho_v0=rhoP/S0^2 and
# action-level active fraction rhoX/H0^2.
rho_v0=a.rhoP/a.S0**2
H0SI=a.H0*1000./MPC_M
rho_crit0=3.*H0SI**2/(8.*math.pi*G)
Omega_v_raw0=rho_v0/rho_crit0
Omega_X0=float(rhoX[i0]/H[i0]**2)
chi0=Omega_X0/Omega_v_raw0
chi_rel=(rhoX/rhoX[i0])*Srel*Srel
chi=chi0*chi_rel

# Unique pX=1 crossing separating tracker-like no-activation from
# the mature p=1 activation branch, if one exists.
cross=[]
for i in range(len(ln_a)-1):
    y1,y2=pX[i]-1.,pX[i+1]-1.
    if y1*y2<0:
        x=ln_a[i]+(ln_a[i+1]-ln_a[i])*(-y1)/(y2-y1)
        cross.append(float(math.exp(-x)-1.))

def sample(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return dict(z=float(z[j]),pX=float(pX[j]),p=float(p[j]),
                Achi=float(Achi[j]),S_over_S0=float(Srel[j]),
                N_lapse=float(Nlapse[j]),Km_over_Km0=float(Kmrel[j]),
                Kp_over_Kp0=float(Kprel[j]),chi=float(chi[j]))

out=dict(
  closure="p=max(1,pX); Achi=2 max(0,1-pX)",
  status="minimal activation closure hypothesis; action fixes pX only",
  active_density_positive=True,
  z_pX_equals_1=cross,
  present=dict(
    pX=float(pX[i0]),p=float(p[i0]),Achi=float(Achi[i0]),
    N_lapse=float(Nlapse[i0]),chi0=float(chi0),
    Omega_X0=Omega_X0,Omega_v_raw0=float(Omega_v_raw0),
    Km0=float(Km0),Kp0=float(Kp0),fb=float(fb),Kp_over_Km=float(Kp0/Km0)
  ),
  samples=[sample(q) for q in [1e7,1e6,1e5,1e4,3400,1000,300,100,30,10,6.8,5,3,2,1,.7,.5,.3,.1,0]]
)
Path(a.out).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print(json.dumps(out,indent=2,sort_keys=True))
