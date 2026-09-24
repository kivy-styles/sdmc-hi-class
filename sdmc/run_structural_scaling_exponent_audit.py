#!/usr/bin/env python3
"""
Structural scaling-exponent audit.

For the accepted action-level active density rho_X and chronological structural
field sigma=S/S0, define

    p     = d ln S / d ln a = 1/(H t),
    p_X   = -1/2 d ln rho_X / d ln a,
    A_chi = 2(p-p_X),

and the structural activation beta function

    beta_chi = d ln chi / d ln sigma = A_chi/p.

Then the effective exponent of the active density with respect to sigma is

    m_X = - d ln rho_X / d ln sigma = 2 p_X/p
        = 2 - beta_chi.

The mature inverse-square fixed point is m_X -> 2, beta_chi -> 0.
A cosmological-constant-like active source has m_X -> 0, beta_chi -> 2.
"""
from pathlib import Path
import json,re
import numpy as np
from scipy.interpolate import CubicSpline

BG=Path("output/linear_cov_free_00_background.dat")
OUT=Path("output/structural_scaling_exponent_audit.json")

AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22

def read(path):
    lines=path.read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"(?:^|\s)1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    x=np.loadtxt(path)
    return {n:x[:,i] for i,n in enumerate(names)}

d=read(BG)
for k in ["z","proper time [Gyr]","H [1/Mpc]","(.)rho_smg","(.)p_smg"]:
    if k not in d: raise RuntimeError(f"missing {k}")

z=np.asarray(d["z"]); N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d["(.)rho_smg"])
pre=np.asarray(d["(.)p_smg"])
tg=np.asarray(d["proper time [Gyr]"])
o=np.argsort(N)
N,H,rho,pre,tg,z=[x[o] for x in (N,H,rho,pre,tg,z)]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre,tg,z=[x[keep] for x in (N,H,rho,pre,tg,z)]

lnH=CubicSpline(N,np.log(H)); h=lnH(N,1); dotH=H*H*h
X=.5*H*H
Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH)
tt=np.tanh(xx); uu=1.-tt*tt
ST=.5*(1.+tt); S1=uu/(4.*WIDTH); S2=-tt*uu/(4.*WIDTH*WIDTH)
F=np.exp(AF*ST); F1=F*(AF*S1); F2=F*((AF*S1)**2+AF*S2)

g=-F1/(H*H)
g1=CubicSpline(N,g)(N,1)
Dtarget=DFLOOR+D0*ST**POWER
alphaM=F1/F; alphaB=-2.*alphaM
alphaK=Dtarget-1.5*alphaB*alphaB

R=3.*rho+2.*g1*X*X+6.*F1*H*H+3.*(F-1.)*H*H
P=3.*pre+2.*g1*X*X-2.*X*F2 \
  -3.*(F-1.)*H*H-2.*(F-1.)*dotH \
  -2.*F1*(H*H+dotH)
Cbg=(R+P)/(2.*X)
k2=(F*alphaK-Cbg+4.*g1*X-6.*g*H*H)/(4.*X)
k1=Cbg-2.*k2*X
V=.5*(R-P)-k2*X*X

# Positive action-level source, Eq. (455), with phi=N so dot(phi)=H.
rhoX=k1*X+3.*k2*X*X+V+6.*H*H*X*g-2.*X*X*g1-3.*H*H*F1
if np.any(rhoX<=0):
    raise RuntimeError("rho_X must remain positive for logarithmic slope audit")

pX=-.5*CubicSpline(N,np.log(rhoX))(N,1)

t=tg*SEC_PER_GYR*C_MS/MPC_M
p=1./(H*t)
Achi=2.*(p-pX)
beta=Achi/p
mX=2.*pX/p
wX=-1.+2.*pX/3.

identity=np.max(np.abs(mX-(2.-beta)))

def sample(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return dict(
      z=float(z[j]),
      p=float(p[j]),pX=float(pX[j]),Achi=float(Achi[j]),
      beta_chi=float(beta[j]),mX_structural=float(mX[j]),
      wX=float(wX[j]),rhoX=float(rhoX[j])
    )

targets=[1e6,1e5,1e4,3466.5,1090,100,30,18.3,16.17,12.7,10,5,3.43,3.01,2,1,.5,.3,.1,0]
out={
  "definitions":{
    "beta_chi":"Achi/p = d ln chi / d ln sigma",
    "mX":"2*pX/p = -d ln rhoX/d ln sigma = 2-beta_chi",
    "wX":"-1 + 2*pX/3"
  },
  "identity_max_abs_error":float(identity),
  "present":sample(0.),
  "rows":[sample(x) for x in targets]
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("STRUCTURAL_SCALING_EXPONENT_AUDIT")
print("IDENTITY_ERR",identity)
print("PRESENT",json.dumps(out["present"],sort_keys=True))
for r in out["rows"]:
    print("SCALING_SAMPLE",json.dumps(r,sort_keys=True))
