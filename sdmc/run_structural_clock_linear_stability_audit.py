#!/usr/bin/env python3
"""
Reduced homogeneous linear-stability audit of the SDMC structural-clock orbit.

Starting from the accepted free linear-G3 background, reconstruct the same
covariant coefficients, field-redefine first to psi=ln(t/t0) and then to
sigma=t/t0=S/S0, and evaluate the reduced minisuperspace equation

    K_hom * dot(v) + C(sigma,v;t) = 0,   v = dot(sigma).

The chronological orbit has v=v_* = 1/t0 and dot(v_*)=0.  Holding the metric
background fixed for a local scalar-only diagnostic gives

    delta_dot_sigma = delta_v
    delta_dot_v = - M_*^2 delta_sigma - Gamma_* delta_v

with
    Gamma_* = C_,v / K_hom
    M_*^2   = C_,sigma / K_hom.

This is NOT a substitute for the full hi_class basin test, because scalar
perturbations also backreact on H, F and the matter fractions.  It is a
Lagrangian-level diagnostic of the constant-kinetic orbit.
"""
from pathlib import Path
import json, re
import numpy as np
from scipy.interpolate import CubicSpline

BG=Path("output/linear_cov_free_00_background.dat")
OUT=Path("output/structural_clock_linear_stability_audit.json")

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
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

d=read(BG)
for k in ["z","proper time [Gyr]","H [1/Mpc]","(.)rho_smg","(.)p_smg"]:
    if k not in d: raise RuntimeError(f"missing {k}")

z=np.asarray(d["z"])
N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d["(.)rho_smg"])
pre=np.asarray(d["(.)p_smg"])
tg=np.asarray(d["proper time [Gyr]"])
o=np.argsort(N)
N,H,rho,pre,tg,z=[x[o] for x in (N,H,rho,pre,tg,z)]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre,tg,z=[x[keep] for x in (N,H,rho,pre,tg,z)]

# Reconstruct accepted linear-G3 action in phi ~= ln a.
lnH=CubicSpline(N,np.log(H))
h=lnH(N,1)
dotH=H*H*h
X=.5*H*H

Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH)
tt=np.tanh(xx); uu=1.-tt*tt
ST=.5*(1.+tt)
S1=uu/(4.*WIDTH)
S2=-tt*uu/(4.*WIDTH*WIDTH)
F=np.exp(AF*ST)
F1=F*(AF*S1)
F2=F*((AF*S1)**2+AF*S2)

g=-F1/(H*H)
g1=CubicSpline(N,g)(N,1)
Dtarget=DFLOOR+D0*ST**POWER
alphaM=F1/F
alphaB=-2.*alphaM
alphaK=Dtarget-1.5*alphaB*alphaB

R=3.*rho+2.*g1*X*X+6.*F1*H*H+3.*(F-1.)*H*H
P=3.*pre+2.*g1*X*X-2.*X*F2 \
  -3.*(F-1.)*H*H-2.*(F-1.)*dotH \
  -2.*F1*(H*H+dotH)
Cbg=(R+P)/(2.*X)
k2=(F*alphaK-Cbg+4.*g1*X-6.*g*H*H)/(4.*X)
k1=Cbg-2.*k2*X
V=.5*(R-P)-k2*X*X

# Cosmic time in Mpc (c=1), then psi=ln(t/t0), sigma=t/t0.
t=tg*SEC_PER_GYR*C_MS/MPC_M
i0=int(np.argmin(np.abs(N)))
t0=float(t[i0])
psi=np.log(t/t0)
sigma=t/t0
A=H*t
Apsi=A+h*A*A

# phi -> psi.
gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4+2.*g*A**2*Apsi

# psi -> sigma, dpsi/dsigma=1/sigma.
gs=gt/sigma**3
k1s=k1t/sigma**2
k2s=(k2t-2.*gt)/sigma**4
Vs=V
Fs=F

# Smooth functions of sigma and derivatives.
spg=CubicSpline(sigma,gs)
spk1=CubicSpline(sigma,k1s)
spk2=CubicSpline(sigma,k2s)
spV=CubicSpline(sigma,Vs)
spF=CubicSpline(sigma,Fs)

g0=spg(sigma); g1s=spg(sigma,1); g2s=spg(sigma,2)
q1=spk1(sigma); q1s=spk1(sigma,1); q1ss=spk1(sigma,2)
q2=spk2(sigma)
V1=spV(sigma,1); V2=spV(sigma,2)
F1s=spF(sigma,1); F2s=spF(sigma,2)

B=.25*q2-(1./6.)*g1s
spB=CubicSpline(sigma,B)
B1=spB(sigma,1); B2=spB(sigma,2)

v=1./t0

Khom=q1+12.*B*v*v+6.*H*g0*v

Cforce=(
  V1
  +.5*q1s*v*v
  +3.*B1*v**4
  +3.*H*q1*v
  +12.*H*B*v**3
  +2.*H*g1s*v**3
  +3.*(dotH+3.*H*H)*g0*v*v
  -3.*(dotH+2.*H*H)*F1s
)

Cv=(
  q1s*v
  +12.*B1*v**3
  +3.*H*q1
  +36.*H*B*v*v
  +6.*H*g1s*v*v
  +6.*(dotH+3.*H*H)*g0*v
)

Cs=(
  V2
  +.5*q1ss*v*v
  +3.*B2*v**4
  +3.*H*q1s*v
  +12.*H*B1*v**3
  +2.*H*g2s*v**3
  +3.*(dotH+3.*H*H)*g1s*v*v
  -3.*(dotH+2.*H*H)*F2s
)

Gamma=Cv/Khom
M2=Cs/Khom
gammaH=Gamma/H
m2H2=M2/(H*H)

disc=Gamma*Gamma-4.*M2
sqrt_disc=np.sqrt(disc.astype(complex))
lam_plus=(-Gamma+sqrt_disc)/2.
lam_minus=(-Gamma-sqrt_disc)/2.
max_re_over_H=np.maximum(np.real(lam_plus),np.real(lam_minus))/H

# Constant-orbit equation residual relative to the absolute term budget.
terms=np.vstack([
  V1,
  .5*q1s*v*v,
  3.*B1*v**4,
  3.*H*q1*v,
  12.*H*B*v**3,
  2.*H*g1s*v**3,
  3.*(dotH+3.*H*H)*g0*v*v,
  -3.*(dotH+2.*H*H)*F1s
])
budget=np.sum(np.abs(terms),axis=0)
closure_rel=np.abs(Cforce)/np.maximum(budget,1e-300)

# Locate Gamma/H = 0 crossing in z <= 100.
m100=(z<=100.)
zs=z[m100]; gg=gammaH[m100]
oo=np.argsort(zs); zs=zs[oo]; gg=gg[oo]
cross=[]
for i in range(len(zs)-1):
    if gg[i]==0.:
        cross.append(float(zs[i]))
    elif gg[i]*gg[i+1] < 0.:
        x=zs[i]+(zs[i+1]-zs[i])*(-gg[i])/(gg[i+1]-gg[i])
        cross.append(float(x))

def sample(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return dict(
      z=float(z[j]),
      sigma=float(sigma[j]),
      K_hom=float(Khom[j]),
      Gamma_over_H=float(gammaH[j]),
      M2_over_H2=float(m2H2[j]),
      max_Re_lambda_over_H=float(max_re_over_H[j]),
      constant_orbit_closure_rel=float(closure_rel[j])
    )

sel=[1e6,1e5,1e4,3466.5,1090,100,57,30,10,5,3,2,1,.5,.3,.1,0]
rows=[sample(x) for x in sel]

late=m100
out={
  "status":"reduced scalar-only frozen-background stability diagnostic; full homogeneous basin test remains decisive",
  "t0_Mpc":t0,
  "v_star_per_Mpc":v,
  "z_Gamma_crossings_z100":cross,
  "z100":{
    "K_hom_min":float(np.min(Khom[late])),
    "K_hom_max":float(np.max(Khom[late])),
    "Gamma_over_H_min":float(np.min(gammaH[late])),
    "Gamma_over_H_max":float(np.max(gammaH[late])),
    "M2_over_H2_min":float(np.min(m2H2[late])),
    "M2_over_H2_max":float(np.max(m2H2[late])),
    "max_constant_orbit_closure_rel":float(np.max(closure_rel[late]))
  },
  "present":sample(0.),
  "samples":rows
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("STRUCTURAL_CLOCK_LINEAR_STABILITY")
print("Z100",json.dumps(out["z100"],sort_keys=True))
print("GAMMA_CROSS",json.dumps(cross))
print("PRESENT",json.dumps(out["present"],sort_keys=True))
for row in rows:
    print("STABILITY_SAMPLE",json.dumps(row,sort_keys=True))
