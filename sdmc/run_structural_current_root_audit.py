#!/usr/bin/env python3
"""
Audit whether the accepted structural-field reconstruction is already close to
a shift-current driven nonzero kinetic root.

For the structural field sigma=S/S0 and v=dot(sigma), define

    J = k1 v + 4 B v^3 + 3 H g v^2 - 3 H F_,sigma,
    B = k2/4 - g_,sigma/6,

with the exact homogeneous equation

    dot(J) + 3 H J = ell_,sigma.

A genuinely shift-current selected mature attractor would require the explicit
field source ell_,sigma to become negligible, so that expansion damps J toward
zero.  This script measures:

  * the source-to-Hubble-current ratio |ell_,sigma|/(3 H |J|),
  * the positive algebraic J=0 root relative to the chronological
    v_star=1/t0,
  * the formal H->0, F_,sigma->0 root
        v_asym^2=-k1/(4B),
  * the epoch at which k1 changes sign.

The J=0 roots are diagnostics only when ell_,sigma is non-negligible.
"""
from pathlib import Path
import json, math, re
import numpy as np
from scipy.interpolate import CubicSpline

BG=Path("output/linear_cov_free_00_background.dat")
OUT=Path("output/structural_current_root_audit.json")

AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
TP=5.391247e-44
S0=2.72e61

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

t=tg*SEC_PER_GYR*C_MS/MPC_M
i0=int(np.argmin(np.abs(N)))
t0=float(t[i0])
sigma=t/t0
A=H*t
Apsi=A+h*A*A

# phi -> psi=ln(t/t0)
gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4+2.*g*A**2*Apsi

# psi -> sigma=t/t0
gs=gt/sigma**3
k1s=k1t/sigma**2
k2s=(k2t-2.*gt)/sigma**4

spg=CubicSpline(sigma,gs)
spk1=CubicSpline(sigma,k1s)
spk2=CubicSpline(sigma,k2s)
spV=CubicSpline(sigma,V)
spF=CubicSpline(sigma,F)

g0=spg(sigma); g1s=spg(sigma,1)
q1=spk1(sigma); q1s=spk1(sigma,1)
q2=spk2(sigma)
V1=spV(sigma,1)
F1s=spF(sigma,1); F2s=spF(sigma,2)

B=.25*q2-(1./6.)*g1s
B1=CubicSpline(sigma,B)(sigma,1)

vstar=1./t0
J=q1*vstar+4.*B*vstar**3+3.*H*g0*vstar**2-3.*H*F1s

ell_sigma=(
  .5*q1s*vstar**2
  +B1*vstar**4
  +H*g1s*vstar**3
  -V1
  -3.*H*F2s*vstar
  -3.*H*H*F1s
)

Jdot=H*CubicSpline(N,J)(N,1)
eom=Jdot+3.*H*J-ell_sigma
eombudget=np.abs(Jdot)+3.*H*np.abs(J)+np.abs(ell_sigma)

Jterms=np.vstack([
  q1*vstar,
  4.*B*vstar**3,
  3.*H*g0*vstar**2,
  -3.*H*F1s
])
Jbudget=np.sum(np.abs(Jterms),axis=0)
J_fraction=np.abs(J)/np.maximum(Jbudget,1e-300)
source_to_damping=np.abs(ell_sigma)/np.maximum(3.*H*np.abs(J),1e-300)

# Formal current roots.
root_pos=np.full(len(z),np.nan)
root_asym=np.full(len(z),np.nan)
for i in range(len(z)):
    roots=np.roots([4.*B[i],3.*H[i]*g0[i],q1[i],-3.*H[i]*F1s[i]])
    pos=[float(r.real) for r in roots if abs(r.imag)<1e-11 and r.real>0.]
    if pos:
        root_pos[i]=min(pos,key=lambda x:abs(x-vstar))
    if q1[i]<0. and B[i]>0.:
        root_asym[i]=math.sqrt(-q1[i]/(4.*B[i]))

# k1=0 crossing.
ordz=np.argsort(z)
zz=z[ordz]; kk=q1[ordz]
k1cross=[]
for i in range(len(zz)-1):
    if kk[i]==0.:
        k1cross.append(float(zz[i]))
    elif kk[i]*kk[i+1]<0.:
        k1cross.append(float(zz[i]+(zz[i+1]-zz[i])*(-kk[i])/(kk[i+1]-kk[i])))

conv=C_MS/MPC_M
Nstar=TP*S0*vstar*conv

def sample(zt):
    j=int(np.argmin(np.abs(z-zt)))
    rp=root_pos[j]
    ra=root_asym[j]
    return dict(
      z=float(z[j]),
      sigma=float(sigma[j]),
      k1_sigma=float(q1[j]),
      B_sigma=float(B[j]),
      J_fraction_of_term_budget=float(J_fraction[j]),
      explicit_source_to_3HJ=float(source_to_damping[j]),
      eom_relative_residual=float(abs(eom[j])/max(eombudget[j],1e-300)),
      Jzero_positive_root_over_vstar=(float(rp/vstar) if np.isfinite(rp) else None),
      Hasym_root_over_vstar=(float(ra/vstar) if np.isfinite(ra) else None),
      N_if_Jzero_root=(float(Nstar*rp/vstar) if np.isfinite(rp) else None),
      N_if_Hasym_root=(float(Nstar*ra/vstar) if np.isfinite(ra) else None)
    )

targets=[100,57.46,30,18,17.4,16,10,5,3,1,.5,.1,0]
rows=[sample(x) for x in targets]
out={
  "status":"current-root diagnostics; J=0 is not an attractor condition until the explicit field source becomes negligible",
  "chronological":{"t0_Mpc":t0,"vstar_per_Mpc":vstar,"Nstar":Nstar},
  "k1_zero_crossings":k1cross,
  "present":sample(0.),
  "rows":rows
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("STRUCTURAL_CURRENT_ROOT_AUDIT")
print("K1_ZERO",json.dumps(k1cross))
print("PRESENT",json.dumps(out["present"],sort_keys=True))
for row in rows:
    print("CURRENT_ROOT_SAMPLE",json.dumps(row,sort_keys=True))
