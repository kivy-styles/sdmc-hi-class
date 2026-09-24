#!/usr/bin/env python3
"""
Mature-coasting normal-form audit for the SDMC structural scalar.

Candidate late structural action (after F -> constant and No-Slip braiding
decays) is

    G2 = sigma^-2 [ kappa1 Z + kappa2 Z^2 - U0 ],
    G3 -> 0,
    G4 -> F_inf/2,

with sigma=S/S0 and Z=dot(sigma)^2/2.

For a constant-velocity synchronized orbit sigma proportional to a proportional
to t, the scalar equation implies

    U0 = 2 kappa1 Z_* + 3 kappa2 Z_*^2.

The stress tensor then obeys exactly

    rho = 3 Z_* (kappa1 + 2 kappa2 Z_*) sigma^-2,
    p   =  -Z_* (kappa1 + 2 kappa2 Z_*) sigma^-2,
    w   = -1/3.

If this sector dominates a flat constant-F future,

    kappa1 + 2 kappa2 Z_* = 2 F_inf,

and

    c_s^2 = F_inf / (F_inf + 2 kappa2 Z_*).

This script checks how close the present reconstructed action is to the
sigma^-2 normal form and records the future lapse/activation consistency
condition

    N_inf = Xi_P sqrt(chi_inf/F_inf)

for p=1 and scalar domination.
"""
from pathlib import Path
import json, math, re
import numpy as np
from scipy.interpolate import CubicSpline

BG=Path("output/linear_cov_free_00_background.dat")
OUT=Path("output/mature_coasting_normal_form_audit.json")

AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
G=6.67430e-11
TP=5.391247e-44
RHOP=5.155e96
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

# Accepted phi -> psi -> sigma field redefinition.
gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4+2.*g*A**2*Apsi
gs=gt/sigma**3
k1s=k1t/sigma**2
k2s=(k2t-2.*gt)/sigma**4

# Proximity to sigma^-2 scale-covariant normal form.
bar_k1=sigma**2*k1s
bar_k2=sigma**2*k2s
bar_V=sigma**2*V
lnsig=np.log(sigma)

def slope(y):
    return CubicSpline(lnsig,np.log(np.abs(y)))(lnsig,1)

s_k1=slope(bar_k1)
s_k2=slope(bar_k2)
s_V=slope(bar_V)
s_F=CubicSpline(lnsig,np.log(F))(lnsig,1)

# Present structural closure anchors from the accepted action.
p0=1.034215658562453
chi0=0.9410895692368776
OmegaX0=0.7237401509603739
N0=3.4135084646534866
XiP=math.sqrt((8.*math.pi*G/3.)*RHOP*TP*TP)
N_density_now=p0*XiP*math.sqrt(chi0/OmegaX0)

# Future scalar-dominated constant-F relation.
F0=float(F[i0])
required_chi_F1=(N0/XiP)**2
required_chi_F0=F0*(N0/XiP)**2
Ninf_if_chi1_F1=XiP
Ninf_if_chi1_F0=XiP/math.sqrt(F0)

def row(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return dict(
      z=float(z[j]),sigma=float(sigma[j]),
      sigma2_k1=float(bar_k1[j]),dln_sigma2k1_dlnsigma=float(s_k1[j]),
      sigma2_k2=float(bar_k2[j]),dln_sigma2k2_dlnsigma=float(s_k2[j]),
      sigma2_V=float(bar_V[j]),dln_sigma2V_dlnsigma=float(s_V[j]),
      F=float(F[j]),dlnF_dlnsigma=float(s_F[j]),
      alphaM=float(alphaM[j]),alphaB=float(alphaB[j])
    )

out={
  "candidate_normal_form":{
    "G2":"sigma^-2*(kappa1*Z+kappa2*Z^2-U0)",
    "G3":"0 asymptotically",
    "G4":"F_inf/2",
    "constant_velocity_condition":"U0=2*kappa1*Zstar+3*kappa2*Zstar^2",
    "equation_of_state":"w=-1/3 exactly on synchronized constant-Z orbit",
    "flat_scalar_dominated_closure":"kappa1+2*kappa2*Zstar=2*F_inf",
    "sound_speed":"cs2=F_inf/(F_inf+2*kappa2*Zstar)"
  },
  "present_identity":{
    "XiP":XiP,"N0":N0,"p0":p0,"chi0":chi0,"OmegaX0":OmegaX0,
    "N_from_density_identity":N_density_now
  },
  "future_lapse_consistency":{
    "if_Finf_1_and_chiinf_1_Ninf":Ninf_if_chi1_F1,
    "if_Finf_equals_present_F_and_chiinf_1_Ninf":Ninf_if_chi1_F0,
    "chiinf_required_to_keep_N0_if_Finf_1":required_chi_F1,
    "chiinf_required_to_keep_N0_if_Finf_present_F":required_chi_F0,
    "fractional_drop_N_if_chiinf1_Finf1":1.-Ninf_if_chi1_F1/N0
  },
  "present_proximity":row(0.),
  "rows":[row(x) for x in [10,5,3,2,1,.5,.3,.1,0]]
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("MATURE_COASTING_NORMAL_FORM_AUDIT")
print("PRESENT_IDENTITY",json.dumps(out["present_identity"],sort_keys=True))
print("FUTURE_LAPSE",json.dumps(out["future_lapse_consistency"],sort_keys=True))
print("PRESENT_PROXIMITY",json.dumps(out["present_proximity"],sort_keys=True))
