#!/usr/bin/env python3
"""
Construct C2 future coefficient-tail designs that join the accepted structural
action at sigma=1 and approach the mature inverse-square normal form.

This is a design audit, not yet a native hi_class evolution.

For Q in {k1,k2,V}, define x=ln(sigma) and Qbar=e^(2x) Q.  The tail

  Qbar = Qinf + exp(-mu x) [A0 + A1 x + A2 x^2/2]

matches Qbar, dQbar/dx and d2Qbar/dx2 at x=0 exactly when

  A0 = Qbar0-Qinf
  A1 = Qbar0' + mu A0
  A2 = Qbar0'' + 2 mu Qbar0' + mu^2 A0.

F is stitched directly to constant Finf with the same C2 construction.

The target kinetic branch is joined with
  Z = Zinf + (Z0-Zinf) exp(-muZ x)(1+muZ x+(muZ x)^2/2),
which preserves Z, Z' and Z'' at the present constant-Z endpoint.

No-Slip is kept exact on the target by defining
  g = -F_,sigma/(2Z).
"""
from pathlib import Path
import json,re,math
import numpy as np
from scipy.interpolate import CubicSpline

BG=Path("output/linear_cov_free_00_background.dat")
OUT=Path("output/future_action_stitch_design.json")

AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366
CS2_REF=0.138335

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
XI=math.sqrt(8.*math.pi/3.)
N0=3.4135084646534866
F_INF=math.exp(AF)
R=(1./CS2_REF-1.)/2.

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
z=np.asarray(d["z"]); N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d["(.)rho_smg"]); pre=np.asarray(d["(.)p_smg"])
tg=np.asarray(d["proper time [Gyr]"])
o=np.argsort(N)
N,H,rho,pre,tg,z=[x[o] for x in (N,H,rho,pre,tg,z)]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre,tg,z=[x[keep] for x in (N,H,rho,pre,tg,z)]

lnH=CubicSpline(N,np.log(H)); h=lnH(N,1); dotH=H*H*h
X=.5*H*H
Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH); tt=np.tanh(xx); uu=1.-tt*tt
ST=.5*(1.+tt); S1=uu/(4.*WIDTH); S2=-tt*uu/(4.*WIDTH*WIDTH)
F=np.exp(AF*ST); F1=F*(AF*S1); F2=F*((AF*S1)**2+AF*S2)
g=-F1/(H*H)
g1=CubicSpline(N,g)(N,1)
Dtarget=DFLOOR+D0*ST**POWER
alphaM=F1/F; alphaB=-2.*alphaM
alphaK=Dtarget-1.5*alphaB*alphaB

RR=3.*rho+2.*g1*X*X+6.*F1*H*H+3.*(F-1.)*H*H
PP=3.*pre+2.*g1*X*X-2.*X*F2 \
   -3.*(F-1.)*H*H-2.*(F-1.)*dotH-2.*F1*(H*H+dotH)
Cbg=(RR+PP)/(2.*X)
k2=(F*alphaK-Cbg+4.*g1*X-6.*g*H*H)/(4.*X)
k1=Cbg-2.*k2*X
V=.5*(RR-PP)-k2*X*X

t=tg*SEC_PER_GYR*C_MS/MPC_M
i0=int(np.argmin(np.abs(N)))
t0=float(t[i0])
sigma=t/t0
x=np.log(sigma)
A=H*t
Apsi=A+h*A*A

# phi -> psi -> sigma
gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4+2.*g*A**2*Apsi
k1s=k1t/sigma**2
k2s=(k2t-2.*gt)/sigma**4
Vs=V

Z0=1./(2.*t0*t0)

def endpoint_triplet(y):
    sp=CubicSpline(x,y)
    return [float(sp(0.0,nu)) for nu in [0,1,2]]

bars={
  "k1":endpoint_triplet(sigma*sigma*k1s),
  "k2":endpoint_triplet(sigma*sigma*k2s),
  "V":endpoint_triplet(sigma*sigma*Vs),
  "F":endpoint_triplet(F)
}

def tail_coeff(q0,q1,q2,qinf,mu):
    A0=q0-qinf
    A1=q1+mu*A0
    A2=q2+2.*mu*q1+mu*mu*A0
    return [A0,A1,A2]

def eval_tail(xv,qinf,abc,mu,rescaled=True):
    A0,A1,A2=abc
    qbar=qinf+math.exp(-mu*xv)*(A0+A1*xv+.5*A2*xv*xv)
    return qbar*math.exp(-2.*xv) if rescaled else qbar

def dFdx(xv,Finf,abc,mu):
    A0,A1,A2=abc
    p=A0+A1*xv+.5*A2*xv*xv
    pp=A1+A2*xv
    return math.exp(-mu*xv)*(pp-mu*p)

scenarios=[]
for name,chiinf in [("bounded_chi",1.0),("matched_chi_F",F_INF)]:
    Ninf=XI*math.sqrt(chiinf/F_INF)
    zratio=(Ninf/N0)**2
    Zinf=Z0*zratio
    kap1=2.*F_INF*(1.-R)
    kap2=R*F_INF/Zinf
    U0=F_INF*Zinf*(4.-R)
    # Use the activation-flow gamma as one concrete transition rate.
    m0=0.06774426095348086
    C=(2.-m0)/m0
    gamma=math.log1p(C)/math.log(chiinf/0.9410895692368776)
    mu=gamma
    muZ=gamma
    asym={"k1":kap1,"k2":kap2,"V":U0,"F":F_INF}
    abc={q:tail_coeff(*bars[q],asym[q],mu) for q in ["k1","k2","V","F"]}
    samples=[]
    for xv in [0.,.01,.02,.03,.05,.1,.2,.5,1.]:
        Ez=math.exp(-muZ*xv)*(1.+muZ*xv+.5*(muZ*xv)**2)
        Zt=Zinf+(Z0-Zinf)*Ez
        sig=math.exp(xv)
        Ft=eval_tail(xv,F_INF,abc["F"],mu,rescaled=False)
        Fx=dFdx(xv,F_INF,abc["F"],mu)
        Fs=Fx/sig
        gtarg=-Fs/(2.*Zt)
        q1=eval_tail(xv,kap1,abc["k1"],mu,True)
        q2=eval_tail(xv,kap2,abc["k2"],mu,True)
        vv=eval_tail(xv,U0,abc["V"],mu,True)
        csproxy=(q1+2.*q2*Zt)/(q1+6.*q2*Zt)
        samples.append({
          "x":xv,"sigma":sig,"Z_over_Z0":Zt/Z0,"N_over_N0":math.sqrt(Zt/Z0),
          "F":Ft,"g_target":gtarg,"k1":q1,"k2":q2,"V":vv,
          "kessence_cs2_proxy":csproxy
        })
    scenarios.append({
      "name":name,"chi_inf":chiinf,"F_inf":F_INF,
      "N_inf":Ninf,"N_inf_over_N0":Ninf/N0,
      "Z_inf_over_Z0":zratio,"r":R,"gamma":gamma,
      "asymptotic_coefficients":{"kappa1":kap1,"kappa2":kap2,"U0":U0},
      "endpoint_rescaled_value_derivatives":bars,
      "tail_A_coefficients":abc,
      "samples":samples
    })

out={
  "status":"C2 coefficient-tail design only; full Horndeski D and cs2 must be checked by native evolution",
  "present":{"t0_Mpc":t0,"Z0":Z0,"N0":N0,"F0":float(F[i0])},
  "scenarios":scenarios
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_ACTION_STITCH_DESIGN")
for s in scenarios:
    print("STITCH",json.dumps({
      "name":s["name"],"N_inf":s["N_inf"],"Z_inf_over_Z0":s["Z_inf_over_Z0"],
      "gamma":s["gamma"],"r":s["r"],
      "cs2_proxy_x1":s["samples"][-1]["kessence_cs2_proxy"]
    },sort_keys=True))
