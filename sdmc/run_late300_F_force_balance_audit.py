#!/usr/bin/env python3
from pathlib import Path
import re,json,math
import numpy as np
from scipy.interpolate import CubicSpline

ROOT=Path("recovered"); OUT=Path("output"); OUT.mkdir(exist_ok=True)
ini=next(ROOT.rglob("linear_cov_target.ini"))
bgp=next(ROOT.rglob("linear_cov_target_00_background.dat"))

def get_list(key):
    for raw in ini.read_text().splitlines():
        q=raw.split("#",1)[0].strip()
        if q.startswith(key) and "=" in q:
            return [float(x.strip()) for x in q.split("=",1)[1].split(",")]
    raise RuntimeError(key)

p=get_list("parameters_smg")
AF,ZC,W,D0,POW,DF=p[:6]

ls=bgp.read_text().splitlines()
hdr=[x for x in ls if x.startswith("#") and re.search(r"(?:^|\s)1\s*:",x)][-1].lstrip("#").strip()
ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
for i,m in enumerate(ms):
    e=ms[i+1].start() if i+1<len(ms) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(bgp)
if a.ndim==1:a=a[None,:]
d={n:a[:,i] for i,n in enumerate(names)}
N=-np.log1p(np.asarray(d["z"],float)); order=np.argsort(N)
z=np.asarray(d["z"],float)[order]; N=N[order]
H=np.asarray(d["H [1/Mpc]"],float)[order]
rho=np.asarray(d["(.)rho_smg"],float)[order]
pre=np.asarray(d["(.)p_smg"],float)[order]
keep=np.r_[True,np.diff(N)>1e-12]
z,N,H,rho,pre=z[keep],N[keep],H[keep],rho[keep],pre[keep]

lnH=CubicSpline(N,np.log(H)); h=lnH(N,1)
Hpa=H*H*h
X=.5*H*H
Nc=-math.log1p(ZC)
xx=(N-Nc)/(2*W); tt=np.tanh(xx); uu=1-tt*tt
S=.5*(1+tt); S1=uu/(4*W); S2=-tt*uu/(4*W*W)
F=np.exp(AF*S); F1=F*(AF*S1); F2=F*((AF*S1)**2+AF*S2)
g=-F1/(H*H)
gsp=CubicSpline(N,g); gp=gsp(N,1)

D=DF+D0*S**POW
alphaM=F1/F; alphaB=-2*alphaM
alphaK=D-1.5*alphaB*alphaB

# Reconstruct covariant functions exactly as accepted audit.
Raux=3*rho + 2*gp*X*X + 6*F1*H*H + 3*(F-1)*H*H
Paux=3*pre + 2*gp*X*X - 2*X*F2 - 3*(F-1)*H*H - 2*(F-1)*Hpa - 2*F1*(H*H+Hpa)
C=(Raux+Paux)/(2*X)
k2=(F*alphaK-C+4*gp*X-6*g*H*H)/(4*X)
k1=C-2*k2*X
V=.5*(Raux-Paux)-k2*X*X

k1p=CubicSpline(N,k1)(N,1)
k2p=CubicSpline(N,k2)(N,1)
Vp=CubicSpline(N,V)(N,1)

G2phi=k1p*X+k2p*X*X-Vp
G3phi=gp*X
G4phi=.5*F1
phi2_over_a2=H*H*(1+h)

# Exact surviving pieces of hi_class Shift for this restricted Horndeski model.
curv=12*G4phi*H*H + 6*G4phi*Hpa
intrinsic=G2phi + 2*G3phi*H*H + G3phi*phi2_over_a2
shift=curv+intrinsic
CR=np.abs(curv)/np.maximum(np.abs(intrinsic),1e-300)

def roots_level(y, level, mask):
    yy=np.log(np.maximum(y,1e-300))-math.log(level)
    out=[]
    ids=np.where(mask)[0]
    for i in ids[:-1]:
        if not mask[i+1]: continue
        y0,y1=yy[i],yy[i+1]
        if not(np.isfinite(y0) and np.isfinite(y1)): continue
        if y0==0: NN=N[i]
        elif y0*y1<0:
            f=-y0/(y1-y0); NN=N[i]+f*(N[i+1]-N[i])
        else: continue
        out.append({"ln_a":float(NN),"z":float(math.exp(-NN)-1)})
    return out

m=(z>=0.05)&(z<=30)&np.isfinite(CR)&(CR>0)
levels={str(v):roots_level(CR,v,m) for v in [0.25,0.5,1.0,2.0,4.0]}
jzc=int(np.argmin(np.abs(z-ZC)))
jpk=int(np.argmax(np.where(m,alphaM,-np.inf)))
# local extrema / slope of CR
logcr=np.log(np.maximum(CR,1e-300))
dcr=np.gradient(logcr,N,edge_order=2)
jsteep=np.where(m)[0][np.argmax(np.abs(dcr[m]))]

# Evaluate action force pieces at accepted center.
out={
 "status":"exact hi_class Shift decomposition of late300 curvature loading",
 "formula":{
   "curvature":"12 G4_phi H^2 + 6 G4_phi Hprime/a = G4_phi R4",
   "intrinsic":"G2_phi + 2 G3_phi H phi_prime/a + G3_phi phi_prime_prime/a^2",
   "CR":"abs(curvature)/abs(intrinsic)"
 },
 "accepted":{"A_F":AF,"z_c":ZC,"width":W},
 "at_accepted_zc":{
   "z":float(z[jzc]),"CR":float(CR[jzc]),"curvature":float(curv[jzc]),
   "intrinsic":float(intrinsic[jzc]),"shift":float(shift[jzc]),
   "alphaM":float(alphaM[jzc]),"F":float(F[jzc])
 },
 "alphaM_peak":{"z":float(z[jpk]),"ln_a":float(N[jpk]),"alphaM":float(alphaM[jpk]),"CR":float(CR[jpk])},
 "CR_level_crossings":levels,
 "steepest_logCR":{"z":float(z[jsteep]),"ln_a":float(N[jsteep]),"dlnCR_dlnA":float(dcr[jsteep]),"CR":float(CR[jsteep])}
}
# choose the CR=1 crossing closest to accepted center only as a diagnostic, not derivation
r1=levels["1.0"]
if r1:
    best=min(r1,key=lambda r:abs(r["z"]-ZC))
    out["CR1_closest_to_accepted_center"]=best
    out["CR1_fractional_z_difference"]=(best["z"]-ZC)/ZC
    Nc_force=float(best["ln_a"])
    slope_force=float(np.interp(Nc_force,N,dcr))
    out["CR1_local_slope_dlnCR_dlnA"]=slope_force
    out["CR1_logistic_width_inverse_slope"]=1.0/abs(slope_force) if slope_force!=0 else None
    out["accepted_width_fractional_difference_from_force_width"]=(1.0/abs(slope_force)-W)/W if slope_force!=0 else None
Path("output/F_force_balance_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("F_FORCE_BALANCE_AUDIT",json.dumps(out,sort_keys=True))
