#!/usr/bin/env python3
from pathlib import Path
import numpy as np,re,json
from scipy.interpolate import CubicSpline

TARGET=Path("output/linear_cov_target_00_background.dat")
FREE=Path("output/linear_cov_free_00_background.dat")
INI=Path("output/linear_cov_target.ini")

def ini_value(path,key):
    for raw in Path(path).read_text().splitlines():
        line=raw.split("#",1)[0].strip()
        if not line or "=" not in line: continue
        k,v=line.split("=",1)
        if k.strip()==key: return v.strip()
    raise RuntimeError(key)
kin=[float(x.strip()) for x in ini_value(INI,"parameters_smg").split(",")]
AF,ZC,W,D0,POWER,DFLOOR=kin

def read(path):
    ls=Path(path).read_text().splitlines()
    hdr=[l for l in ls if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); ns=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        ns.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(ns)}

t=read(TARGET)
N=-np.log1p(t["z"]); o=np.argsort(N)
N=N[o]; H=np.asarray(t["H [1/Mpc]"])[o]
rho=np.asarray(t["(.)rho_smg"])[o]; pre=np.asarray(t["(.)p_smg"])[o]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre=N[keep],H[keep],rho[keep],pre[keep]
z=np.exp(-N)-1.

lnH=CubicSpline(N,np.log(H)); h=lnH(N,1); Hpa=H*H*h; X=.5*H*H
Nc=-np.log1p(ZC); xx=(N-Nc)/(2.*W); tt=np.tanh(xx); uu=1.-tt*tt
S=.5*(1.+tt); S1=uu/(4.*W); S2=-tt*uu/(4.*W*W)
F=np.exp(AF*S); F1=F*(AF*S1); F2=F*((AF*S1)**2+AF*S2)
g=-F1/(H*H); g1=CubicSpline(N,g)(N,1)
Dtar=DFLOOR+D0*S**POWER
aM=F1/F; aB=-2*aM; aK=Dtar-1.5*aB*aB
R=3*rho+2*g1*X*X+6*F1*H*H+3*(F-1.)*H*H
P=3*pre+2*g1*X*X-2*X*F2-3*(F-1.)*H*H-2*(F-1.)*Hpa-2*F1*(H*H+Hpa)
C=(R+P)/(2*X)
k2=(F*aK-C+4*g1*X-6*g*H*H)/(4*X)
k1=C-2*k2*X
V=.5*(R-P)-k2*X*X
if not np.all(V>0): raise RuntimeError(f"V non-positive min={V.min()}")
pV=-.5*CubicSpline(N,np.log(V))(N,1)

f=read(FREE)
Nf=-np.log1p(f["z"]); of=np.argsort(Nf)
Nf=Nf[of]; Hf=np.asarray(f["H [1/Mpc]"])[of]
Ff=np.asarray(f["M*^2_smg"])[of]; aMf=np.asarray(f["M2_running_smg"])[of]
rn=np.asarray(f["rho_tot_wo_smg_dbg"])[of]
keepf=np.r_[True,np.diff(Nf)>1e-12]
Nf,Hf,Ff,aMf,rn=[x[keepf] for x in (Nf,Hf,Ff,aMf,rn)]
zf=np.exp(-Nf)-1
rhoX=Hf*Hf*Ff*(1+aMf)-rn
pX=-.5*CubicSpline(Nf,np.log(rhoX))(Nf,1)
pXi=np.interp(N,Nf,pX)
rhoXi=np.interp(N,Nf,rhoX)

ratio=rhoXi/V
dP=pV-pXi
Achi_X=2*(1-pXi)
Achi_V=2*(1-pV)

def win(zlo,zhi):
    m=(z>=zlo)&(z<=zhi)
    return {
      "z_window":[zlo,zhi],"n":int(m.sum()),
      "pX_median":float(np.median(pXi[m])),
      "pV_median":float(np.median(pV[m])),
      "pV_minus_pX_median":float(np.median(dP[m])),
      "pV_minus_pX_rms":float(np.sqrt(np.mean(dP[m]**2))),
      "rhoX_over_V_median":float(np.median(ratio[m])),
      "rhoX_over_V_min":float(np.min(ratio[m])),
      "rhoX_over_V_max":float(np.max(ratio[m])),
      "Achi_rawp1_from_X_median":float(np.median(Achi_X[m])),
      "Achi_rawp1_from_V_median":float(np.median(Achi_V[m])),
    }

def at(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return {
      "z":float(z[j]),"pX":float(pXi[j]),"pV":float(pV[j]),
      "rhoX_over_V":float(ratio[j]),
      "Achi_X_if_rawp1":float(Achi_X[j]),
      "Achi_V_if_rawp1":float(Achi_V[j]),
      "F":float(F[j]),"V_over_H2":float(V[j]/H[j]**2),
      "k1X_over_V":float(k1[j]*X[j]/V[j]),
      "k2X2_over_V":float(k2[j]*X[j]*X[j]/V[j]),
    }

out={
 "status":"potential-factorization diagnostic; V is not assumed identical to rhoX",
 "radiation":win(1e5,1e7),
 "matter":win(30,300),
 "transition":win(1,5),
 "late":win(0,0.2),
 "today":at(0),
 "samples":[at(x) for x in [300,100,30,10,5,3,2,1,.5,.1,0]]
}
Path("output/late266_activation_factorization_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("ACTIVATION_FACTORIZATION_AUDIT")
for k in ["radiation","matter","transition","late"]:
    print(k,json.dumps(out[k],sort_keys=True))
print("TODAY",json.dumps(out["today"],sort_keys=True))
for row in out["samples"]: print("ACTIVATION_SAMPLE",json.dumps(row,sort_keys=True))
