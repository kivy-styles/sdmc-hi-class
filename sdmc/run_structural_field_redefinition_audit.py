#!/usr/bin/env python3
"""
Audit a monotonic field redefinition from the accepted covariant scalar
phi ~= ln a to the SDMC structural clock field

    sigma = S/S0 = t_c/t0

on the classical cosmological branch.

Original action:
  G4 = F(phi)/2
  G3 = g(phi) X_phi
  G2 = k1(phi) X_phi + k2(phi) X_phi^2 - V(phi)

For phi=f(sigma), A=df/dsigma, X_phi=A^2 X_sigma,

  G4_t = F/2
  G3_t = [g A^3] X_sigma
  G2_t = [k1 A^2] X_sigma
       + [k2 A^4 + 2 g A^2 A_sigma] X_sigma^2 - V.

Thus the same polynomial/linear Horndeski ansatz is closed under this
field redefinition.  On sigma=t/t0, X_sigma=1/(2 t0^2) is constant.

This is an on-branch structural-field reconstruction, not a proof of a
unique Planck-boundary completion.
"""
from pathlib import Path
import re, json
import numpy as np
from scipy.interpolate import CubicSpline

PATH=Path("output/linear_cov_target_00_background.dat")

def read(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}
d=read(PATH)

AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366

z=np.asarray(d["z"])
N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d["(.)rho_smg"])
pre=np.asarray(d["(.)p_smg"])
t_gyr=np.asarray(d["proper time [Gyr]"])
order=np.argsort(N)
N,H,rho,pre,t_gyr,z=[x[order] for x in (N,H,rho,pre,t_gyr,z)]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre,t_gyr,z=[x[keep] for x in (N,H,rho,pre,t_gyr,z)]

lnHsp=CubicSpline(N,np.log(H))
h=lnHsp(N,1)
Hpa=H*H*h
Xphi=.5*H*H

Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH)
tt=np.tanh(xx); uu=1.-tt*tt
ST=.5*(1.+tt)
S1=uu/(4.*WIDTH)
S2=-tt*uu/(4.*WIDTH*WIDTH)
F=np.exp(AF*ST)
F1=F*(AF*S1)
F2=F*((AF*S1)**2+AF*S2)

# Accepted linear-G3 action reconstruction.
g=-F1/(H*H)
g1=CubicSpline(N,g)(N,1)
Dtarget=DFLOOR+D0*ST**POWER
alphaM=F1/F
alphaB=-2.*alphaM
alphaK_target=Dtarget-1.5*alphaB*alphaB
R=3.*rho + 2.*g1*Xphi*Xphi + 6.*F1*H*H + 3.*(F-1.)*H*H
P=3.*pre + 2.*g1*Xphi*Xphi - 2.*Xphi*F2 \
  - 3.*(F-1.)*H*H - 2.*(F-1.)*Hpa - 2.*F1*(H*H+Hpa)
C=(R+P)/(2.*Xphi)
k2=(F*alphaK_target-C+4.*g1*Xphi-6.*g*H*H)/(4.*Xphi)
k1=C-2.*k2*Xphi
V=.5*(R-P)-k2*Xphi*Xphi

# sigma = t/t0, using proper time converted to Mpc (c=1).
SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
t_mpc=t_gyr*SEC_PER_GYR*C_MS/MPC_M
i0=int(np.argmin(np.abs(N)))
t0=t_mpc[i0]
sigma=t_mpc/t0

# A=dphi/dsigma=t0 H because phi=ln a and dphi/dt=H.
A=t0*H
# dA/dsigma=t0^2 dH/dt = t0^2 H^2 h.
Asigma=t0*t0*H*H*h
Xsigma=np.full_like(H,0.5/(t0*t0))

gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4 + 2.*g*A**2*Asigma
Vt=V
Ft=F

# Algebraic checks.
Xphi_from_sigma=A*A*Xsigma
x_map_rel=np.max(np.abs(Xphi_from_sigma/Xphi-1.))
noslip_old=Xphi*g + 0.5*F1
# F_sigma = F_phi * dphi/dsigma = F1*A.
Fsigma=F1*A
noslip_new=Xsigma*gt + 0.5*Fsigma

G3old=g*Xphi
G3new=gt*Xsigma
G2old=k1*Xphi+k2*Xphi*Xphi-V
G2new=k1t*Xsigma+k2t*Xsigma*Xsigma-V
# Under field redefinition:
# G3_t = A G3_old; G2_t = G2_old + 2 A_sigma X_sigma G3_old.
g3_map=np.max(np.abs(G3new-A*G3old)/(1.+np.abs(A*G3old)))
g2_rhs=G2old+2.*Asigma*Xsigma*G3old
g2_map=np.max(np.abs(G2new-g2_rhs)/(1.+np.abs(g2_rhs)))

# Background structural identities.
p=1./(H*t_mpc)
Nlap_const=5.391247e-44*2.72e61/(t_gyr[i0]*SEC_PER_GYR)

# Compare scalar redefinition directly to free/target phi=ln a convention.
# sigma should be monotonic and A positive.
out={
 "status":"on-branch monotonic field redefinition audit",
 "t0_Gyr":float(t_gyr[i0]),
 "sigma_min":float(np.min(sigma)),
 "sigma_max":float(np.max(sigma)),
 "A_min":float(np.min(A)),
 "A_max":float(np.max(A)),
 "Xsigma_t0sq":float(Xsigma[i0]*t0*t0),
 "X_map_max_rel":float(x_map_rel),
 "old_noslip_max_abs":float(np.max(np.abs(noslip_old))),
 "new_noslip_max_abs":float(np.max(np.abs(noslip_new))),
 "G3_transform_max_scaled_abs":float(g3_map),
 "G2_transform_max_scaled_abs":float(g2_map),
 "N_lapse_constant":float(Nlap_const),
 "present":{
   "sigma":float(sigma[i0]),"p":float(p[i0]),"A_dphi_dsigma":float(A[i0]),
   "F":float(Ft[i0]),"gtilde_Xsigma":float(G3new[i0]),
   "k1tilde_Xsigma":float(k1t[i0]*Xsigma[i0]),
   "k2tilde_Xsigma2":float(k2t[i0]*Xsigma[i0]**2),
   "V_over_H2":float(Vt[i0]/H[i0]**2)
 }
}
Path("output/structural_field_redefinition_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("STRUCTURAL_FIELD_REDEFINITION",json.dumps(out,sort_keys=True))
