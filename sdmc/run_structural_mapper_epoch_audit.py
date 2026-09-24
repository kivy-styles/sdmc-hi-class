#!/usr/bin/env python3
"""
Epoch-dependent SDMC matter/baryon mapper audit for the chronological closure.

Definitions inherited from the manuscript:
    S/S0 = t/t0
    Sm = (rhoP/rho_m)^(1/3),   Km = S/Sm
    Sb = (rhoP/rho_b)^(1/3),   Kp = S/Sb

For conserved nonrelativistic matter and baryons,
    rho_m = rho_m0 a^-3, rho_b = rho_b0 a^-3,
so
    Sm/Sm0 = Sb/Sb0 = a
and therefore
    Km/Km0 = Kp/Kp0 = (S/S0)/a.

This script evaluates the exact accepted-action background at selected epochs,
including last scattering and baryon drag, and checks the density and
S-based forms against each other.
"""
from pathlib import Path
import argparse, json, math, re
import numpy as np

P=argparse.ArgumentParser()
P.add_argument("background")
P.add_argument("--H0",type=float,default=69.45160505326464)
P.add_argument("--S0",type=float,default=2.72e61)
P.add_argument("--rhoP",type=float,default=5.155e96)
P.add_argument("--omega-b",type=float,default=0.02208511574370519)
P.add_argument("--omega-cdm",type=float,default=0.12298509428428875)
P.add_argument("--zstar",type=float,default=1090.0)
P.add_argument("--zdrag",type=float,default=1059.0)
P.add_argument("--zeq",type=float,default=3466.5021)
P.add_argument("--out",default="output/structural_mapper_epoch_audit.json")
a=P.parse_args()

def read(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    x=np.loadtxt(path)
    return {n:x[:,i] for i,n in enumerate(names)}

d=read(a.background)
for k in ["z","proper time [Gyr]","H [1/Mpc]"]:
    if k not in d: raise RuntimeError(f"missing required column {k}")

z=np.asarray(d["z"])
t=np.asarray(d["proper time [Gyr]"])
H=np.asarray(d["H [1/Mpc]"])
o=np.argsort(z)
z,t,H=z[o],t[o],H[o]

# Unique ascending z grid.
keep=np.r_[True,np.diff(z)>1e-12]
z,t,H=z[keep],t[keep],H[keep]

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
G=6.67430e-11
tP=5.391247e-44
H100=100.*1000./MPC_M
rho_c100=3.*H100*H100/(8.*math.pi*G)
rho_m0=(a.omega_b+a.omega_cdm)*rho_c100
rho_b0=a.omega_b*rho_c100
Sm0=(a.rhoP/rho_m0)**(1./3.)
Sb0=(a.rhoP/rho_b0)**(1./3.)
Km0=a.S0/Sm0
Kp0=a.S0/Sb0
fb=a.omega_b/(a.omega_b+a.omega_cdm)

# Present accepted-action age from the table.
t0=float(np.interp(0.0,z,t))
Nchrono=tP*a.S0/(t0*SEC_PER_GYR)

def sample(zt):
    tt=float(np.interp(zt,z,t))
    HH=float(np.interp(zt,z,H))
    aa=1./(1.+zt)
    Srel=tt/t0
    S=a.S0*Srel
    fac=Srel/aa
    rhom=rho_m0/aa**3
    rhob=rho_b0/aa**3

    Km_anchor=Km0*fac
    Kp_anchor=Kp0*fac
    Km_direct=S*(rhom/a.rhoP)**(1./3.)
    Kp_direct=S*(rhob/a.rhoP)**(1./3.)

    t_mpc=tt*SEC_PER_GYR*C_MS/MPC_M
    p=1./(HH*t_mpc)

    km16=Km_anchor**(-1./6.)
    kp16=Kp_anchor**(-1./6.)
    lam=km16/Nchrono
    return dict(
        z=float(zt), a=float(aa), age_Gyr=tt, S_over_S0=float(Srel),
        p=float(p), mapper_relative=float(fac),
        Km=float(Km_anchor), Kp=float(Kp_anchor),
        Km_over_Km0=float(Km_anchor/Km0),
        Kp_over_Kp0=float(Kp_anchor/Kp0),
        Kp_over_Km=float(Kp_anchor/Km_anchor),
        fb_one_third=float(fb**(1./3.)),
        Km_minus_one_sixth=float(km16),
        Kp_minus_one_sixth=float(kp16),
        lambda_source=float(lam),
        source_factor_lambda3=float(lam**3),
        Km_direct_from_S_rhom=float(Km_direct),
        Kp_direct_from_S_rhob=float(Kp_direct),
        Km_identity_relerr=float(Km_direct/Km_anchor-1.),
        Kp_identity_relerr=float(Kp_direct/Kp_anchor-1.)
    )

targets=[0.0,0.3,1.0,3.0,10.0,100.0,a.zdrag,a.zstar,a.zeq,1e4,1e5,1e6]
rows=[sample(x) for x in targets]

out=dict(
    definition=dict(
        chronological="S/S0=t/t0",
        Km="S*(rho_m/rhoP)^(1/3)",
        Kp="S*(rho_b/rhoP)^(1/3)",
        conserved_form="Km/Km0 = Kp/Kp0 = (S/S0)/a",
        evolution="d ln Km/d ln a = d ln Kp/d ln a = p-1",
        direct_ratio="Kp/Km = f_b^(1/3)"
    ),
    anchors=dict(
        age0_Gyr=t0,N_lapse=Nchrono,S0=a.S0,
        Km0=Km0,Kp0=Kp0,fb=fb,Kp_over_Km=Kp0/Km0
    ),
    selected=dict(
        last_scattering=sample(a.zstar),
        drag_epoch=sample(a.zdrag),
        equality=sample(a.zeq)
    ),
    rows=rows
)
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
Path(a.out).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("STRUCTURAL_MAPPER_EPOCH_AUDIT")
print("ANCHORS",json.dumps(out["anchors"],sort_keys=True))
for key in ["last_scattering","drag_epoch","equality"]:
    print(key.upper(),json.dumps(out["selected"][key],sort_keys=True))
