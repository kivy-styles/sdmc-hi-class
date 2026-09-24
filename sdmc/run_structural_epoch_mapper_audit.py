#!/usr/bin/env python3
"""
Epoch-dependent SDMC matter/baryon mapper audit for the accepted free
covariant background.

Chronological closure:
    S/S0 = t/t0

Conserved matter:
    Km/Km0 = (S/S0)/a
    Kp/Kp0 = (S/S0)/a
    Kp/Km   = f_b^(1/3)

The script records exact interpolated values at equality, last scattering,
drag, and representative late epochs.  It is intentionally separate from any
acoustic-ruler projection: these are structural/domain mappers.
"""
from pathlib import Path
import json,re
import numpy as np

BG=Path("output/linear_cov_free_00_background.dat")
OUT=Path("output/structural_epoch_mappers.json")
S0=2.72e61
Km0=2.1992619194015187e20
Kp0=1.1743136282444456e20
fb=0.15223742861779455

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
z=np.asarray(d["z"])
t=np.asarray(d["proper time [Gyr]"])
o=np.argsort(z)
z,t=z[o],t[o]
t0=float(np.interp(0.,z,t))

def epoch(zt,label):
    tt=float(np.interp(zt,z,t))
    aa=1./(1.+zt)
    srel=tt/t0
    krel=srel/aa
    km=Km0*krel
    kp=Kp0*krel
    return dict(label=label,z=float(zt),a=float(aa),t_Gyr=tt,
                S_over_S0=float(srel),S=float(S0*srel),
                Km_over_Km0=float(krel),Kp_over_Kp0=float(krel),
                Km=float(km),Kp=float(kp),
                Kp_over_Km=float(kp/km))

targets=[
    (3466.50,"accepted equality"),
    (1090.0,"last scattering"),
    (1059.0,"representative drag"),
    (1000.0,"post recombination reference"),
    (100.0,"matter tracker"),
    (10.0,"tracker handoff vicinity"),
    (3.4328776988,"Planck-mass transition center"),
    (0.0,"present"),
]
rows=[epoch(zv,l) for zv,l in targets]
out={
  "definition":{
    "Km":"S*(rho_m/rho_mP)^(1/3)",
    "common_primordial_reference":"rho_mP=rho_P",
    "density_only":"Km=(rho_P*rho_m^2/rho_v^3)^(1/6)",
    "chronological":"Km/Km0=(S/S0)/a",
    "baryon":"Kp=Km*f_b^(1/3)"
  },
  "anchors":{"S0":S0,"Km0":Km0,"Kp0":Kp0,"fb":fb,
             "fb_one_third":fb**(1./3.)},
  "epochs":rows,
  "note":"Structural/domain mappers only; do not automatically reapply the absolute Kp to a sound horizon already computed in the local photon-baryon thermal domain."
}
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("STRUCTURAL_EPOCH_MAPPERS")
for r in rows:
    print(json.dumps(r,sort_keys=True))
