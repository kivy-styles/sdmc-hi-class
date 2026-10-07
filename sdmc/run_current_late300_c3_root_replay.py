#!/usr/bin/env python3
"""
Replay the manuscript C3 Planck/braiding tail-root equation on the CURRENT
late300 accepted frozen action.

This intentionally reuses the historical structural-field map and the exact
Balanced-Identity kinetic normalization recovered from the historical
successful C3 audit.  It changes only the accepted late-time background/action
jet.

Historical equation:
  mu_F A_F3(mu_F) = 2 Z_BI A_g3(mu_F+1)
with
  A_F3 = F3 + 3 mu F2 + 3 mu^2 F1 + mu^3(F0-Finf)
  A_g3 = g3 + 3(mu+1)g2 + 3(mu+1)^2 g1 + (mu+1)^3 g0.

The structural conversion is the same one used in
run_future_homogeneous_covariant_stitch_audit.py:
  sigma=t/t0, A=H t, g_sigma = g_phi A^3 / sigma^3.
"""
from pathlib import Path
import json, math, re
import numpy as np
from scipy.interpolate import CubicSpline

ROOT=Path("inputs/current_action")
OUT=Path("output/current_c3_root_replay")
OUT.mkdir(parents=True,exist_ok=True)

TARGET=next(ROOT.rglob("linear_cov_target.ini"))
BG=next(ROOT.rglob("linear_cov_free_00_background.dat"))

# Exact value recovered from the historical successful C3 audit.
ZBI=2.063631290689503e-08
XI=math.sqrt(8.0*math.pi/3.0)
OLD_ROOTS=[2.8610674139117047,5.573218997214239]
OLD_FINF=1.02388542769

SEC_PER_GYR=1e9*365.25*86400.0
C_MS=299792458.0
MPC_M=3.085677581491367e22

def ini_value(key):
    for raw in TARGET.read_text().splitlines():
        line=raw.split("#",1)[0].strip()
        if not line or "=" not in line: continue
        k,v=line.split("=",1)
        if k.strip()==key: return v.strip()
    raise RuntimeError(f"missing {key}")

kin=[float(x.strip()) for x in ini_value("parameters_smg").split(",")]
if len(kin)!=6: raise RuntimeError(f"expected six parameters_smg values, got {kin}")
AF,ZC,WIDTH,D0,POWER,DFLOOR=kin

def read_bg(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"(?:^|\s)1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

d=read_bg(BG)
for key in ["z","proper time [Gyr]","H [1/Mpc]"]:
    if key not in d: raise RuntimeError(f"background missing {key}")

N=-np.log1p(np.asarray(d["z"],float))
H=np.asarray(d["H [1/Mpc]"],float)
tg=np.asarray(d["proper time [Gyr]"],float)
order=np.argsort(N)
N,H,tg=N[order],H[order],tg[order]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,tg=N[keep],H[keep],tg[keep]

t=tg*SEC_PER_GYR*C_MS/MPC_M
i0=int(np.argmin(np.abs(N)))
t0=float(t[i0]); H0=float(H[i0]); N_today=float(N[i0])
sigma=t/t0

# Current accepted action's Planck-factor trajectory and exact on-trajectory
# No-Slip linear-G3 coefficient.
lnH=CubicSpline(N,np.log(H))
hN=lnH(N,1)
Nc=-math.log1p(ZC)
xx=(N-Nc)/(2.0*WIDTH)
tt=np.tanh(xx); uu=1.0-tt*tt
S=0.5*(1.0+tt)
S1=uu/(4.0*WIDTH)
F=np.exp(AF*S)
F1=F*(AF*S1)
gphi=-F1/(H*H)

# Historical structural-field map: phi -> psi=ln sigma -> sigma.
A=H*t
gt=gphi*A**3
gs=gt/sigma**3

# Remove any duplicated sigma nodes.
oo=np.argsort(sigma)
sig=sigma[oo]; FF=F[oo]; GG=gs[oo]
kk=np.r_[True,np.diff(sig)>1e-14]
sig,FF,GG=sig[kk],FF[kk],GG[kk]

spF=CubicSpline(sig,FF)
spg=CubicSpline(sig,GG)

def xjet(sp):
    # x=ln sigma; at sigma=1, d/dx = sigma d/dsigma.
    q0=float(sp(1.0)); q1=float(sp(1.0,1)); q2=float(sp(1.0,2)); q3=float(sp(1.0,3))
    return [q0,q1,q1+q2,q1+3.0*q2+q3]

Fjet=xjet(spF)
gjet=xjet(spg)

def root_audit(Finf):
    F0,F1x,F2x,F3x=Fjet
    g0,g1x,g2x,g3x=gjet
    mu=np.poly1d([1.0,0.0]); u=np.poly1d([1.0,1.0])
    AF3=np.poly1d([F0-Finf,3.0*F1x,3.0*F2x,F3x])
    Ag3=g0*u**3+3.0*g1x*u**2+3.0*g2x*u+g3x
    P=mu*AF3-2.0*ZBI*Ag3
    roots=np.roots(P)
    rows=[]
    for z in roots:
        zr=float(np.real(z)); zi=float(np.imag(z))
        row={"real":zr,"imag":zi,"is_real":abs(zi)<1e-9,
             "is_positive_real":abs(zi)<1e-9 and zr>0.0}
        if row["is_positive_real"]:
            af=float(AF3(zr)); ag=float(Ag3(zr))
            row.update(mu_F=zr,mu_g=zr+1.0,A_F3=af,A_g3=ag,
                       root_residual=float(P(zr)),
                       Z_from_tail_ratio=float(zr*af/(2.0*ag)))
        rows.append(row)
    pos=sorted([r for r in rows if r["is_positive_real"]],key=lambda r:r["real"])
    return {
      "F_inf":Finf,
      "polynomial_coefficients_descending":[float(x) for x in P.c],
      "roots":rows,
      "positive_real_roots":[float(r["real"]) for r in pos],
      "positive_root_count":len(pos)
    }

source_Finf=math.exp(AF)
manuscript=root_audit(OLD_FINF)
source_sat=root_audit(source_Finf)

Z0=0.5/t0**2
N0_struct_recalibrated=XI*math.sqrt(Z0/ZBI)
p0=1.0/(H0*t0)

out={
 "status":"current late300 accepted-action C3 Planck/braiding tail-root replay",
 "method":{
   "structural_map":"sigma=t/t0; A=H*t; g_sigma=g_phi*A^3/sigma^3",
   "root_equation":"mu A_F3(mu)=2 Z_BI A_g3(mu+1)",
   "Z_BI":ZBI,
   "historical_positive_roots":OLD_ROOTS,
   "historical_F_inf":OLD_FINF,
   "field_normalization":"identical to historical structural-clock attractor audit"
 },
 "current_action":{
   "H0_km_s_Mpc":float(ini_value("H0")),
   "A_F":AF,"z_c":ZC,"width":WIDTH,"D0":D0,"power":POWER,"D_floor":DFLOOR,
   "source_saturation_F_inf":source_Finf,
   "today_table_N":N_today,
   "t0_Mpc_over_c":t0,
   "H0_1_Mpc":H0,
   "p0_from_chronological_structural_map":p0,
   "Z0":Z0,
   "N0_struct_recalibrated_from_retained_ZBI":N0_struct_recalibrated,
   "F_x_jet":Fjet,
   "g_x_jet":gjet
 },
 "retained_manuscript_mature_endpoint":manuscript,
 "diagnostic_source_saturation_endpoint":source_sat,
 "verdict":{
   "old_roots_can_be_reused_without_recalculation":False,
   "retained_mature_endpoint_has_positive_C3_routes":manuscript["positive_root_count"]>0,
   "route_count_retained_mature_endpoint":manuscript["positive_root_count"],
   "interpretation":"The current accepted late300 present jet has been inserted into the same historical C3 root equation. The retained manuscript mature endpoint is the scientific test; the source-saturation endpoint is shown only as a diagnostic because the source future-domain audit already proved that source parameterization is not the mature completion."
 }
}
(OUT/"current_late300_c3_root_replay.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("CURRENT_LATE300_C3_ROOT_REPLAY",json.dumps(out,sort_keys=True))
