#!/usr/bin/env python3
"""
SDMC structural composite closure audit.

Input: a hi_class background table from the accepted covariant solution.

The covariant action determines the active structural density
  rho_X = 3 Mpl^2 H0^2 [ E^2(F+F') - Omega_m a^-3 - Omega_r a^-4 ].
The original SDMC inverse-square law applies to a bare density
  rho_v = rho_P / S^2,
with rho_X = chi rho_v.

Hence
  p_X = -1/2 d ln rho_X / d ln a = p - A_chi/2.

This script evaluates the minimal-activation closure
  p = max(1,p_X),
  A_chi = 2 max(0,1-p_X),
which is the unique pointwise choice satisfying:
  (i) mature structural floor p>=1,
  (ii) nondecreasing activation A_chi>=0,
  (iii) zero activation whenever the action density can be explained
       by bare structural dilution alone.

This is a model closure hypothesis, not yet a theorem of the Horndeski action.
"""
from pathlib import Path
import argparse, json, math, re
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid

P=argparse.ArgumentParser()
P.add_argument("background")
P.add_argument("--H0",type=float,default=69.45160505326464)
P.add_argument("--S0",type=float,default=2.72e61)
P.add_argument("--rhoP",type=float,default=5.155e96)
P.add_argument("--omega-b",type=float,default=0.02208511574370519)
P.add_argument("--omega-cdm",type=float,default=0.12298509428428875)
P.add_argument("--out",default="structural_composite_closure.json")
a=P.parse_args()

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    x=np.loadtxt(path)
    return {n:x[:,i] for i,n in enumerate(names)}

d=table(a.background)
req=["z","proper time [Gyr]","H [1/Mpc]","M*^2_smg","M2_running_smg",
     "(.)rho_g","(.)rho_b","(.)rho_cdm","(.)rho_ur"]
miss=[k for k in req if k not in d]
if miss: raise RuntimeError(f"missing columns {miss}")

z=np.asarray(d["z"])
ln_a=-np.log1p(z)
proper_gyr=np.asarray(d["proper time [Gyr]"])
H=np.asarray(d["H [1/Mpc]"])
F=np.asarray(d["M*^2_smg"])
alphaM=np.asarray(d["M2_running_smg"])
rho_ord=(np.asarray(d["(.)rho_g"])+np.asarray(d["(.)rho_b"]) +
         np.asarray(d["(.)rho_cdm"])+np.asarray(d["(.)rho_ur"]))
rhoX=H*H*F*(1.+alphaM)-rho_ord

o=np.argsort(ln_a)
z,ln_a,proper_gyr,H,F,alphaM,rhoX=[x[o] for x in (z,ln_a,proper_gyr,H,F,alphaM,rhoX)]
keep=np.r_[True,np.diff(ln_a)>1e-12]
z,ln_a,proper_gyr,H,F,alphaM,rhoX=[x[keep] for x in (z,ln_a,proper_gyr,H,F,alphaM,rhoX)]
if np.any(rhoX<=0):
    raise RuntimeError("active structural density becomes non-positive")

pX=-0.5*CubicSpline(ln_a,np.log(rhoX))(ln_a,1)
p=np.maximum(1.,pX)
Achi=2.*(p-pX)

# S/S0 from d ln S / d ln a = p, normalized to S(0)=S0.
integ=cumulative_trapezoid(p,ln_a,initial=0.)
i0=int(np.argmin(np.abs(ln_a)))
Srel=np.exp(integ-integ[i0])
ascale=np.exp(ln_a)
Kmrel=Srel/ascale
Kprel=Kmrel

# Lapse N = t_P S H p.
tP=5.391247e-44
MPC_KM=3.085677581491367e19
H0_s=a.H0/MPC_KM
pref=tP*a.S0*H0_s
E=H/H[i0]
Nlapse=pref*Srel*E*p

# Absolute present structural/matter/baryon mapping anchors.
G=6.67430e-11
MPC_M=3.085677581491367e22
H100=100.*1000./MPC_M
rho_c100=3.*H100*H100/(8.*math.pi*G)
rho_m0=(a.omega_b+a.omega_cdm)*rho_c100
rho_b0=a.omega_b*rho_c100
Sm0=(a.rhoP/rho_m0)**(1./3.)
Sb0=(a.rhoP/rho_b0)**(1./3.)
Km0=a.S0/Sm0
Kp0=a.S0/Sb0
fb=a.omega_b/(a.omega_b+a.omega_cdm)

# Present activation normalization from rho_v0=rhoP/S0^2 and
# action-level active fraction rhoX/H0^2.
rho_v0=a.rhoP/a.S0**2
H0SI=a.H0*1000./MPC_M
rho_crit0=3.*H0SI**2/(8.*math.pi*G)
Omega_v_raw0=rho_v0/rho_crit0
Omega_X0=float(rhoX[i0]/H[i0]**2)
chi0=Omega_X0/Omega_v_raw0
chi_rel=(rhoX/rhoX[i0])*Srel*Srel
chi=chi0*chi_rel

# Chronological closure inherited from Part III: S proportional to cosmic
# time after the nonclassical Planck-to-radiation transition.  Then
# p=d ln S/d ln a=1/(H t_c), S/S0=t_c/t0, and the lapse is constant.
SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
t_mpc=proper_gyr*SEC_PER_GYR*C_MS/MPC_M
p_chrono=1./(H*t_mpc)
Srel_chrono=t_mpc/t_mpc[i0]
Kmrel_chrono=Srel_chrono/ascale
N_chrono=pref*Srel_chrono*E*p_chrono
Achi_chrono=2.*(p_chrono-pX)
chi_rel_chrono=(rhoX/rhoX[i0])*Srel_chrono*Srel_chrono
chi_chrono=chi0*chi_rel_chrono

# Unique pX=1 crossing separating tracker-like no-activation from
# the mature p=1 activation branch, if one exists.
cross=[]
for i in range(len(ln_a)-1):
    y1,y2=pX[i]-1.,pX[i+1]-1.
    if y1*y2<0:
        x=ln_a[i]+(ln_a[i+1]-ln_a[i])*(-y1)/(y2-y1)
        cross.append(float(math.exp(-x)-1.))

def sample(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return dict(
      z=float(z[j]),pX=float(pX[j]),
      minimal=dict(p=float(p[j]),Achi=float(Achi[j]),S_over_S0=float(Srel[j]),
                   N_lapse=float(Nlapse[j]),Km_over_Km0=float(Kmrel[j]),
                   Kp_over_Kp0=float(Kprel[j]),chi=float(chi[j])),
      chronological=dict(p=float(p_chrono[j]),Achi=float(Achi_chrono[j]),
                         S_over_S0=float(Srel_chrono[j]),N_lapse=float(N_chrono[j]),
                         Km_over_Km0=float(Kmrel_chrono[j]),
                         Kp_over_Kp0=float(Kmrel_chrono[j]),chi=float(chi_chrono[j]))
    )

# Radiation-matter analytic bridge from Part III, using the actual equality epoch.
Om=np.asarray(d["Omega_m(z)"])[o][keep] if "Omega_m(z)" in d else None
Or=np.asarray(d["Omega_r(z)"])[o][keep] if "Omega_r(z)" in d else None
rm_diag=None
if Om is not None and Or is not None:
    q=Om-Or
    z_eq=None
    for ii in range(len(ln_a)-1):
        if q[ii]*q[ii+1] <= 0.:
            xe=ln_a[ii]+(ln_a[ii+1]-ln_a[ii])*(-q[ii])/(q[ii+1]-q[ii])
            z_eq=math.exp(-xe)-1.
            break
    if z_eq is not None:
        aeq=1./(1.+z_eq)
        xrm=np.exp(ln_a)/aeq
        urm=np.sqrt(1.+xrm)
        p_rm=3.*(urm+1.)**2/(2.*urm*(urm+2.))
        mm=(z>=100.)&(z<=1.e6)
        rm_diag=dict(z_eq=float(z_eq),
                     max_abs_delta_p=float(np.max(np.abs(p_chrono[mm]-p_rm[mm]))),
                     rms_delta_p=float(np.sqrt(np.mean((p_chrono[mm]-p_rm[mm])**2))),
                     median_delta_p=float(np.median(p_chrono[mm]-p_rm[mm])))

# If the legacy present benchmark N0=3.37 is retained exactly, solve for
# the S0 implied by the actual covariant cosmic age.
N0_benchmark=3.37
t0_s=proper_gyr[i0]*SEC_PER_GYR
S0_for_N337=N0_benchmark*t0_s/tP
S0_ratio_N337=S0_for_N337/a.S0
Km0_N337=Km0*S0_ratio_N337
Kp0_N337=Kp0*S0_ratio_N337
R0_N337_m=1.616255e-35*S0_for_N337

# Current source-derived hierarchy, retained only as a domain-projection
# diagnostic, not literal rest-mass evolution.
lambda0=Km0**(-1./6.)/N_chrono[i0]
source_factor0=lambda0**3

out=dict(
  status="action fixes pX=p-Achi/2; two explicit closures are audited",
  active_density_positive=True,
  z_pX_equals_1=cross,
  anchors=dict(chi0=float(chi0),Omega_X0=Omega_X0,
               Omega_v_raw0=float(Omega_v_raw0),Km0=float(Km0),Kp0=float(Kp0),
               fb=float(fb),Kp_over_Km=float(Kp0/Km0),
               lambda0=float(lambda0),source_factor0=float(source_factor0)),
  radiation_matter_bridge=rm_diag,
  N337_calibration=dict(S0=float(S0_for_N337),S0_ratio_to_2p72e61=float(S0_ratio_N337),
                        R0_m=float(R0_N337_m),Km0=float(Km0_N337),Kp0=float(Kp0_N337)),
  minimal_activation=dict(
    definition="p=max(1,pX); Achi=2 max(0,1-pX)",
    status="pointwise minimal nondecreasing-activation closure with p>=1 floor",
    present=dict(pX=float(pX[i0]),p=float(p[i0]),Achi=float(Achi[i0]),
                 N_lapse=float(Nlapse[i0]),chi=float(chi[i0]))
  ),
  chronological=dict(
    definition="S/S0=t_c/t0; p=1/(H t_c); N=tP*S0/t0",
    status="Part-III chronological-scaling closure; excludes the unresolved Planck-to-classical transition",
    age0_Gyr=float(proper_gyr[i0]),
    present=dict(pX=float(pX[i0]),p=float(p_chrono[i0]),Achi=float(Achi_chrono[i0]),
                 N_lapse=float(N_chrono[i0]),chi=float(chi_chrono[i0])),
    N_lapse_min=float(np.min(N_chrono)),N_lapse_max=float(np.max(N_chrono)),
    Achi_min=float(np.min(Achi_chrono)),
    Achi_max=float(np.max(Achi_chrono)),
    z_Achi_min=float(z[int(np.argmin(Achi_chrono))]),
    z_Achi_max=float(z[int(np.argmax(Achi_chrono))]),
    chi_min=float(np.min(chi_chrono)),
    chi_max=float(np.max(chi_chrono)),
    z_chi_min=float(z[int(np.argmin(chi_chrono))]),
    z_chi_max=float(z[int(np.argmax(chi_chrono))]),
    chi_bounded_0_1=bool(np.all((chi_chrono>=0.)&(chi_chrono<=1.)))
  ),
  samples=[sample(q) for q in [1e7,1e6,1e5,1e4,3400,1090,1000,300,100,30,20,10,6.8,5,3,2,1,.7,.5,.3,.1,0]]
)
# Comparison with the historical Part-II numerical hierarchy.
legacy=dict(N0=3.37,Km0=2.19e20,Kp0=1.18e20,lambda0=1.209e-4,source_factor0=1.766e-12)
out["legacy_comparison_percent"]=dict(
 N0=100.*(N_chrono[i0]/legacy["N0"]-1.),
 Km0=100.*(Km0/legacy["Km0"]-1.),
 Kp0=100.*(Kp0/legacy["Kp0"]-1.),
 lambda0=100.*(lambda0/legacy["lambda0"]-1.),
 source_factor0=100.*(source_factor0/legacy["source_factor0"]-1.),
 lambda0_if_N337=100.*((Km0_N337**(-1./6.)/N0_benchmark)/legacy["lambda0"]-1.),
 source_factor0_if_N337=100.*(((Km0_N337**(-1./6.)/N0_benchmark)**3)/legacy["source_factor0"]-1.)
)

Path(a.out).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print(json.dumps(out,indent=2,sort_keys=True))
