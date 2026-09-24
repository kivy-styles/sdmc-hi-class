#!/usr/bin/env python3
from pathlib import Path
import numpy as np, re, json
from scipy.interpolate import CubicSpline

SRC=Path("output/linear_cov_free_00_background.dat")
if not SRC.exists():
    raise SystemExit("missing free covariant background")

lines=SRC.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
for i,m in enumerate(ms):
    e=ms[i+1].start() if i+1<len(ms) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(SRC)
d={n:a[:,i] for i,n in enumerate(names)}

required=["z","H [1/Mpc]","M*^2_smg","M2_running_smg","rho_tot_wo_smg_dbg","(.)rho_smg"]
missing=[x for x in required if x not in d]
if missing:
    raise RuntimeError(f"missing columns {missing}; available={names}")

z=np.asarray(d["z"])
N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
F=np.asarray(d["M*^2_smg"])
alphaM=np.asarray(d["M2_running_smg"])
rho_non=np.asarray(d["rho_tot_wo_smg_dbg"])
rho_smg_book=np.asarray(d["(.)rho_smg"])

# IMPORTANT:
# hi_class rho_smg is an effective Jordan-frame bookkeeping density and is
# not the raw SDMC inverse-square reservoir.  The positive structural-source
# combination used in the manuscript's scalar-tensor Friedmann reconstruction is
#
#   rho_X = H^2 (F + F') - rho_nonSMG
#         = H^2 F (1+alpha_M) - rho_nonSMG .
#
# This is the ACTIVE source rho_X = chi rho_v.  It must not be used directly as
# rho_v when reconstructing the raw SDMC stretch.
rhoX=H*H*F*(1.+alphaM)-rho_non

order=np.argsort(N)
N,H,F,alphaM,rho_non,rho_smg_book,rhoX,z=[
    np.asarray(x)[order] for x in (N,H,F,alphaM,rho_non,rho_smg_book,rhoX,z)
]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,F,alphaM,rho_non,rho_smg_book,rhoX,z=[
    x[keep] for x in (N,H,F,alphaM,rho_non,rho_smg_book,rhoX,z)
]

if not np.all(rhoX>0):
    raise RuntimeError(f"active structural source not positive: min={rhoX.min()}")

lnrho=CubicSpline(N,np.log(rhoX))
pX=-0.5*lnrho(N,1)  # half-slope of ACTIVE density, not raw p

i0=int(np.argmin(np.abs(N)))
H0=H[i0]
rhoX0=rhoX[i0]
OmegaX0=rhoX0/(H0*H0)

# Raw present SDMC normalization from the manuscript's S0.
S0_work=2.72e61
rhoP_si=5.15500e96
G_si=6.67430e-11
MPC_M=3.085677581491367e22
H0_kms=69.45160505326464
H0_si=H0_kms*1000./MPC_M
rho_crit_si=3.*H0_si*H0_si/(8.*np.pi*G_si)
rho_raw0_si=rhoP_si/(S0_work*S0_work)
OmegaRaw0=rho_raw0_si/rho_crit_si
chi0=OmegaX0/OmegaRaw0

# Exact identity from rho_X = chi rho_v and rho_v proportional S^-2:
#
#   p_X = p_raw - A_chi/2
#   p_raw = p_X + A_chi/2.
#
# The action determines p_X.  A separate activation law/dynamical identity is
# still needed to determine p_raw, and hence N, Km and Kp, globally.
#
# As a DIAGNOSTIC only, test the manuscript's mature raw branch p_raw=1.
praw_late=1.0
Achi_if_praw1=2.*(praw_late-pX)

# Under p_raw=1, S/S0=a and chi can be reconstructed directly from rho_X.
a_scale=np.exp(N)
Xnorm=rhoX/(H0*H0)
chi_if_praw1=Xnorm*a_scale*a_scale/OmegaRaw0

# For exact radiation and matter raw trackers, compare the action-derived active
# half-slope with p_raw=2 and 3/2.  The inferred A_chi should be near zero if the
# structural sector is merely tracking rather than activating.
Achi_if_rad=2.*(2.0-pX)
Achi_if_mat=2.*(1.5-pX)

def stats(zlo,zhi,rawp=None):
    m=(z>=zlo)&(z<=zhi)&np.isfinite(pX)
    if not np.any(m): return None
    out={
      "z_window":[zlo,zhi],"n":int(np.sum(m)),
      "pX_median":float(np.median(pX[m])),
      "pX_mean":float(np.mean(pX[m])),
      "pX_min":float(np.min(pX[m])),
      "pX_max":float(np.max(pX[m])),
    }
    if rawp is not None:
      achi=2.*(rawp-pX[m])
      out["raw_p_hypothesis"]=float(rawp)
      out["Achi_inferred_median"]=float(np.median(achi))
      out["Achi_inferred_mean"]=float(np.mean(achi))
    return out

def sample(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return {
      "z":float(z[j]),"Nlog":float(N[j]),"rhoX_active":float(rhoX[j]),
      "rho_smg_bookkeeping":float(rho_smg_book[j]),
      "pX_active_half_slope":float(pX[j]),
      "Achi_if_raw_p1":float(Achi_if_praw1[j]),
      "chi_if_raw_p1":float(chi_if_praw1[j]),
      "F":float(F[j]),"alphaM":float(alphaM[j])
    }

out={
 "status":"rhoX is active scalar-tensor source, not raw rho_v",
 "rho_smg_bookkeeping_positive_fraction":float(np.mean(rho_smg_book>0)),
 "rhoX_positive_fraction":float(np.mean(rhoX>0)),
 "rhoX_min":float(np.min(rhoX)),
 "OmegaX0_from_action":float(OmegaX0),
 "OmegaRaw0_from_S0":float(OmegaRaw0),
 "chi0_from_action_over_raw":float(chi0),
 "identity":"p_raw = pX + Achi/2",
 "today":sample(0.0),
 "radiation_tracker_test":stats(1e5,1e7,2.0),
 "pre_equality_test":stats(5e3,5e4,None),
 "matter_tracker_test":stats(30,300,1.5),
 "late_raw_p1_test":stats(0,0.2,1.0),
 "transition_raw_p1_diagnostic":stats(1,5,1.0),
 "samples":[sample(x) for x in [1e7,1e6,1e5,1e4,3400,1000,300,100,30,10,3,1,.5,.1,0]]
}
Path("output/late266_structural_composite_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("STRUCTURAL_COMPOSITE_AUDIT_V2")
print("rho_smg_bookkeeping_positive_fraction",out["rho_smg_bookkeeping_positive_fraction"])
print("rhoX_positive_fraction",out["rhoX_positive_fraction"])
print("OmegaX0_from_action",OmegaX0)
print("OmegaRaw0_from_S0",OmegaRaw0)
print("chi0_from_action_over_raw",chi0)
for key in ["radiation_tracker_test","pre_equality_test","matter_tracker_test","late_raw_p1_test","transition_raw_p1_diagnostic"]:
    print(key,json.dumps(out[key],sort_keys=True))
print("TODAY",json.dumps(out["today"],sort_keys=True))
for row in out["samples"]:
    print("STRUCTURAL_SAMPLE",json.dumps(row,sort_keys=True))
