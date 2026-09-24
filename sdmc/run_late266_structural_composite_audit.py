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

rho_name="(.)rho_smg"
if rho_name not in d:
    raise RuntimeError(f"{rho_name} missing; available={names}")

z=np.asarray(d["z"])
N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d[rho_name])
order=np.argsort(N)
N,H,rho,z=N[order],H[order],rho[order],z[order]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,z=N[keep],H[keep],rho[keep],z[keep]

positive=rho>0
frac_positive=float(np.mean(positive))
if not np.all(positive):
    print("STRUCTURAL_RHO_NONPOSITIVE",int(np.sum(~positive)),"of",len(rho),flush=True)
    # Use only positive connected range containing today for logarithmic reconstruction.
    i0=int(np.argmin(np.abs(N)))
    lo=i0
    while lo>0 and rho[lo-1]>0: lo-=1
    N,H,rho,z=N[lo:],H[lo:],rho[lo:],z[lo:]

lnrho=CubicSpline(N,np.log(rho))
p=-0.5*lnrho(N,1)

i0=int(np.argmin(np.abs(N)))
rho0=rho[i0]
H0=H[i0]
Srel=np.sqrt(rho0/rho)
a_scale=np.exp(N)
Kmrel=Srel/a_scale
Kprel=Kmrel
E=H/H0
nu=Srel*E*p  # N_lapse /(tP*S0*H0)

# Absolute lapse uses the manuscript working present stretch S0.
S0_work=2.72e61
tP_s=5.391247e-44
MPC_KM=3.085677581491367e19
H0_kms=69.45160505326464
H0_s=H0_kms/MPC_KM
pref=tP_s*S0_work*H0_s
Nlapse=pref*nu

# Effective equation of state inferred from rho scaling if separately conserved.
w_eff=2.*p/3.-1.

def at_z(zt):
    j=int(np.argmin(np.abs(z-zt)))
    return {
      "z":float(z[j]),"Nlog":float(N[j]),"rho_smg":float(rho[j]),
      "p_struct":float(p[j]),"w_from_scaling":float(w_eff[j]),
      "S_over_S0":float(Srel[j]),"Km_over_Km0":float(Kmrel[j]),
      "Kp_over_Kp0":float(Kprel[j]),"E":float(E[j]),
      "nu_lapse":float(nu[j]),"N_lapse_S0_2p72e61":float(Nlapse[j])
    }

targets=[1e7,1e6,1e5,1e4,3400,1000,300,100,30,10,3,1,0.5,0.1,0.0]
rows=[at_z(x) for x in targets]

def stats(zlo,zhi):
    m=(z>=zlo)&(z<=zhi)&np.isfinite(p)
    if not np.any(m): return None
    return {
      "z_window":[zlo,zhi],"n":int(np.sum(m)),
      "p_median":float(np.median(p[m])),
      "p_mean":float(np.mean(p[m])),
      "p_min":float(np.min(p[m])),
      "p_max":float(np.max(p[m])),
      "N_lapse_median":float(np.median(Nlapse[m])),
      "Km_log_slope_median":float(np.median(p[m]-1.))
    }

out={
 "rho_positive_fraction":frac_positive,
 "lapse_prefactor_tP_S0_H0":pref,
 "today":at_z(0.0),
 "radiation_window":stats(1e5,1e7),
 "pre_equality_window":stats(5e3,5e4),
 "matter_window":stats(30,300),
 "transition_window":stats(1,5),
 "late_window":stats(0,0.2),
 "samples":rows
}
Path("output/late266_structural_composite_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("STRUCTURAL_COMPOSITE_AUDIT")
print("rho_positive_fraction",frac_positive)
print("lapse_prefactor_tP_S0_H0",pref)
for key in ["radiation_window","pre_equality_window","matter_window","transition_window","late_window"]:
    print(key,json.dumps(out[key],sort_keys=True))
print("TODAY",json.dumps(out["today"],sort_keys=True))
for row in rows:
    print("STRUCTURAL_SAMPLE",json.dumps(row,sort_keys=True))
