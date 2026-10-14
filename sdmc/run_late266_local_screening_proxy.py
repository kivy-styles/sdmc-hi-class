#!/usr/bin/env python3
from pathlib import Path
import math, re

p=Path("output/linear_cov_reconstruction_summary.txt")
d={}
for line in p.read_text().splitlines():
    if "=" not in line or line.strip().endswith("]"): continue
    k,v=line.split("=",1)
    try: d[k.strip()]=float(v.strip())
    except ValueError: pass

H0_kms=d["H0_kms"]
F0=d["F0_target"]
alphaM=d["alphaM0_target"]
alphaB=d["alphaB0_target"]
D0=d["D0_target"]
ghat=abs(d["g0_H0sq"])

# From the scalar quadratic action used by hi_class:
# Q_s/Mp^2 = 2 F D/(2-alpha_B)^2.
# With phi ~= ln a, delta phi = -zeta on the background trajectory.
# Canonical normalization in the time-kinetic sector therefore uses
# Z_t = 2 Q_s/Mp^2.
Zt=4.*F0*D0/(2.-alphaB)**2

# Read the free-evolution present c_s^2 if the comparison file exists.
cs2=0.1383349697896
j=Path("output/linear_cov_replay_comparison.json")
if j.exists():
    import json
    cs2=float(json.loads(j.read_text())["cs2_0"])
Zr=Zt*cs2

# Einstein-frame conformal coupling proxies after canonical normalization.
Qt=abs(alphaM)/(2.*math.sqrt(Zt))
Qr=abs(alphaM)/(2.*math.sqrt(Zr))

# Linear G3 coefficient satisfies g H^2 = -F alpha_M on the target.
# For a canonical field, Lambda^3/(Mp H0^2) = Z^(3/2)/|g H0^2|.
Lam3_t=Zt**1.5/ghat
Lam3_r=Zr**1.5/ghat

c_kms=299792.458
H0=H0_kms/c_kms   # 1/Mpc in c=1 units
MPC_M=3.085677581491367e22
AU_M=149597870700.0
PC_AU=206264.80624709636
YEAR_S=365.25*86400.
MPC_KM=3.085677581491367e19

# Schwarzschild radii.
rs_sun_m=2953.25008
rs_earth_m=0.008870056

def rV_pc(rs_m,Q,lam3ratio):
    rs=rs_m/MPC_M
    # Convention: rV^3 = 2 Q r_s/(Lambda^3/Mp).
    # Lambda^3/Mp = lam3ratio*H0^2.
    return (2.*Q*rs/(lam3ratio*H0*H0))**(1./3.)*1e6

rv_sun_t=rV_pc(rs_sun_m,Qt,Lam3_t)
rv_sun_r=rV_pc(rs_sun_m,Qr,Lam3_r)
rv_earth_t=rV_pc(rs_earth_m,Qt,Lam3_t)
rv_earth_r=rV_pc(rs_earth_m,Qr,Lam3_r)

def deep_screened_eps(rs_m,r_m):
    # In the deep-Vainshtein cubic limit, combining
    # Q=alphaM/(2 sqrt Z), Lambda^3=Mp H0^2 Z^(3/2)/(F alphaM)
    # cancels Z:
    # eps_phi ~= 2Q^2 (r/rV)^(3/2)
    #          = alphaM H0 r^(3/2)/(2 sqrt(F r_s)).
    rs=rs_m/MPC_M
    r=r_m/MPC_M
    return abs(alphaM)*H0*r**1.5/(2.*math.sqrt(F0*rs))

eps_sun_1au=deep_screened_eps(rs_sun_m,AU_M)
eps_earth_moon=deep_screened_eps(rs_earth_m,384400e3)
eps_earth_surface=deep_screened_eps(rs_earth_m,6371e3)

# Unscreened scalar-force and simple conformal PPN proxies.
force_uns_t=2.*Qt*Qt
force_uns_r=2.*Qr*Qr
gamma_minus1_t=-4.*Qt*Qt
gamma_minus1_r=-4.*Qr*Qr

# Naive unscreened cosmological Gdot/G = -alpha_M H0.
H0_yr=(H0_kms/MPC_KM)*YEAR_S
gdot=-alphaM*H0_yr

lines=[
"STATUS=decoupling-limit proxy; not an exact PPN solution",
f"F0={F0:.15g}",
f"alphaM0={alphaM:.15g}",
f"alphaB0={alphaB:.15g}",
f"D0={D0:.15g}",
f"cs2_0={cs2:.15g}",
f"gH0sq={ghat:.15g}",
f"Z_time={Zt:.15g}",
f"Z_radial_proxy={Zr:.15g}",
f"Q_time={Qt:.15g}",
f"Q_radial_proxy={Qr:.15g}",
f"Lambda3_over_MpH0sq_time={Lam3_t:.15g}",
f"Lambda3_over_MpH0sq_radial={Lam3_r:.15g}",
f"rV_sun_pc_time={rv_sun_t:.15g}",
f"rV_sun_pc_radial={rv_sun_r:.15g}",
f"rV_earth_pc_time={rv_earth_t:.15g}",
f"rV_earth_pc_radial={rv_earth_r:.15g}",
f"unscreened_force_fraction_time={force_uns_t:.15g}",
f"unscreened_force_fraction_radial={force_uns_r:.15g}",
f"gamma_minus1_proxy_time={gamma_minus1_t:.15g}",
f"gamma_minus1_proxy_radial={gamma_minus1_r:.15g}",
f"deep_screened_force_sun_1AU={eps_sun_1au:.15g}",
f"deep_screened_force_earth_moon={eps_earth_moon:.15g}",
f"deep_screened_force_earth_surface={eps_earth_surface:.15g}",
f"naive_Gdot_over_G_per_yr={gdot:.15g}",
]
out=Path("output/late266_local_screening_proxy.txt")
out.write_text("\n".join(lines)+"\n")
print("\n".join(lines))
