#!/usr/bin/env python3
"""
Mature inverse-square fixed-point basin audit.

For the scalar-dominated asymptotic normal form

    G2 = sigma^-2 f(Z),  G3 -> 0,  G4 -> F/2,

with Z=dot(sigma)^2/2 and constant F, the Friedmann equation gives

    H*sigma = sqrt[(2 Z f_Z - f)/(3 F)].

Using ln(sigma) as time, the exact homogeneous scalar equation becomes

    dZ/dlnsigma = -B/K,

    K = f_Z + 2 Z f_ZZ,

    B = 3(H sigma) f_Z sqrt(2Z) + 2f - 4Z f_Z.

At the coasting fixed point

    f + Z f_Z = 0,
    f_Z = 2F,

one finds universally

    d(delta Z)/dlnsigma = -2 delta Z,

provided K_* is finite and positive. Hence delta Z ~ sigma^-2 and
delta N/N ~ sigma^-2.

This script checks the exact nonlinear basin for the quadratic family

    f/(F Z*) = 2(1-r)y + r y^2 - (4-r),  y=Z/Z*,

for representative r values, including the value that reproduces the present
accepted scalar sound speed as an asymptotic reference.
"""
from pathlib import Path
import json, math
import numpy as np
from scipy.integrate import solve_ivp

OUT=Path("output/mature_fixed_point_basin.json")
CS2_PRESENT=0.138335
R_PRESENT_MATCH=(1./CS2_PRESENT-1.)/2.

RVALUES=[1.5,R_PRESENT_MATCH,3.9]
Y0S=[0.70,0.80,0.90,1.10,1.50,2.00]
XMAX=3.0

def fbar(y,r):
    return 2.*(1.-r)*y+r*y*y-(4.-r)

def fy(y,r):
    return 2.*(1.-r)+2.*r*y

def fyy(y,r):
    return 2.*r

def rho_bar(y,r):
    return 2.*y*fy(y,r)-fbar(y,r)

def Kbar(y,r):
    return fy(y,r)+2.*y*fyy(y,r)

def cs2(y,r):
    return fy(y,r)/Kbar(y,r)

def rhs(x,yy,r):
    y=float(yy[0])
    rb=rho_bar(y,r)
    kb=Kbar(y,r)
    if y<=0. or rb<=0. or kb==0.:
        return [np.nan]
    B=(math.sqrt(6.*y*rb)*fy(y,r)
       +2.*fbar(y,r)-4.*y*fy(y,r))
    return [-B/kb]

def run(r,y0):
    sol=solve_ivp(lambda x,y:rhs(x,y,r),(0.,XMAX),[y0],
                  rtol=1e-11,atol=1e-13,dense_output=True,max_step=.02)
    xs=np.array([0.,.25,.5,1.,2.,3.])
    ys=sol.sol(xs)[0] if sol.success else np.full_like(xs,np.nan)
    grad_floor=(r-1.)/r if r>1. else 0.
    kinetic_floor=(r-1.)/(3.*r) if r>1. else 0.
    return {
      "r":r,"y0":y0,"success":bool(sol.success),
      "gradient_stability_floor_y":grad_floor,
      "homogeneous_kinetic_floor_y":kinetic_floor,
      "initial_cs2":cs2(y0,r),
      "final_y":float(sol.y[0,-1]),
      "final_delta_y":float(sol.y[0,-1]-1.),
      "samples":[{"lnsigma":float(x),"sigma":float(math.exp(x)),
                  "y":float(y),"delta_y":float(y-1.)}
                 for x,y in zip(xs,ys)]
    }

runs=[run(r,y0) for r in RVALUES for y0 in Y0S]
out={
  "status":"exact homogeneous scalar-dominated normal-form basin; not yet the full stitched late266-to-future Horndeski evolution",
  "universal_linear_result":{
    "d_deltaZ_dlnsigma":"-2 deltaZ",
    "deltaZ_scaling":"sigma^-2",
    "deltaN_over_N_scaling":"sigma^-2"
  },
  "quadratic_family":{
    "r_values":RVALUES,
    "r_from_present_cs2_reference":R_PRESENT_MATCH,
    "cs2_fixed_point_formula":"1/(1+2r)"
  },
  "runs":runs
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("MATURE_FIXED_POINT_BASIN")
print("R_PRESENT_MATCH",R_PRESENT_MATCH)
for rr in runs:
    print("BASIN",json.dumps({
      "r":rr["r"],"y0":rr["y0"],"success":rr["success"],
      "initial_cs2":rr["initial_cs2"],"final_y":rr["final_y"],
      "gradient_floor":rr["gradient_stability_floor_y"]
    },sort_keys=True))
