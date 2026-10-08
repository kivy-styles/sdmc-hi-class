#!/usr/bin/env python3
"""SDMC kinetic D0 *correction* to Dfloor: conditional derivation and action degeneracy.

Strict scope: the accepted linear-G3 covariant reconstruction
  G2=k1(phi)X+k2(phi)X^2-V(phi), G3=g(phi)X, G4=F(phi)/2.
This proves algebraic properties on an arbitrary known homogeneous trajectory;
it does not independently construct a microscopic action or compute a new
cosmological likelihood. It adds no free parameter to the accepted benchmark.
"""
import json
import math
from pathlib import Path

p=json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]
D_floor=p["D_floor"]
D0=p["D0"]
power=p["kinetic_power"]
wD=p["width"];zD=p["z_c"]
# Derived-bundle F differs from the historic D-window: preserved by its
# current successful covariant replay workflow.
AF=0.020520;zcF=4.034077502899133;wF=0.32925377073762596

def S(N,z_c,w):
    arg=(N+math.log1p(z_c))/w
    return 1/(1+math.exp(-arg))

def dS_dN(N,z_c,w):
    ss=S(N,z_c,w)
    return ss*(1-ss)/w

def F(N):
    return math.exp(AF*S(N,zcF,wF))

def Xbg(N):
    # Illustrative strictly positive homogeneous X(φ) profile, not a
    # claimed physical SDMC background or an independent H prediction.
    return 0.5*(0.70+0.30*math.exp(-3*N))

def Xbg_prime(N):
    return -0.45*math.exp(-3*N)

def G2delta(N,X,mu):
    dx=X-Xbg(N)
    return mu*dx*dx

def G2delta_X(N,X,mu):
    return 2*mu*(X-Xbg(N))

def G2delta_phi(N,X,mu,mu_prime):
    dx=X-Xbg(N)
    return mu_prime*dx*dx-2*mu*dx*Xbg_prime(N)

def need(ok,msg):
    if not ok: raise AssertionError(msg)

def close(a,b,tol=2e-12):
    return math.isclose(a,b,rel_tol=tol,abs_tol=tol)

def main():
    ncenter=-math.log1p(zD)
    S0=S(0,zD,wD)
    correct_today=D0*S0**power
    D_today=D_floor+correct_today
    Dlate=D_floor+D0
    slope_mid=D0*power*(0.5**(power-1))*dS_dN(ncenter,zD,wD)
    D0_deriv_mid=slope_mid*(2**(power+1))*wD/power
    need(close(D0_deriv_mid,D0),"D0 recovered from derivative at transition centre")
    need(close(Dlate-D_floor,D0),"D0 recovered from true late-minus-early support")
    need(close(D_today,0.3913854938989103),"accepted D(today) check")
    # This local derivative is a conditional extraction formula only: its
    # value is computed FROM the accepted correction, not predicted.
    results={"status":"verified correction identities and background-preserving covariant G2 degeneracy",
      "historical_corrective_interpretation":{
        "Dfloor":D_floor,"D0":D0,"kinetic_power":power,"D0_is_correction_to_Dfloor":True,
        "S_D_today":S0,"D_correction_today":correct_today,"D_today":D_today,
        "D_late_asymptote":Dlate,
        "D0_from_midpoint_slope_if_the_slope_were_independently_given":D0_deriv_mid,
        "midpoint_slope_dD_dN_not_independently_predicted":slope_mid,
        "midpoint_N":ncenter,
        "conditional_identity":"D0=D_late-Dfloor if S(+infty)=1",
      },
      "background_preserving_deformation":{
        "definition":"delta G2 = mu(phi)*(X-Xbg(phi))^2",
        "on_shell_trajectory":"phi=N, Xbg(phi)>0",
        "background_vanishing":["deltaG2","deltaG2_X","deltaG2_phi"],
        "perturbation_nonzero":"deltaG2_XX=2mu; deltaalphaK=deltaD=4*mu*Xbg/F in current linear-G3 conventions",
        "recipe":"mu(phi)=F(phi)*deltaD0*S_D(phi)^pD/(4*Xbg(phi))",
        "implication":"For any modest deltaD0, this changes D0 by deltaD0 while preserving G2, G2_X, G2_phi on the selected background. It changes the action and linear perturbations; it does not establish observational viability for each deformation.",
      },
      "cases":[],
      "checks":[],
    }
    # Different allowed corrections preserve the same chosen background
    # G2 and G2_X on-shell; their kinetic second derivative shifts.
    for deltaD0 in [-0.04,0.,0.04]:
      case={"deltaD0":deltaD0,"D0_result":D0+deltaD0,"nodes":[]}
      for z in [0.,.5,1.,4.,16.74228966,100.]:
        N=-math.log1p(z);xb=Xbg(N); f=F(N); ss=S(N,zD,wD)
        target=deltaD0*ss**power
        mu=f*target/(4*xb)
        # Same potential shifts: delta k1=-2mu Xbg, delta k2=mu,
        # delta V=-mu Xbg^2. delta G2(Xbg)=0, G2_X(Xbg)=0.
        dk1=-2*mu*xb
        dk2=mu
        dV=-mu*xb*xb
        same_action_value=dk1*xb+dk2*xb*xb-dV
        same_1st_X=dk1+2*dk2*xb
        shifted_kinetic=(dk1+6*dk2*xb)/f
        need(close(same_action_value,0,1e-9),"background G2 is unchanged")
        need(close(same_1st_X,0,1e-9),"background G2_X unchanged")
        need(close(G2delta(N,xb,mu),0) and close(G2delta_X(N,xb,mu),0),
             "background path is a double root")
        # phi derivative also zero on background even if mu'(phi)≠0.
        need(close(G2delta_phi(N,xb,mu,1.5),0),"background G2_phi unchanged")
        need(close(shifted_kinetic,target,1e-12),"kinetic correction match")
        case["nodes"].append({"z":z,"S_D":ss,"F":f,"X_bg_illustrative":xb,
                              "delta_D_at_node":shifted_kinetic,
                              "mu_illustrative":mu,
                              "delta_G2_X_on_bg":same_1st_X,
                              "delta_G2_on_bg":same_action_value})
      results["cases"].append(case)
    results["checks"]=[
      "verified accepted D0 correction and D0=Dlate-Dfloor",
      "verified D0 from midpoint slope D'(Nc)*2^(p+1)*w/p",
      "verified G2 and G2_X unchanged on chosen homogeneous background for 3 different corrections",
      "verified G2_phi vanishes on chosen background",
      "verified deltaD=4mu Xbg/F reproduces prescribed correction at 6 redshifts",
    ]
    results["physical_limits"]=[
      "The whole selected background is unchanged, but the deformed action differs at quadratic order and alters the scalar sound speed and perturbed observables",
      "This is a mathematical non-uniqueness construction, NOT a new preferred SDMC action or a free-evolution test of the alternatives",
      "A prediction for D0 requires independent UV coefficients, scalar initial conditions, a predicted late kinetic endpoint, or an additional physically justified selection principle",
      "Positivity of gradient sound speed alone cannot fix D0 when it enters the strictly positive kinetic denominator at fixed background",
    ]
    out=Path("output/late300_D0_correction_origin")
    out.mkdir(parents=True,exist_ok=True)
    (out/"audit.json").write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
    print("SDMC_D0_CORRECTION_DEGENERACY_SUCCESS",
          json.dumps({"D_today":D_today,"D_late":Dlate,
                      "midpoint_derivative":slope_mid,
                      "midpoint_D0_recovered":D0_deriv_mid,
                      "variants":[v["D0_result"] for v in results["cases"]],
                      "tests":len(results["checks"]),
                      "conclusion":"correction status verified, value underdetermined by background/no-slip equations"},sort_keys=True))

if __name__=="__main__":main()
