#!/usr/bin/env python3
"""
Future scalar-curvature perturbation envelope audit for the C3-refined SDMC
tails.

This is an action-level mode audit, not a replacement for the full hi_class
multi-species hierarchy.  It uses the same released C3 No-Slip-refined action
as the homogeneous future audit and evaluates

    Q_s = 2 F D / (2-alpha_B)^2

together with the exact hi_class/Bellini-Sawicki background scalar sound speed.

The curvature mode equation is

    zeta_NN + [3 + dlnH/dN + dlnQ_s/dN] zeta_N
             + c_s^2 [k/(aH)]^2 zeta = 0.

At the mature inverse-square fixed point,

    alpha_M -> alpha_B -> 0,
    D -> 2(1+2r),
    c_s^2 -> 1/(1+2r),
    Q_s -> F(1+2r),
    a proportional to t,

so aH is constant and

    zeta ~ a^s,
    s = -1 +/- sqrt(1-c_s^2 [k/(aH)]^2).

Thus the mature super-horizon pair is a constant mode plus a^-2, while
sufficiently sub-horizon modes oscillate with envelope a^-1.  No mature
growing scalar-curvature mode is present.
"""
from pathlib import Path
import importlib.util, json, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE=Path(__file__).resolve().parent
BASE=HERE/"run_future_homogeneous_covariant_stitch_audit.py"
OUT=Path("output/future_perturbation_envelope_audit.json")

spec=importlib.util.spec_from_file_location("future_stitch",BASE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

NMAX=6.0
NPTS=1201
NU0=[0.0,0.1,1.0,10.0]

def background(model):
    sol=solve_ivp(lambda ne,y:m.rhs_diag(ne,y,model)[0],
                  (0.,NMAX),[1.,m.v0],
                  rtol=3e-9,atol=[1e-11,1e-14],
                  max_step=.01,dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)

    NN=np.linspace(0.,NMAX,NPTS)
    yy=sol.sol(NN)
    diag=[m.rhs_diag(float(ne),[float(s),float(v)],model)[1]
          for ne,(s,v) in zip(NN,yy.T)]
    H=np.asarray([d["H"] for d in diag])
    F=np.asarray([d["F"] for d in diag])
    aB=np.asarray([d["alphaB"] for d in diag])
    aM=np.asarray([d["alphaM"] for d in diag])
    D=np.asarray([d["D_full"] for d in diag])
    aBN=np.asarray([d["alphaB_N"] for d in diag])

    cs2=np.empty_like(NN)
    for j,d in enumerate(diag):
        b=aB[j]; mm=aM[j]; Hn=H[j]; Fv=F[j]
        xp=d["rho_smg_class"]+d["p_smg_class"]
        matter_class=(d["rho_m_action"]+4.*d["rho_r_action"]/3.)/3.
        num=((2.-b)*(b+2.*mm)/2.
             +1.5*(2.-b)*xp/(Hn*Hn)
             -1.5*(2.-2.*Fv+b*Fv)*matter_class/(Hn*Hn*Fv)
             +aBN[j])
        cs2[j]=num/D[j]

    Qs=2.*F*D/(2.-aB)**2
    h=CubicSpline(NN,np.log(H))(NN,1)
    dlnQs=CubicSpline(NN,np.log(Qs))(NN,1)
    friction=3.+h+dlnQs

    return NN,H,cs2,Qs,friction,D,aB,aM

def mode_audit(NN,H,cs2,friction,nu0):
    if nu0==0.:
        return dict(nu0=0.,max_abs_zeta=1.,final_zeta=1.,
                    final_abs_zeta=1.,max_growth_over_initial=1.,
                    success=True)

    sH=CubicSpline(NN,np.log(H))
    sc=CubicSpline(NN,cs2)
    sf=CubicSpline(NN,friction)
    H0=H[0]

    def rhs(ne,y):
        Hn=math.exp(float(sH(ne)))
        c2=float(sc(ne))
        fr=float(sf(ne))
        nu=nu0*H0/(math.exp(ne)*Hn)
        return [y[1],-fr*y[1]-c2*nu*nu*y[0]]

    sol=solve_ivp(rhs,(NN[0],NN[-1]),[1.,0.],
                  rtol=2e-8,atol=2e-10,max_step=.01)
    z=np.asarray(sol.y[0])
    return dict(
      nu0=nu0,
      max_abs_zeta=float(np.max(np.abs(z))),
      max_growth_over_initial=float(np.max(np.abs(z))),
      final_zeta=float(z[-1]),
      final_abs_zeta=float(abs(z[-1])),
      success=bool(sol.success)
    )

out={
  "status":(
    "future scalar-curvature envelope audit on the C3 No-Slip-refined "
    "homogeneous background. It tests the action-level scalar mode envelope "
    "through the transition but is not the complete multi-species hi_class "
    "perturbation hierarchy beyond a=1."
  ),
  "analytic_mature_result":{
    "Qs":"F_inf*(1+2r)",
    "mode_equation":"zeta_NN+2 zeta_N+cs2*(k/(aH))^2 zeta=0",
    "mode_exponents":"s=-1 +/- sqrt(1-cs2*(k/(aH))^2)",
    "superhorizon":"constant mode plus a^-2",
    "subhorizon":"oscillatory with a^-1 envelope",
    "growing_mode":False
  },
  "candidates":{}
}

for cspec in m.CANDIDATES:
    model=m.build_refined_noslip_candidate(
      cspec,iterations=2,blend_rate=500.
    )
    NN,H,cs2,Qs,friction,D,aB,aM=background(model)
    modes=[mode_audit(NN,H,cs2,friction,x) for x in NU0]
    out["candidates"][cspec["name"]]={
      "N_inf":model["N_inf"],
      "min_Qs":float(np.min(Qs)),
      "max_Qs":float(np.max(Qs)),
      "min_friction":float(np.min(friction)),
      "max_friction":float(np.max(friction)),
      "min_D":float(np.min(D)),
      "min_cs2":float(np.min(cs2)),
      "max_cs2":float(np.max(cs2)),
      "max_abs_alphaB_plus_2alphaM":
        float(np.max(np.abs(aB+2.*aM))),
      "modes":modes
    }

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_PERTURBATION_ENVELOPE_AUDIT")
for name,r in out["candidates"].items():
    print("PERTURB",name,json.dumps(r,sort_keys=True))
