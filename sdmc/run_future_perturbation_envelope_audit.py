#!/usr/bin/env python3
"""
Future scalar-perturbation envelope audit for the released SDMC tails.

This imports the exact homogeneous future stitch audit, reconstructs the full
Bellini-Sawicki scalar kinetic coefficient

    Q_s = 2 F D / (2-alpha_B)^2,

and integrates the curvature mode equation

    zeta_NN + [3 + dlnH/dN + dlnQ_s/dN] zeta_N
             + c_s^2 [k/(aH)]^2 zeta = 0

for representative present-horizon ratios nu0=k/H0.

At the mature inverse-square fixed point:
    alpha_M -> alpha_B -> 0,
    D -> 2(1+2r),
    c_s^2 -> 1/(1+2r),
    Q_s -> F(1+2r),
    a proportional to t.

Hence k/(aH) is constant and
    zeta ~ a^s,
    s = -1 +/- sqrt(1 - c_s^2 [k/(aH)]^2).

There is therefore no growing mature scalar mode: the super-horizon pair is a
constant mode plus a^-2, while sub-horizon modes oscillate with envelope a^-1.
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

    aM=[]; aB=[]; aK=[]; D=[]; matter=[]
    for ne,(sig,v),dd in zip(NN,yy.T,diag):
        aq=model["action"](float(sig))
        Hn=dd["H"]; Z=.5*v*v; F=aq["F"]
        am=v*aq["Fs"]/(Hn*F)
        ab=2.*v*(Z*aq["g"]-.5*aq["Fs"])/(Hn*F)
        ak=(2.*Z*(aq["k1"]+6.*aq["k2"]*Z-4.*Z*aq["gs"])
            +12.*v*Z*Hn*aq["g"])/(Hn*Hn*F)
        ddyn=ak+1.5*ab*ab
        rm=m.rho_m0*math.exp(-3.*ne)
        rr=m.rho_r0*math.exp(-4.*ne)
        aM.append(am); aB.append(ab); aK.append(ak); D.append(ddyn)
        matter.append((rm+4.*rr/3.)/(Hn*Hn*F))

    aM=np.asarray(aM); aB=np.asarray(aB); aK=np.asarray(aK)
    D=np.asarray(D); matter=np.asarray(matter)
    aBp=CubicSpline(NN,aB)(NN,1)
    h=CubicSpline(NN,np.log(H))(NN,1)
    cs2=((2.-aB)*(-h+.5*aB+aM)-matter+aBp)/D

    F=np.asarray([model["action"](float(sig))["F"] for sig in yy[0]])
    Qs=2.*F*D/(2.-aB)**2
    dlnQs=CubicSpline(NN,np.log(Qs))(NN,1)
    friction=3.+h+dlnQs

    return NN,H,cs2,Qs,friction,D,aB,aM

def mode_audit(NN,H,cs2,Qs,friction,nu0):
    if nu0==0.:
        return dict(nu0=0.,max_abs_zeta=1.,final_zeta=1.,
                    final_abs_zeta=1.,success=True)

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
    z=sol.y[0]
    return dict(
      nu0=nu0,
      max_abs_zeta=float(np.max(np.abs(z))),
      final_zeta=float(z[-1]),
      final_abs_zeta=float(abs(z[-1])),
      success=bool(sol.success)
    )

out={
  "status":(
    "future scalar-curvature envelope audit on the released homogeneous "
    "background. It tests linear scalar-mode damping through the transition "
    "but is not a replacement for the complete multi-species hi_class "
    "perturbation hierarchy beyond a=1."
  ),
  "analytic_mature_result":{
    "Qs":"F_inf*(1+2r)",
    "mode_exponents":"s=-1 +/- sqrt(1-nu^2), nu^2=cs2*(k/(aH))^2",
    "superhorizon":"constant mode plus a^-2",
    "subhorizon":"oscillatory with a^-1 envelope",
    "growing_mode":False
  },
  "candidates":{}
}

for cspec in m.CANDIDATES:
    model=m.build_candidate(cspec)
    NN,H,cs2,Qs,friction,D,aB,aM=background(model)
    modes=[mode_audit(NN,H,cs2,Qs,friction,x) for x in NU0]
    out["candidates"][cspec["name"]]={
      "N_inf":model["N_inf"],
      "min_Qs":float(np.min(Qs)),
      "max_Qs":float(np.max(Qs)),
      "min_friction":float(np.min(friction)),
      "max_friction":float(np.max(friction)),
      "min_D":float(np.min(D)),
      "min_cs2":float(np.min(cs2)),
      "max_cs2":float(np.max(cs2)),
      "modes":modes
    }

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_PERTURBATION_ENVELOPE_AUDIT")
for name,r in out["candidates"].items():
    print("PERTURB",name,json.dumps(r,sort_keys=True))
