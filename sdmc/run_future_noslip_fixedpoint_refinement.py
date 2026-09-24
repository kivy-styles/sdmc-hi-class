#!/usr/bin/env python3
"""
C2 No-Slip fixed-point refinement for the released future SDMC tails.

The baseline future action already preserves the accepted z>=0 action and
approaches the inverse-square mature fixed point, but its independently
stitched linear-G3 coefficient leaves a small transient
alpha_B+2 alpha_M ~ 5e-5.

This script iteratively reconstructs the future-only g(sigma) profile from

    g_required(sigma) = -F_,sigma / (2 Z_native(sigma)),

using the native homogeneous trajectory from the previous iteration.  The
correction is multiplied by the C2 switch

    W(x)=1-exp(-lambda x)(1+lambda x+(lambda x)^2/2), x=ln sigma,

so the accepted present value and first two derivatives of g are left
unchanged exactly.  The default lambda=100 turns on the correction rapidly
but smoothly in the unobserved future.

The refinement is an inverse action reconstruction.  The final profile is
validated by a fresh native homogeneous evolution and by the exact Horndeski
D and c_s^2 audit from run_future_homogeneous_covariant_stitch_audit.py.
"""
from pathlib import Path
import json, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/future_noslip_fixedpoint_refinement.json")
LAMBDA_SWITCH=100.0
NITER=4
NMAX=10.0
GRID=6001

def trajectory(model):
    sol=solve_ivp(
      lambda Ne,y: aud.rhs_diag(Ne,y,model)[0],
      (0.,NMAX),[1.,aud.v0],
      rtol=2e-9,atol=[1e-11,1e-14],
      max_step=.004,dense_output=True
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    NN=np.linspace(0.,NMAX,GRID)
    yy=sol.sol(NN)
    diag=[
      aud.rhs_diag(float(ne),[float(sig),float(v)],model)[1]
      for ne,(sig,v) in zip(NN,yy.T)
    ]
    return NN,yy,diag

def refine(spec):
    base=aud.build_candidate(spec)
    model=base
    history=[]

    for iteration in range(NITER+1):
        NN,yy,diag=trajectory(model)
        ns=np.array([x["noslip_alpha_combo"] for x in diag])
        qq=np.array([x["q"] for x in diag])
        nstruct=np.array([x["N_struct"] for x in diag])

        imax=int(np.argmax(np.abs(ns)))
        history.append({
          "iteration":iteration,
          "max_abs_alphaB_plus_2alphaM":float(abs(ns[imax])),
          "at_ln_a":float(NN[imax]),
          "N_at_final":float(nstruct[-1]),
          "q_max":float(np.max(qq)),
        })

        if iteration==NITER:
            break

        sigarr=yy[0]
        varr=yy[1]
        xarr=np.log(sigarr)
        zarr=.5*varr*varr

        keep=np.r_[True,np.diff(xarr)>1e-12]
        xarr=xarr[keep]
        sigarr=sigarr[keep]
        zarr=zarr[keep]

        greq=[]
        for sig,z in zip(sigarr,zarr):
            aq=base["action"](float(sig))
            greq.append(-aq["Fs"]/(2.*z))
        greq=np.asarray(greq)
        sp=CubicSpline(xarr,greq,bc_type="natural")

        base_action=base["action"]
        lam=LAMBDA_SWITCH
        xlo=float(xarr[0]); xhi=float(xarr[-1])

        def action_refined(sig,sp=sp,base_action=base_action,
                           lam=lam,xlo=xlo,xhi=xhi):
            out=base_action(sig)
            x=math.log(sig)
            xc=min(max(x,xlo),xhi)

            # C2 switch: W(0)=W'(0)=W''(0)=0, W(infinity)=1.
            u=lam*x
            eu=math.exp(-u)
            W=1.-eu*(1.+u+.5*u*u)
            Wx=.5*lam*eu*u*u
            Wxx=.5*lam*lam*eu*(2.*u-u*u)

            gr=float(sp(xc))
            grx=float(sp(xc,1))
            grxx=float(sp(xc,2))

            g0=out["g"]
            g0x=out["gs"]*sig
            g0xx=out["gss"]*sig*sig+g0x

            dg=gr-g0
            dgx=grx-g0x
            dgxx=grxx-g0xx

            gv=g0+W*dg
            gx=g0x+Wx*dg+W*dgx
            gxx=g0xx+Wxx*dg+2.*Wx*dgx+W*dgxx

            out["g"]=gv
            out["gs"]=gx/sig
            out["gss"]=(gxx-gx)/(sig*sig)
            return out

        model={**base,"action":action_refined}

    # Fresh end-to-end validation of the refined fixed action.
    validation=aud.evolve(model,Nmax=NMAX,npts=2001)

    return {
      "lambda_switch":LAMBDA_SWITCH,
      "iterations":NITER,
      "history":history,
      "validation":validation,
    }

results={}
for spec in aud.CANDIDATES:
    results[spec["name"]]=refine(spec)

out={
  "status":(
    "future-only inverse reconstruction of the linear-G3 coefficient; "
    "the accepted z>=0 action is unchanged through second derivative at "
    "sigma=1.  The refined profile is validated by a fresh native homogeneous "
    "evolution, but full future perturbation propagation inside hi_class "
    "remains outstanding."
  ),
  "results":results,
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_NOSLIP_FIXEDPOINT_REFINEMENT")
for name,r in results.items():
    v=r["validation"]
    print("REFINED",name,json.dumps({
      "history":r["history"],
      "health":v["health_diagnostics"],
      "future":v["future_kinematics"],
      "final":v["final"],
    },sort_keys=True))
