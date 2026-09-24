#!/usr/bin/env python3
"""
Focused refined search around the split-rate legacy-canonical stitch.

The preceding split-rate audit found an unrefined strict point near
(mu_k1,mu_k2,mu_V)=(40,10,120), while the two-iteration C3 No-Slip refinement
missed the c_s^2<=1.0001 gate by only ~9e-6.  This script searches a narrow
future-only neighborhood after the full two-iteration refinement.

Accepted z>=0 coefficients are unchanged.  The mature target is fixed to
    r=0, chi_inf=F_inf, N_inf=sqrt(8*pi/3).
"""
from pathlib import Path
import json,math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/legacy_canonical_refined_local_search.json")
MU1=[36.0,40.0,44.0]
MU2=[8.0,10.0,12.0]
MUV=[110.0,120.0,130.0]
TOL=1e-4

def evaluate(mu1,mu2,muv):
    spec={
      "name":f"canon_refined_{mu1:g}_{mu2:g}_{muv:g}",
      "chi_inf":aud.FINF,"r":0.0,"muQ":mu1,
      "mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
    }
    try:
        model=aud.build_refined_noslip_candidate(
            spec,iterations=2,blend_rate=500.)
        res=aud.evolve(model,Nmax=10.,npts=2001)
        h=res["health_diagnostics"]
        gm=float(h["Gamma_action_max"])
        cmin=float(h["min_cs2_horndeski"])
        cmax=float(h["max_cs2_horndeski"])
        dmin=float(h["min_D_horndeski"])
        ns=float(h["max_abs_alphaB_plus_2alphaM"])
        strict=(gm<=1.+TOL and cmin>0. and cmax<=1.+TOL and dmin>0.)
        excess=max(0.,gm-(1.+TOL))+10.*max(0.,cmax-(1.+TOL))
        score=excess+ns+1e-3*abs(float(res["future_kinematics"]["q_max"]))
        return {
          "status":"ok","mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "Gamma_max":gm,
          "Gamma_final":float(h["Gamma_action_final"]),
          "min_D":dmin,"min_cs2":cmin,"max_cs2":cmax,
          "max_abs_alphaB_plus_2alphaM":ns,
          "N_final":float(res["final"]["N_struct"]),
          "p_final":float(res["final"]["p"]),
          "q_final":float(res["final"]["q"]),
          "q_max":float(res["future_kinematics"]["q_max"]),
          "strict_feasible":strict,"score":score,
        }
    except Exception as e:
        return {
          "status":"failure","mu_k1":mu1,"mu_k2":mu2,"mu_V":muv,
          "strict_feasible":False,"score":1e99,"error":str(e)
        }

rows=[]
for mu1 in MU1:
    for mu2 in MU2:
        for muv in MUV:
            rows.append(evaluate(mu1,mu2,muv))

ok=[r for r in rows if r["status"]=="ok"]
strict=[r for r in ok if r["strict_feasible"]]
best=min(strict,key=lambda r:r["score"]) if strict else min(ok,key=lambda r:r["score"])

# Curvature-envelope check only for the selected refined point.
spec={
  "name":"canon_refined_best","chi_inf":aud.FINF,"r":0.0,
  "muQ":best["mu_k1"],"mu_k1":best["mu_k1"],
  "mu_k2":best["mu_k2"],"mu_V":best["mu_V"],
}
model=aud.build_refined_noslip_candidate(spec,iterations=2,blend_rate=500.)

NMAX=6.; NPTS=1201
sol=solve_ivp(lambda ne,y:aud.rhs_diag(ne,y,model)[0],
              (0.,NMAX),[1.,aud.v0],
              rtol=3e-9,atol=[1e-11,1e-14],max_step=.01,dense_output=True)
NN=np.linspace(0.,NMAX,NPTS); yy=sol.sol(NN)
diag=[aud.rhs_diag(float(ne),[float(s),float(v)],model)[1]
      for ne,(s,v) in zip(NN,yy.T)]
H=np.asarray([d["H"] for d in diag])
F=np.asarray([d["F"] for d in diag])
aB=np.asarray([d["alphaB"] for d in diag])
aM=np.asarray([d["alphaM"] for d in diag])
D=np.asarray([d["D_full"] for d in diag])
aBN=np.asarray([d["alphaB_N"] for d in diag])
c2=np.empty_like(NN)
for j,d in enumerate(diag):
    b=aB[j]; mm=aM[j]; Hn=H[j]; Fv=F[j]
    xp=d["rho_smg_class"]+d["p_smg_class"]
    matter=(d["rho_m_action"]+4.*d["rho_r_action"]/3.)/3.
    num=((2.-b)*(b+2.*mm)/2.
         +1.5*(2.-b)*xp/(Hn*Hn)
         -1.5*(2.-2.*Fv+b*Fv)*matter/(Hn*Hn*Fv)
         +aBN[j])
    c2[j]=num/D[j]
Qs=2.*F*D/(2.-aB)**2
fr=3.+CubicSpline(NN,np.log(H))(NN,1)+CubicSpline(NN,np.log(Qs))(NN,1)

def mode(nu0):
    if nu0==0.:
        return {"nu0":0.,"success":True,"max_abs_zeta":1.,"final_abs_zeta":1.}
    sH=CubicSpline(NN,np.log(H)); sc=CubicSpline(NN,c2); sf=CubicSpline(NN,fr)
    H0=H[0]
    def rhs(ne,y):
        hn=math.exp(float(sH(ne)))
        nu=nu0*H0/(math.exp(ne)*hn)
        return [y[1],-float(sf(ne))*y[1]-float(sc(ne))*nu*nu*y[0]]
    ss=solve_ivp(rhs,(NN[0],NN[-1]),[1.,0.],
                 rtol=2e-8,atol=2e-10,max_step=.01)
    z=np.asarray(ss.y[0])
    return {
      "nu0":nu0,"success":bool(ss.success),
      "max_abs_zeta":float(np.max(np.abs(z))),
      "final_abs_zeta":float(abs(z[-1]))
    }

modes=[mode(x) for x in [0.,.1,1.,10.]]
out={
  "status":(
    "focused two-iteration C3 No-Slip search around the split-rate canonical "
    "future stitch; accepted z>=0 action unchanged"
  ),
  "strict_limits":{"Gamma_max":1.+TOL,"cs2_max":1.+TOL,
                   "D_min":0.,"cs2_min":0.},
  "rows":rows,
  "strict_count":len(strict),
  "best":best,
  "best_modes":modes,
  "best_all_envelope_modes_non_growing":
    all(m["success"] and m["max_abs_zeta"]<=1.000001 for m in modes),
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("LEGACY_CANONICAL_REFINED_LOCAL_SEARCH")
print("STRICT_COUNT",len(strict))
print("BEST",json.dumps(best,sort_keys=True))
print("MODES",json.dumps(modes,sort_keys=True))
