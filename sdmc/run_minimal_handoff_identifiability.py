#!/usr/bin/env python3
from pathlib import Path
import json, math
import numpy as np
from scipy.optimize import least_squares, minimize_scalar

OUT=Path("output/minimal_handoff_identifiability"); OUT.mkdir(parents=True,exist_ok=True)

# Accepted present endpoint from the supported-domain canonical-approach / asymptotic-match audits.
Om0=0.2985060132835426
q0=-0.5272467393408948
qN0=-0.7430775900807528
sqrt6=math.sqrt(6.0)

# Manuscript Section 18:
# u=a^-1
# Omega_m = m1*u + 3*m1^2*u^2 + O(u^3 ln u)
# q = -(3/2)m1*u + u^2[(9/2)m1^2 ln u -(9/2)m1^2 -2sqrt(6)s2] + ...
def section18(u,m1,s2):
    if not (u>0): return (np.nan,np.nan,np.nan)
    Om=m1*u+3.0*m1*m1*u*u
    B=4.5*m1*m1*math.log(u)-4.5*m1*m1-2.0*sqrt6*s2
    q=-1.5*m1*u+u*u*B
    dqdu=-1.5*m1+2*u*B+4.5*m1*m1*u
    qN=-u*dqdu
    return Om,q,qN

def residual(x):
    u,m1,s2=x
    Om,q,qN=section18(u,m1,s2)
    return np.array([Om-Om0,q-q0,qN-qN0])

# Search for a direct shifted-clock embedding of today's three invariants.
best=None
for u0 in np.geomspace(0.03,0.95,12):
  for m0 in (0.05,0.1,0.2,0.4,0.8,1.5,3.0):
    for s0 in (-3,-1,-.3,0,.3,1,3):
      r=least_squares(residual,[u0,m0,s0],bounds=([1e-4,1e-8,-20],[0.999999,20,20]),
                      xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=20000)
      score=float(np.linalg.norm(r.fun))
      if best is None or score<best["score"]:
        best={"score":score,"x":r.x.tolist(),"residual":r.fun.tolist(),"success":bool(r.success)}

# Stronger one-dimensional diagnostic: for every u, force Omega_m and q to match
# exactly, then inspect the implied q_N.  Omega_m gives a positive quadratic m1.
def m1_from_Om(u):
    # 3 u^2 m1^2 + u m1 - Om0 = 0
    return (-u+math.sqrt(u*u+12*u*u*Om0))/(6*u*u)

def s2_from_q(u,m1):
    no_s=-1.5*m1*u+u*u*(4.5*m1*m1*math.log(u)-4.5*m1*m1)
    return (no_s-q0)/(2*sqrt6*u*u)

def qp_resid_abs(logu):
    u=math.exp(logu)
    m=m1_from_Om(u); s=s2_from_q(u,m)
    return abs(section18(u,m,s)[2]-qN0)

opt=minimize_scalar(qp_resid_abs,bounds=(math.log(1e-3),math.log(0.999999)),method="bounded",
                    options={"xatol":1e-13,"maxiter":1000})
ub=math.exp(opt.x); mb=m1_from_Om(ub); sb=s2_from_q(ub,mb)
omb,qb,qpb=section18(ub,mb,sb)

out={
 "status":"minimal future-handoff identifiability / shifted-mature-branch audit",
 "accepted_present":{"Omega_m0":Om0,"q0":q0,"dq_dN0":qN0},
 "manuscript_basis":{
   "section17_eigenmodes":{"matter":-1.0,"scalar":-2.0},
   "section18_variables":["u=a^-1","m1","s2"],
   "equations_used":["18.1","18.3","18.5"],
   "scope":"second-order mature asymptotic expansion only"
 },
 "direct_three_invariant_fit":{
   "variables":["u","m1","s2"],
   "best_u_m1_s2":best["x"],
   "residual_Omega_q_qprime":best["residual"],
   "euclidean_residual":best["score"],
   "exact_embedding_found":bool(best["score"]<1e-8)
 },
 "forced_Omega_q_match_scan":{
   "best_u":ub,"best_m1":mb,"best_s2":sb,
   "matched_Omega_m":omb,"matched_q":qb,
   "implied_dq_dN":qpb,
   "accepted_dq_dN":qN0,
   "absolute_dq_dN_residual":abs(qpb-qN0),
   "zero_qprime_residual_found":bool(abs(qpb-qN0)<1e-8)
 },
 "verdict":{
   "present_state_is_rephased_section18_point":False if best["score"]>=1e-8 else True,
   "interpretation":"Allowing the asymptotic phase u and both Section-18 amplitudes to float is insufficient to reproduce the accepted present Omega_m, q and dq/dN simultaneously." if best["score"]>=1e-8 else "A shifted asymptotic embedding exists.",
   "what_is_ruled_out":"A mere clock rephasing/amplitude relabeling as the explanation of the present-to-mature mismatch." if best["score"]>=1e-8 else "Nothing beyond the tested embedding.",
   "what_remains":"A genuine future transverse/canonicalization dynamics (or an additional derived matching principle) must carry the accepted action from the present overshoot into the mature fixed-point manifold.",
   "parameter_warning":"This audit does not license inserting an arbitrary interpolation or treating a handoff time as a prediction; any new observable coupling/time scale must be derived or included in model-family complexity."
 }
}
(OUT/"minimal_handoff_identifiability.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("MINIMAL_HANDOFF_IDENTIFIABILITY",json.dumps(out,sort_keys=True))
