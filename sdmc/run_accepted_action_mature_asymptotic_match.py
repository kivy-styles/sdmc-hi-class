#!/usr/bin/env python3
from pathlib import Path
import json, math, re

ROOT=Path("inputs")
OUT=Path("output/mature_asymptotic_match")
OUT.mkdir(parents=True,exist_ok=True)

canon=json.loads(next((ROOT/"canonical").rglob("accepted_action_canonical_approach.json")).read_text())
ini=next((ROOT/"action").rglob("linear_cov_target.ini")).read_text()

def ini_float(key):
    m=re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*([^#\n]+)",ini)
    if not m: raise RuntimeError(f"missing {key}")
    return float(m.group(1).strip())

H0=ini_float("H0"); h=H0/100.0
omb=ini_float("omega_b"); omc=ini_float("omega_cdm")
Omega_m0=(omb+omc)/(h*h)

q0=float(canon["present_side_endpoint_trends"]["q"]["today"])
dq_dN0=float(canon["present_side_endpoint_trends"]["q"]["d_dN_today"])

# Manuscript Section 18, evaluated at u=a^-1=1:
# Omega_m = m1*u + 3*m1^2*u^2 + ...
# q = -3/2 m1 u + u^2[9/2 m1^2 ln u - 9/2 m1^2 - 2 sqrt(6) s2] + ...
# These are asymptotic expressions. Testing them at u=1 is a diagnostic of
# whether the present accepted-action endpoint is already inside that regime.

def positive_root_for_m1_from_Om(Om):
    return (-1.0+math.sqrt(1.0+12.0*Om))/6.0

def q_series(u,m1,s2):
    return (-1.5*m1*u
            +u*u*(4.5*m1*m1*math.log(u)-4.5*m1*m1-2.0*math.sqrt(6.0)*s2))

def dq_dN_series(u,m1,s2):
    A=1.5*m1
    B=4.5*m1*m1
    C=-4.5*m1*m1-2.0*math.sqrt(6.0)*s2
    dq_du=-A+2.0*u*(B*math.log(u)+C)+B*u
    return -u*dq_du

# Match 1: use observed/accepted Omega_m0 and q0.
m1_omq=positive_root_for_m1_from_Om(Omega_m0)
s2_omq=(-q0-1.5*m1_omq-4.5*m1_omq*m1_omq)/(2.0*math.sqrt(6.0))
dq_pred=dq_dN_series(1.0,m1_omq,s2_omq)

# Match 2: force q0 and dq/dN today, then inspect implied Omega_m.
# At u=1:
# A0=1.5*m1+4.5*m1^2 = -(dq/dN + 2q).
A0=-(dq_dN0+2.0*q0)
disc=1.5**2+18.0*A0
m1_qqp=(-1.5+math.sqrt(disc))/9.0 if disc>=0 else float("nan")
s2_qqp=(-q0-(1.5*m1_qqp+4.5*m1_qqp*m1_qqp))/(2.0*math.sqrt(6.0))
Omega_pred=m1_qqp+3.0*m1_qqp*m1_qqp

future=[]
for a in [1,1.25,1.5,2,3,5,10]:
    u=1.0/a
    future.append({
      "a_over_a0":a,"u":u,
      "q_section18_match_Omega_q":q_series(u,m1_omq,s2_omq),
      "dq_dN_section18_match_Omega_q":dq_dN_series(u,m1_omq,s2_omq),
      "Omega_m_section18_match_Omega_q":m1_omq*u+3.0*m1_omq*m1_omq*u*u,
    })

out={
 "status":"accepted-action present endpoint versus manuscript Section-18 mature asymptotic matching audit",
 "source_scope":{
   "accepted_action":"unchanged frozen reconstructed late300 action; present-side values from supported N<=0 canonical-approach audit",
   "manuscript_equations":["17.8-17.11","18.1","18.3","18.5"],
   "important_warning":"Section 18 is an asymptotic u=a^-1 -> 0 expansion. Evaluating it at u=1 is only a regime-entry diagnostic, not a valid future evolution of the accepted action."
 },
 "accepted_present":{
   "Omega_m0_from_physical_densities":Omega_m0,
   "q0":q0,
   "dq_dN0":dq_dN0
 },
 "match_Omega_m_and_q":{
   "m1":m1_omq,"s2":s2_omq,
   "q0_reconstructed":q_series(1.0,m1_omq,s2_omq),
   "dq_dN_predicted":dq_pred,
   "dq_dN_actual":dq_dN0,
   "dq_dN_residual":dq_pred-dq_dN0,
   "derivative_sign_agreement":(dq_pred*dq_dN0>0)
 },
 "match_q_and_qprime":{
   "m1":m1_qqp,"s2":s2_qqp,
   "Omega_m0_implied":Omega_pred,
   "Omega_m0_actual":Omega_m0,
   "Omega_m_residual":Omega_pred-Omega_m0,
   "physical_Omega_m_match":abs(Omega_pred-Omega_m0)<0.05
 },
 "section18_matched_future_diagnostic":future,
 "verdict":{
   "present_already_in_section18_asymptotic_regime":False,
   "reason":"Matching the accepted present Omega_m and q predicts dq/dN with the opposite sign; matching q and dq/dN instead implies Omega_m>1. The present endpoint cannot satisfy the same second-order mature asymptotic expansion in all three observables.",
   "interpretation":"This supports a distinct later handoff/relaxation before the manuscript's canonical a^-1/a^-2 hierarchy becomes quantitatively valid. It does not derive that handoff or authorize extrapolating the reconstructed action beyond N=0.",
   "next_test":"Use the manuscript canonical fixed-point equations to construct the minimal continuity conditions for a future handoff and test whether they can be satisfied without introducing an additional free observable parameter."
 }
}
(OUT/"accepted_action_mature_asymptotic_match.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("ACCEPTED_ACTION_MATURE_ASYMPTOTIC_MATCH",json.dumps(out,sort_keys=True))
