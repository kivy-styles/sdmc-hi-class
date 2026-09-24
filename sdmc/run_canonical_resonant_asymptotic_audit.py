#!/usr/bin/env python3
"""
Audit the resonant second-order asymptotics of the closed canonical SDMC
endpoint.

For the r=0, lambda_l=sqrt(2) mature action the autonomous eigenvalues are
-1 (dust) and -2 (intrinsic scalar).  The quadratic self-coupling of the dust
mode is therefore resonant with the scalar mode and produces a^-2 ln(a)
corrections.

This script evaluates the exact homogeneous structural action for the refined
legacy-canonical future stitch and tests the analytic expansion over the
mature interval.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_resonant_asymptotic_audit.json")
XI=math.sqrt(8.*math.pi/3.)

SPEC={
  "name":"legacy_canonical_resonant",
  "chi_inf":aud.FINF,
  "r":0.0,
  "muQ":50.0,
  "mu_k1":50.0,
  "mu_k2":5.0,
  "mu_V":120.0,
}

model=aud.build_refined_noslip_candidate(SPEC,iterations=2,blend_rate=500.)
NN,yy,diag=aud.release_trajectory(model,Nmax=10.,npts=4001)

N=np.asarray(NN)
u=np.exp(-N)
H=np.asarray([d["H"] for d in diag])
F=np.asarray([d["F"] for d in diag])
p=np.asarray([d["p"] for d in diag])
q=np.asarray([d["q"] for d in diag])
Nstruct=np.asarray([d["N_struct"] for d in diag])
rm=np.asarray([d["rho_m_action"] for d in diag])
rr=np.asarray([d["rho_r_action"] for d in diag])
sigma=np.asarray(yy[0])
mapper=sigma/np.exp(N)

Omega_m=rm/(3.*F*H*H)
Omega_r=rr/(3.*F*H*H)
x=p/math.sqrt(3.)

# Canonical potential fraction y^2=V/(3 F H^2).
V=np.empty_like(N)
for i,sig in enumerate(sigma):
    V[i]=model["action"](float(sig))["V"]
y2=V/(3.*F*H*H)
closure=x*x+y2+Omega_m+Omega_r-1.

# From Omega_m = m1 u + 3 m1^2 u^2 + ...,
# solve the quadratic point-by-point for the asymptotic constant m1.
disc=np.maximum(1.+12.*Omega_m,0.)
m1_local=(-1.+np.sqrt(disc))/(6.*u)

# Use a window well after the coefficient release but before u^2 becomes so
# tiny that extracting the independent scalar eigenmode amplifies roundoff.
fit=(N>=5.0)&(N<=7.5)
m1=float(np.median(m1_local[fit]))

eps=m1*u
logu=np.log(u)

# p = 1 - 3/2 eps
#     + u^2[(9/4)m1^2 log u -45/8 m1^2 -sqrt(6) s2] + ...
forced_p=(1.-1.5*eps
          +u*u*((9./4.)*m1*m1*logu-(45./8.)*m1*m1))
s2_from_p=-(p-forced_p)/(math.sqrt(6.)*u*u)

# q = -3/2 eps
#     + u^2[(9/2)m1^2 log u -9/2 m1^2 -2sqrt(6) s2] + ...
forced_q=(-1.5*eps
          +u*u*((9./2.)*m1*m1*logu-(9./2.)*m1*m1))
s2_from_q=-(q-forced_q)/(2.*math.sqrt(6.)*u*u)

s2=float(np.median(np.r_[s2_from_p[fit],s2_from_q[fit]]))

def predictions(i):
    ui=float(u[i]); li=float(logu[i]); e=m1*ui
    p2=(1.-1.5*e
        +ui*ui*((9./4.)*m1*m1*li-(45./8.)*m1*m1-math.sqrt(6.)*s2))
    q2=(-1.5*e
        +ui*ui*((9./2.)*m1*m1*li-(9./2.)*m1*m1-2.*math.sqrt(6.)*s2))
    nr=(1.-1.5*e
        +ui*ui*((27./8.)*m1*m1*li-(45./8.)*m1*m1
                 -(3.*math.sqrt(6.)/2.)*s2))
    om=e+3.*e*e
    mrel=(1.+1.5*e
          +ui*ui*(-(9./8.)*m1*m1*li+(9./2.)*m1*m1
                   +(math.sqrt(6.)/2.)*s2))
    prod=(1.+ui*ui*((9./4.)*m1*m1*li-(27./8.)*m1*m1
                     -math.sqrt(6.)*s2))
    hsig=(1.+ui*ui*((9./8.)*m1*m1*li
                     -(math.sqrt(6.)/2.)*s2))
    return dict(p=p2,q=q2,Nrel=nr,Omega_m=om,
                mapper_rel=mrel,lapse_mapper_product=prod,
                Hsigma_rel=hsig)

# Estimate M_inf independently from the analytic mapper correction throughout
# the fit window, then use one constant in all comparisons.
Minf_samples=[]
for i in np.where(fit)[0]:
    pr=predictions(i)
    Minf_samples.append(float(mapper[i]/pr["mapper_rel"]))
Minf=float(np.median(Minf_samples))

samples=[]
for nt in [5.,6.,7.,8.,9.,10.]:
    i=int(np.argmin(np.abs(N-nt)))
    pr=predictions(i)
    nrel=float(Nstruct[i]/XI)
    mrel=float(mapper[i]/Minf)
    samples.append({
      "ln_a":float(N[i]),
      "u":float(u[i]),
      "Omega_m":float(Omega_m[i]),
      "Omega_m_pred":pr["Omega_m"],
      "p":float(p[i]),"p_pred":pr["p"],
      "q":float(q[i]),"q_pred":pr["q"],
      "N_over_Ninf":nrel,"N_over_Ninf_pred":pr["Nrel"],
      "M_over_Minf":mrel,"M_over_Minf_pred":pr["mapper_rel"],
      "NM_product":nrel*mrel,
      "NM_product_pred":pr["lapse_mapper_product"],
      "N_over_p_relative":float((Nstruct[i]/p[i])/XI),
      "N_over_p_relative_pred":pr["Hsigma_rel"],
      "Omega_r":float(Omega_r[i]),
      "canonical_closure_residual":float(closure[i]),
    })

# Report scaling of the residual after subtracting the first-order dust term.
# A resonant a^-2 ln a term should make res*a^2 approximately affine in N.
late=(N>=5.)&(N<=10.)
first_p=p-(1.-1.5*m1*u)
first_q=q-(-1.5*m1*u)
first_N=Nstruct/XI-(1.-1.5*m1*u)

out={
  "status":(
    "second-order resonant asymptotic audit of the refined r=0 canonical "
    "structural action"
  ),
  "analytic":{
    "fixed_point":"x=1/sqrt(3), y=sqrt(2/3)",
    "eigenvalues":[-1.0,-2.0],
    "resonance":"2*(-1)=-2 -> a^-2 ln(a) nonlinear correction",
    "Omega_m":"m1*u+3*m1^2*u^2+...",
    "p":"1-3/2*m1*u+u^2[(9/4)m1^2 ln u-45/8 m1^2-sqrt(6)s2]+...",
    "q":"-3/2*m1*u+u^2[(9/2)m1^2 ln u-9/2 m1^2-2sqrt(6)s2]+...",
    "Nrel":"1-3/2*m1*u+u^2[(27/8)m1^2 ln u-45/8 m1^2-(3sqrt(6)/2)s2]+...",
    "mapper":"1+3/2*m1*u+u^2[-9/8*m1^2 ln u+9/2 m1^2+(sqrt(6)/2)s2]+..."
  },
  "fit":{
    "window_ln_a":[5.0,7.5],
    "m1":m1,
    "m1_local_median":float(np.median(m1_local[fit])),
    "m1_local_min":float(np.min(m1_local[fit])),
    "m1_local_max":float(np.max(m1_local[fit])),
    "s2":s2,
    "s2_from_p_median":float(np.median(s2_from_p[fit])),
    "s2_from_q_median":float(np.median(s2_from_q[fit])),
    "M_inf_est":Minf,
    "canonical_closure_max_abs_fit":
      float(np.max(np.abs(closure[fit]))),
  },
  "first_order_residual_scaled":{
    "p_a2_at_5":float(first_p[np.argmin(np.abs(N-5.))]*math.exp(10.)),
    "p_a2_at_10":float(first_p[-1]*math.exp(20.)),
    "q_a2_at_5":float(first_q[np.argmin(np.abs(N-5.))]*math.exp(10.)),
    "q_a2_at_10":float(first_q[-1]*math.exp(20.)),
    "N_a2_at_5":float(first_N[np.argmin(np.abs(N-5.))]*math.exp(10.)),
    "N_a2_at_10":float(first_N[-1]*math.exp(20.)),
    "note":"non-constant a^2-scaled residual is expected because the second-order term contains ln(u)=-ln(a)"
  },
  "samples":samples,
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("CANONICAL_RESONANT_ASYMPTOTIC_AUDIT")
print("FIT",json.dumps(out["fit"],sort_keys=True))
for row in samples:
    print("SAMPLE",json.dumps(row,sort_keys=True))
