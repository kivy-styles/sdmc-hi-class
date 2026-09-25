#!/usr/bin/env python3
"""
Quantify how the two positive C3 Planck/braiding tail roots differ only through
finite transition memory while sharing the same canonical infrared fixed point.

For each root we:
  * integrate the same refined r=0 canonical action to ln(a)=10;
  * extract the resonant matter amplitude m1 from Omega_m;
  * infer the frozen mapper M_inf using the second-order resonant formula;
  * predict m1 independently from present matter normalization and M_inf.

This distinguishes universal IR data (fixed point, eigenvalue) from
branch-dependent memory amplitudes.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_tail_root_memory_audit.json")
ROOTS=[2.8610674139117047,5.573218997214239]
BASE={
  "name":"canonical_root_memory",
  "chi_inf":aud.FINF,
  "r":0.0,
  "muQ":50.0,
  "mu_k1":50.0,
  "mu_k2":5.0,
  "mu_V":120.0,
}
XI=aud.XI

def analyze(muF):
    spec={**BASE,"name":f"canonical_root_{muF:.9f}","muF_override":muF}
    model=aud.build_refined_noslip_candidate(spec,iterations=2,blend_rate=500.)
    NN,yy,diag=aud.release_trajectory(model,Nmax=10.,npts=4001)
    N=np.asarray(NN); u=np.exp(-N)
    sig=np.asarray(yy[0])
    mapper=sig/np.exp(N)
    H=np.asarray([d["H"] for d in diag])
    F=np.asarray([d["F"] for d in diag])
    p=np.asarray([d["p"] for d in diag])
    q=np.asarray([d["q"] for d in diag])
    Ns=np.asarray([d["N_struct"] for d in diag])
    rm=np.asarray([d["rho_m_action"] for d in diag])
    Om=rm/(3.*F*H*H)

    disc=np.maximum(1.+12.*Om,0.)
    m1_local=(-1.+np.sqrt(disc))/(6.*u)
    fit=(N>=5.)&(N<=7.5)
    m1=float(np.median(m1_local[fit]))

    logu=np.log(u)
    forced_p=(1.-1.5*m1*u
              +u*u*((9./4.)*m1*m1*logu-(45./8.)*m1*m1))
    forced_q=(-1.5*m1*u
              +u*u*((9./2.)*m1*m1*logu-(9./2.)*m1*m1))
    s2p=-(p-forced_p)/(math.sqrt(6.)*u*u)
    s2q=-(q-forced_q)/(2.*math.sqrt(6.)*u*u)
    s2=float(np.median(np.r_[s2p[fit],s2q[fit]]))

    mrel=(1.+1.5*m1*u
          +u*u*(-(9./8.)*m1*m1*logu+(9./2.)*m1*m1
                 +(math.sqrt(6.)/2.)*s2))
    Minf=float(np.median(mapper[fit]/mrel[fit]))

    Om0=float(aud.rho_m0/(3.*aud.H0*aud.H0))
    m1_pred=(Om0*Minf*Minf/(aud.FINF*aud.p0_bg*aud.p0_bg)
             *(aud.N0_STRUCT/XI)**2)

    i5=int(np.argmin(abs(N-5.)))
    i10=-1
    return {
      "mu_F":muF,
      "mu_g":muF+1.,
      "M_inf":Minf,
      "m1_fit":m1,
      "m1_from_present_mapper_closure":m1_pred,
      "m1_fractional_closure_error":m1/m1_pred-1.,
      "s2_intrinsic_mode":s2,
      "at_ln_a_5":{
        "N_struct":float(Ns[i5]),"p":float(p[i5]),"q":float(q[i5]),
        "mapper":float(mapper[i5]),"Omega_m":float(Om[i5]),
      },
      "at_ln_a_10":{
        "N_struct":float(Ns[i10]),"p":float(p[i10]),"q":float(q[i10]),
        "mapper":float(mapper[i10]),"Omega_m":float(Om[i10]),
      },
    }

rows=[analyze(r) for r in ROOTS]
slow,fast=rows
out={
  "status":(
    "two-root canonical transition-memory audit: same IR fixed point and "
    "eigenvalue, slightly different frozen mapper and matter-mode amplitude"
  ),
  "branches":rows,
  "differences_fast_minus_slow":{
    "delta_M_inf":fast["M_inf"]-slow["M_inf"],
    "fractional_M_inf":fast["M_inf"]/slow["M_inf"]-1.,
    "delta_m1_fit":fast["m1_fit"]-slow["m1_fit"],
    "fractional_m1_fit":fast["m1_fit"]/slow["m1_fit"]-1.,
    "delta_s2":fast["s2_intrinsic_mode"]-slow["s2_intrinsic_mode"],
    "delta_N_ln_a_10":
      fast["at_ln_a_10"]["N_struct"]-slow["at_ln_a_10"]["N_struct"],
    "delta_p_ln_a_10":
      fast["at_ln_a_10"]["p"]-slow["at_ln_a_10"]["p"],
    "delta_q_ln_a_10":
      fast["at_ln_a_10"]["q"]-slow["at_ln_a_10"]["q"],
  },
  "universal":{
    "N_inf":XI,
    "p_inf":1.0,
    "q_inf":0.0,
    "tangent_eigenvalue":1.0,
    "canonical_shape_r":0.0,
    "cs2_inf":1.0,
  },
  "interpretation":{
    "endpoint":"branch independent",
    "infrared_eigenvalue":"branch independent",
    "transition_memory":(
      "M_inf and therefore the coefficient m1 retain a small dependence on "
      "which healthy C3 Planck/braiding root is followed"
    ),
    "consequence":(
      "the roots belong to the same IR universality class but correspond to "
      "slightly different amplitudes along the universal matter eigenmode"
    )
  }
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("CANONICAL_TAIL_ROOT_MEMORY_AUDIT")
for row in rows:
    print("BRANCH_MEMORY",json.dumps(row,sort_keys=True))
print("DIFFERENCES",json.dumps(out["differences_fast_minus_slow"],sort_keys=True))
