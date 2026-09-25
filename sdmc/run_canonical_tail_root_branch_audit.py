#!/usr/bin/env python3
"""
Compare the two positive C3 No-Slip Planck/braiding tail roots for the
legacy-canonical endpoint.

The quartic tail-root audit finds
    mu_F ~ 2.8610674
and
    mu_F ~ 5.5732190.
Both satisfy the same asymptotic No-Slip ratio and the same accepted present
C3 jet.  This script builds the full canonical action-level continuation for
each root, applies the same on-trajectory No-Slip refinement, and compares
background/health diagnostics through ln(a)=10.

This tests interpolation-branch uniqueness without changing the mature
canonical endpoint.
"""
from pathlib import Path
import json,math
import numpy as np
import run_future_homogeneous_covariant_stitch_audit as aud

OUT=Path("output/canonical_tail_root_branch_audit.json")

ROOTS=[2.8610674139117047,5.573218997214239]
MU_RATES={"k1":50.0,"k2":5.0,"V":120.0}

def build_for_root(muF):
    chiinf=aud.FINF
    Ninf=aud.XI
    Zinf=aud.Z0*(Ninf/aud.N0_STRUCT)**2
    r=0.0
    kappa1=2.*aud.FINF
    kappa2=0.0
    U0=4.*aud.FINF*Zinf
    mug=muF+1.

    AFc=aud.coeffs_match(aud.bars["F"],aud.FINF,muF)
    Agc=aud.coeffs_match(aud.g_xder,0.,mug)
    asym={"k1":kappa1,"k2":kappa2,"V":U0}
    Aq={q:aud.coeffs_match(aud.bars[q],asym[q],MU_RATES[q]) for q in asym}

    def action(sig):
        x=max(0.,math.log(sig))
        out={}
        for q in ["k1","k2","V"]:
            qb=[aud.tail_der(Aq[q],asym[q],MU_RATES[q],x,j) for j in range(3)]
            qx=[]
            for k in range(3):
                qx.append(math.exp(-2.*x)*sum(
                    math.comb(k,j)*((-2.)**(k-j))*qb[j]
                    for j in range(k+1)))
            out[q]=qx[0]
            out[q+"s"]=qx[1]/sig
            out[q+"ss"]=(qx[2]-qx[1])/(sig*sig)

        fx=[aud.tail_der(AFc,aud.FINF,muF,x,j) for j in range(3)]
        out["F"]=fx[0]; out["Fs"]=fx[1]/sig
        out["Fss"]=(fx[2]-fx[1])/(sig*sig)

        gx=[aud.tail_der(Agc,0.,mug,x,j) for j in range(3)]
        out["g"]=gx[0]; out["gs"]=gx[1]/sig
        out["gss"]=(gx[2]-gx[1])/(sig*sig)
        return out

    return {
      "name":f"canonical_muF_{muF:.9f}",
      "chi_inf":chiinf,"r":r,"muQ":50.0,
      "mu_k1":MU_RATES["k1"],"mu_k2":MU_RATES["k2"],"mu_V":MU_RATES["V"],
      "N_inf":Ninf,"Z_inf":Zinf,"muF":muF,"mug":mug,
      "mu_rates":dict(MU_RATES),
      "kappa1":kappa1,"kappa2":kappa2,"U0":U0,
      "action":action,"_base_action":action,
    }

def refined(model):
    m=model
    for _ in range(2):
        _,yy,_=aud.release_trajectory(m)
        m=aud.refine_noslip_model(m,yy,blend_rate=500.)
    return m

rows=[]
for muF in ROOTS:
    try:
        m=refined(build_for_root(muF))
        res=aud.evolve(m,Nmax=10.,npts=2001)
        h=res["health_diagnostics"]; f=res["final"]; p=res["present_eom_check"]
        rows.append({
          "mu_F":muF,"mu_g":muF+1.,"status":"ok",
          "present":{
            "H_rel_vs_free":p["H_rel_vs_free"],
            "p_minus_free":p["p_minus_free"],
            "q_minus_free":p["q_minus_free"],
            "cs2":p["cs2_horndeski"],
            "D":p["D_horndeski"],
          },
          "health":{
            "min_D":h["min_D_horndeski"],
            "min_cs2":h["min_cs2_horndeski"],
            "max_cs2":h["max_cs2_horndeski"],
            "max_abs_alphaB_plus_2alphaM":h["max_abs_alphaB_plus_2alphaM"],
            "Gamma_final":h["Gamma_action_final"],
            "mapper_final":h["mapper_ratio_final"],
          },
          "future":res["future_kinematics"],
          "final":{
            "N_struct":f["N_struct"],"p":f["p"],"q":f["q"],
            "D":f["D_horndeski"],"cs2":f["cs2_horndeski"],
            "F":f["F"],"Gamma":f["Gamma_action"],
            "alphaM":f["alphaM"],"alphaB":f["alphaB"],
          }
        })
    except Exception as e:
        rows.append({"mu_F":muF,"mu_g":muF+1.,"status":"failed","error":repr(e)})

ok=[r for r in rows if r["status"]=="ok"]
out={
  "status":"action-level comparison of all positive C3 canonical F/g tail roots",
  "branches":rows,
  "summary":{
    "positive_root_count":len(ROOTS),
    "action_level_healthy_count":len(ok),
    "endpoint_same_within_current_tolerances":
      bool(len(ok)==2 and
           abs(ok[0]["final"]["N_struct"]-ok[1]["final"]["N_struct"])<5e-4 and
           abs(ok[0]["final"]["p"]-ok[1]["final"]["p"])<5e-4),
    "interpretation":(
      "If both roots are healthy, the mature endpoint is unique under the "
      "closure conditions but the finite Planck/braiding interpolation is not. "
      "Native perturbation tests of the slower root would then be the next "
      "selection test."
    )
  }
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("CANONICAL_TAIL_ROOT_BRANCH_AUDIT")
for row in rows:
    print("BRANCH",json.dumps(row,sort_keys=True))
print("SUMMARY",json.dumps(out["summary"],sort_keys=True))
