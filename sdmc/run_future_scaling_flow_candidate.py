#!/usr/bin/env python3
"""
Diagnostic future fixed-point flow for the SDMC active structural exponent.

This is not imposed on hi_class.  It is an analytic normal-form candidate for
the future extension after the accepted z>=0 reconstruction.

Define
    m_X = -d ln rho_X / d ln sigma
and choose the minimal autonomous flow
    d m_X / d ln sigma = gamma m_X (2-m_X).

The exact solution is
    m_X = 2/(1+C sigma^(-2 gamma)),
    C=(2-m0)/m0.

The corresponding density multiplier is
    rho_X/rho_X0 = [(1+C)/(sigma^(2 gamma)+C)]^(1/gamma),

and therefore
    chi/chi0 = sigma^2 rho_X/rho_X0.

Specifying chi_inf fixes gamma:
    gamma = ln(1+C)/ln(chi_inf/chi0).

Three scenarios are reported:
  A. canonical bounded completion: F_inf=1, chi_inf=1;
  B. accepted-F matched completion: F_inf=exp(AF), chi_inf=F_inf,
     which still gives N_inf=Xi_P;
  C. bounded chi with accepted F_inf, giving N_inf=Xi_P/sqrt(F_inf).

The purpose is to quantify how rapidly the present near-m_X=0 source must
flow toward the mature inverse-square m_X=2 fixed point.
"""
from pathlib import Path
import json, math

SRC=Path("output/structural_scaling_exponent_audit.json")
OUT=Path("output/future_scaling_flow_candidate.json")

AF=0.023604633340554453
CHI0=0.9410895692368776
N0=3.4135084646534866
XI=math.sqrt(8.*math.pi/3.)

if SRC.exists():
    d=json.loads(SRC.read_text())
    m0=float(d["present"]["mX_structural"])
    beta0=float(d["present"]["beta_chi"])
else:
    m0=0.06774426095348086
    beta0=1.9322557390465191

C=(2.-m0)/m0
Finf_acc=math.exp(AF)

def scenario(name,chiinf,Finf):
    if chiinf<=CHI0:
        raise RuntimeError("future saturation chi must exceed present chi")
    gamma=math.log1p(C)/math.log(chiinf/CHI0)

    def mx(sigma):
        return 2./(1.+C*sigma**(-2.*gamma))
    def rho_rel(sigma):
        return ((1.+C)/(sigma**(2.*gamma)+C))**(1./gamma)
    def chi(sigma):
        return CHI0*sigma*sigma*rho_rel(sigma)

    thresholds={}
    for mt in [0.1,0.25,0.5,1.0,1.5,1.9,1.99]:
        sig=(C/(2./mt-1.))**(1./(2.*gamma))
        thresholds[str(mt)]={
            "sigma":sig,
            "fractional_sigma_growth":sig-1.,
            "chi":chi(sig),
            "beta_chi":2.-mt,
            "w_if_p1":-1.+mt/3.
        }

    samples=[]
    for sig in [1.,1.005,1.01,1.02,1.03,1.05,1.1,1.2,1.5,2.]:
        m=mx(sig)
        samples.append({
            "sigma":sig,
            "mX":m,
            "beta_chi":2.-m,
            "rhoX_over_rhoX0":rho_rel(sig),
            "chi":chi(sig),
            "w_if_p1":-1.+m/3.
        })

    Ninf=XI*math.sqrt(chiinf/Finf)
    return {
        "name":name,
        "chi_inf":chiinf,
        "F_inf":Finf,
        "gamma":gamma,
        "N_inf":Ninf,
        "N_inf_over_N0":Ninf/N0,
        "thresholds":thresholds,
        "samples":samples,
        "asymptotic_checks":{
            "chi_from_closed_form":CHI0*(1.+C)**(1./gamma),
            "mX_inf":2.0,
            "beta_chi_inf":0.0,
            "w_if_p1_inf":-1./3.
        }
    }

out={
    "status":"analytic future normal-form candidate only; not a fitted or imposed cosmological history",
    "present":{
        "mX0":m0,"beta_chi0":beta0,"chi0":CHI0,"N0":N0,
        "C":C,"Xi":XI
    },
    "scenarios":[
        scenario("canonical_bounded",1.0,1.0),
        scenario("accepted_F_matched",Finf_acc,Finf_acc),
        scenario("bounded_chi_accepted_F",1.0,Finf_acc)
    ]
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("FUTURE_SCALING_FLOW_CANDIDATE")
print("PRESENT",json.dumps(out["present"],sort_keys=True))
for s in out["scenarios"]:
    print("SCENARIO",json.dumps({
        "name":s["name"],"chi_inf":s["chi_inf"],"F_inf":s["F_inf"],
        "gamma":s["gamma"],"N_inf":s["N_inf"],
        "N_inf_over_N0":s["N_inf_over_N0"],
        "sigma_at_m1":s["thresholds"]["1.0"]["sigma"],
        "sigma_at_m1p9":s["thresholds"]["1.9"]["sigma"]
    },sort_keys=True))
