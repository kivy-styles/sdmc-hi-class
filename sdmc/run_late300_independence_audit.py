#!/usr/bin/env python3
"""Numerical identifiability of the *selected Late300 target functions*, not
an estimate of fundamental SDMC action constants or likelihood Fisher rank.
Code/parameter source: sdmc/late300_candidate.json and accepted SDMC patches.
"""
import argparse
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/"late300_candidate.json").read_text())["parameters"]
THETA=("lambda_e","z_t","A_late","B_late","A_F","z_c","width",
       "D_floor","D0","kinetic_power")
RADIATION_OMEGA_H2=4.17998772e-5

def flat_density(p=P):
    h=p["H0"]/100.0
    om=(p["omega_b"]+p["omega_cdm"])/h**2
    orad=RADIATION_OMEGA_H2/h**2
    return om,orad,1.0-om-orad

def target(z,p=P):
    """Return logarithms of the published E(z), F(z), and D(z) profiles."""
    z=np.asarray(z,float)
    n=-np.log1p(z)  # e-fold time, NOT the SDMC structural lapse
    om,orad,ox=flat_density(p)
    matter=om*(1+z)**3
    radiation=orad*(1+z)**4
    radiative_fraction=radiation/(matter+radiation)
    nt=-np.log1p(p["z_t"])
    wtr=.5*(1.0-np.tanh((n-nt)/.5))
    ftr=wtr*(3+radiative_fraction)/p["lambda_e"]**2
    deformation=(p["A_late"]*z*np.exp(-z/.25)
                 -p["B_late"]*z*z*np.exp(-z/1.5))
    e2=(matter+radiation+ox)*(1+deformation)**2/(1-ftr)
    nc=-np.log1p(p["z_c"])
    u=.5*(1.0+np.tanh((n-nc)/(2*p["width"])))
    f=np.exp(p["A_F"]*u)
    d=p["D_floor"]+p["D0"]*u**p["kinetic_power"]
    if np.any(e2<=0) or np.any(d<=0):
        raise ValueError("Nonpositive physical target profile")
    return {"lnE":.5*np.log(e2),"lnF":np.log(f),"lnD":np.log(d)}

def jacobian(z):
    columns={"lnE":[],"lnF":[],"lnD":[]}
    for name in THETA:
        val=float(P[name])
        step=max(abs(val)*1e-6,1e-8)
        pp=dict(P);pm=dict(P)
        pp[name]=val+step;pm[name]=val-step
        u=target(z,pp);v=target(z,pm)
        for k in columns:
            columns[k].append(val*(u[k]-v[k])/(2*step))
    return {k:np.column_stack(columns[k]) for k in columns}

def svd_rank(array, names, relative_threshold=1e-7):
    norm=np.linalg.norm(array,axis=0)
    keep=norm>1e-14
    sv=np.linalg.svd(array[:,keep]/norm[keep],compute_uv=False)
    return {"rank":int(np.count_nonzero(sv>relative_threshold*sv[0])),
            "active_parameters":[name for name,yes in zip(names,keep) if yes],
            "singular_values":[float(x) for x in sv],
            "smallest_over_largest":float(sv[-1]/sv[0]),
            "meaning":"rank of differentiable target functions; not an observational likelihood rank"}

def audit():
    z=np.expm1(np.linspace(0,np.log1p(1e5),1001))
    j=jacobian(z)
    stacked=np.vstack([j["lnE"],j["lnF"],j["lnD"]])
    joint=svd_rank(stacked,THETA)
    h=svd_rank(j["lnE"][:,:4],THETA[:4])
    f=svd_rank(j["lnF"][:,4:7],THETA[4:7])
    d=svd_rank(j["lnD"][:,7:],THETA[7:])
    conditioning={}
    for top in (0.5,1.,3.,10.,30.,100.,1100.):
        zz=np.expm1(np.linspace(0,np.log1p(top),801))
        jj=jacobian(zz)
        r=svd_rank(np.vstack([jj["lnE"],jj["lnF"],jj["lnD"]]),THETA)
        conditioning[str(top)]={"rank":r["rank"],"smallest_over_largest":r["smallest_over_largest"]}
    om,orr,ox=flat_density(P)
    candidate=2*np.pi/P["lambda_e"]
    canon=16/P["lambda_e"]**2
    out={
       "scope":"phenomenological target function sensitivity, NOT microscopic parameter count",
       "source":"sdmc/late300_candidate.json",
       "ranks":{"H":h,"F":f,"D_conditional":d,"H_F_D":joint},
       "redshift_conditioning":conditioning,
       "omega_flat_closure":{"omega_m":om,"omega_r":orr,"omega_x":ox,
                             "residual":ox-P["Omega_x0"]},
       "D_analytic_checks":{"canonical_radiation_Dfloor":canon,
          "fitted_Dfloor":P["D_floor"],"canonical_floor_fractional_gap":P["D_floor"]/canon-1,
          "candidate_2pi_over_lambda":candidate,"fitted_D0":P["D0"],
          "candidate_D0_fractional_gap":candidate/P["D0"]-1},
       "model_selection":{"six_shared_ordinary_parameters":["H0","omega_b","omega_cdm",
                "A_s","n_s","tau_reio"],
          "sdmc_specific_listed":11,"omega_x_algebraically_closed":1,
          "locally_distinguishable_remaining_target_coefficients":10,
          "additional_retuning_during_frozen_action_replay":0,
          "physical_family_delta_k":"NOT DERIVED",
          "warning":"The full-family action and priors may impose additional relations. This rank is not delta_k."},
       "source_guardrail":"Accepted action is a covariant reconstruction from fitted target H,F,D. Successful free evolution validates that reconstructed action, not an independent upstream prediction of its coefficients."
    }
    assert abs(out["omega_flat_closure"]["residual"])<1e-12
    assert (h["rank"],f["rank"],d["rank"],joint["rank"])==(4,3,3,10)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=Path("output/late300_local021_independence.json"))
    args=ap.parse_args()
    out=audit()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("TARGET_FUNCTION_RANKS",
          {k:v["rank"] for k,v in out["ranks"].items()})
    print("STRUCTURAL_FAMILY_DELTA_K",out["model_selection"]["physical_family_delta_k"])
if __name__=="__main__":
    main()
