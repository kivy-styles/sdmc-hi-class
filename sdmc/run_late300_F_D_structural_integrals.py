#!/usr/bin/env python3
"""Accepted Late300 F/D history integrals and conditional structural mapping.
Scientific scope:
  - Exact: F, D, alphaM from published accepted Late300 target profiles.
  - Conditional: S proportional to ordinary cosmic age t, constant historic
    structural lapse N_struct. NOT a free prediction of the structural field.
  - Do NOT identify the manuscript's N^-3 local rate-source projector with
    the cosmological effective Planck mass F^-1 without a covariant derivation.
"""
import json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

P=json.loads((Path(__file__).parent/"late300_candidate.json").read_text())["parameters"]
AF,ZC,W,DC,DF=P["A_F"],P["z_c"],P["width"],P["D0"],P["D_floor"]
H0=P["H0"]
h=H0/100.
OM=(P["omega_b"]+P["omega_cdm"])/h**2
OR=(4.17998772e-5)/h**2
OX=1.-OM-OR

def U(z):
    n=-np.log1p(z)
    nc=-math.log1p(ZC)
    return 0.5*(1.+np.tanh((n-nc)/(2.*W)))

def alpha_M(z):
    u=U(z)
    return AF/W*u*(1-u)

def F(z):
    return np.exp(AF*U(z))

def D(z):
    return DF+DC*U(z)**P["kinetic_power"]

def H_over_H0(a):
    z=1/a-1
    xr=OR/(OR+OM*a)
    n=math.log(a); nt=-math.log1p(P["z_t"])
    tracker=0.5*(1-math.tanh((n-nt)/0.5))*(3+xr)/P["lambda_e"]**2
    delta=P["A_late"]*z*math.exp(-z/.25)-P["B_late"]*z*z*math.exp(-z/1.5)
    return math.sqrt((OM/a**3+OR/a**4+OX)*(1+delta)**2/(1-tracker))

def H0_age(a):
    return quad(lambda aa:1/(aa*H_over_H0(aa)),0,a,epsabs=1e-11,epsrel=2e-10,limit=220)[0]

def reconstruct(z,z_initial=1e7):
    a=1/(1+z); ni=-math.log1p(z_initial); n=math.log(a)
    alpha_int=quad(lambda nn:alpha_M(math.expm1(-nn)),ni,n,
                   epsabs=1e-13,epsrel=1e-12)[0]
    return dict(lnF_gain_integrated=alpha_int,
         lnF_gain_closed=float(math.log(F(z)/F(z_initial))),
         D_gain_integrated=float((DC/AF)*alpha_int),
         D_gain_closed=float(D(z)-D(z_initial)))

def table():
    t0=H0_age(1)
    rows=[]
    for z in [1090.,30.,20.,P["z_t"],10.,5.,P["z_c"],2.2,1.,0.5,0.]:
        a=1/(1+z)
        t=H0_age(a)
        p=1/(H_over_H0(a)*t)
        mapper=t/(t0*a)  # only under S proportional to t
        m=float(alpha_M(z))
        rows.append({"redshift":float(z),"F":float(F(z)),"D":float(D(z)),
            "alpha_M":m,"Geff_over_today_NoSlip":float(F(0)/F(z)),
            "p_structural_chronological_approx":p,
            "M_mapper_chronological_approx":mapper,
            "beta_mapper_required":m/(p-1),
            "beta_structural_required":m/p})
    return rows

def audit():
    rows=table()
    integrals=[dict(redshift=z,**reconstruct(z)) for z in [30.,P["z_t"],P["z_c"],2.2,0.]]
    maxf=max(abs(x["lnF_gain_integrated"]-x["lnF_gain_closed"]) for x in integrals)
    maxd=max(abs(x["D_gain_integrated"]-x["D_gain_closed"]) for x in integrals)
    assert maxf<1e-12 and maxd<1e-12
    assert abs((DC/AF)-16.713608919761782)<1e-10
    peak=next(x for x in rows if x["redshift"]==P["z_c"])
    present=next(x for x in rows if x["redshift"]==0.)
    early=next(x for x in rows if x["redshift"]==30.)
    assert abs(peak["alpha_M"]-.0154627816)<1e-9
    assert .029<peak["beta_mapper_required"]<.032
    assert .00040<early["beta_mapper_required"]<.00055
    assert present["p_structural_chronological_approx"]>1
    # This conditional cubic projection requires an ~0.68% lapse change
    # if (and ONLY if) it were identified with cosmological gravity.
    required_n=(F(0)/F(1e7))**(1/3)
    return {"scope":"Exact Late300 finite F/D integrals; hypothetical chronological structural map, not autonomous S evolution.",
        "equations":{"F":"exp[A_F U(n)]","D":"D_floor+D0*U(n)",
         "alphaM":"(A_F/w)*U*(1-U)","integrated_ln_F":"integral alphaM dln(a)",
         "integrated_D":"(D0/A_F)*integral alphaM dln(a)",
         "mapper":"M=S/(S0 a)","p":"d ln(S)/d ln(a)",
         "alphaM_over_p":"d ln(F)/d ln(S), if the structural scalar is explicitly identified",
         "alphaM_over_pminus1":"d ln(F)/d ln(M), when p is not 1"},
        "D0_over_AF":DC/AF,"z_c":ZC,"z_t":P["z_t"],"integral_crosschecks":integrals,
        "max_lnF_integral_error":maxf,"max_D_integral_error":maxd,"row_table":rows,
        "illustrative_N3_projection":{"Nratio_required_from_early_to_today_if_F_propto_N_cubed":float(required_n),
          "fractional_N_shift_percent":float((required_n-1)*100),
          "source_status":"Manuscript dynamical N^-3 local rate/source projection postulate; not an identity with cosmological F."},
        "finding":"Structural mapper and F vary on distinguishable profiles; a single constant beta_M is not predicted and fails local-slope uniformity. Structural N is assumed constant only within the historical chronological approximation. No UV source functions are derived."}

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("output/late300_F_D_structural_integral_audit.json"))
    args=parser.parse_args()
    d=audit();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print("F_D_HISTORY_INTEGRALS_VERIFIED",d["max_lnF_integral_error"],d["max_D_integral_error"])
    print("STRUCTURAL_BETA_NONCONSTANT",[(r["redshift"],r["beta_mapper_required"]) for r in d["row_table"] if r["redshift"] in (30.,P["z_c"],0.)])
    print("NO_NEW_FUNDAMENTAL_ACTION_CLAIMED")
if __name__=="__main__":main()
