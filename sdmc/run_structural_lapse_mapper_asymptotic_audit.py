#!/usr/bin/env python3
"""
Audit the exact structural-lapse / matter-mapper identities and the
matter-forced late-time cancellation implied by the released covariant action.

Definitions:
    N_s = t_P S_0 dot(sigma)
    p   = d ln(sigma)/d ln(a)
    M   = K_m/K_m0 = K_p/K_p0 = sigma/a

Exact identities:
    d ln M / d ln a = p - 1
    N_s = t_P S_0 * p * M * a H

For the mature inverse-square scalar with residual conserved matter:
    A = 3/(2 kappa),  kappa = 1 + 2 r
    N_s/N_inf - 1 = -A m + O(a^-2)
    p - 1           = -A m + O(a^-2)
    q               = -A m + O(a^-2)
    M/M_inf - 1     = +A m + O(a^-2)

Hence the leading O(a^-1) matter correction cancels in N_s M:
    (N_s/N_inf)(M/M_inf) = 1 + O(a^-2).

This script evaluates those relations on the C3 No-Slip-refined candidates
stored by run_future_homogeneous_covariant_stitch_audit.py.
"""
from pathlib import Path
import json, math

SRC=Path("output/future_homogeneous_covariant_stitch_audit.json")
OUT=Path("output/structural_lapse_mapper_asymptotic_audit.json")

d=json.loads(SRC.read_text())
acc=d["accepted_present"]
F_inf=float(acc["F_inf"])
Z0=float(acc["Z0"])
N0=float(acc["N0_struct"])

out={
  "status":(
    "Exact lapse/mapper identities plus first-order matter-forced mature "
    "asymptotics evaluated on the C3 No-Slip-refined future candidates."
  ),
  "exact_relations":{
    "lapse":"N_struct=t_P*S0*dot(sigma)",
    "mapper":"M=Km/Km0=Kp/Kp0=sigma/a",
    "mapper_flow":"dlnM/dlna=p-1",
    "kinematic_identity":"N_struct=t_P*S0*p*M*a*H",
    "lapse_flow":"dlnN_struct/dlna=dot(v)/(H*v)",
    "p_flow":"dlnp/dlna=dlnN_struct/dlna+1+q-p",
  },
  "asymptotic_relations":{
    "A":"3/(2*kappa)",
    "kappa":"1+2*r=1/cs_inf^2",
    "lapse":"N/Ninf-1=-A*m+O(a^-2)",
    "bridge":"p-1=-A*m+O(a^-2)",
    "deceleration":"q=-A*m+O(a^-2)",
    "mapper":"M/Minf-1=+A*m+O(a^-2)",
    "product":"(N/Ninf)*(M/Minf)=1+O(a^-2)"
  },
  "candidates":{}
}

for name,row in d["noslip_refined_candidates"].items():
    r=float(row["tail_parameters"]["r"])
    kappa=1.+2.*r
    A=3./(2.*kappa)
    Ninf=float(row["asymptotic_target"]["N_inf"])
    Zinf=Z0*(Ninf/N0)**2
    Rstar=6.*F_inf*Zinf

    samples=sorted(row["samples"],key=lambda x:x["ln_a"])
    late=[s for s in samples if s["ln_a"]>=8.]
    if not late:
        raise RuntimeError(f"{name}: no late samples")
    anchor=late[-1]
    m_anchor=(float(anchor["rho_m_action"])*float(anchor["sigma"])**2/Rstar)
    M_anchor=float(anchor["mapper_ratio_to_present"])

    # First-order integrated mapper law:
    # ln(M/Minf)=A*m + O(a^-2).  Use the latest sample to estimate Minf.
    Minf_est=M_anchor*math.exp(-A*m_anchor)

    rows=[]
    for s in late:
        Na=float(s["ln_a"])
        sig=float(s["sigma"])
        m=float(s["rho_m_action"])*sig*sig/Rstar
        am=A*m
        Nrel=float(s["N_struct"])/Ninf
        Mrel=float(s["mapper_ratio_to_present"])/Minf_est
        predN=math.exp(-am)
        predM=math.exp(+am)
        rows.append({
          "ln_a":Na,
          "a":math.exp(Na),
          "matter_loading_m":m,
          "A_m":am,
          "N_over_Ninf":Nrel,
          "N_minus_first_order_prediction":Nrel-predN,
          "p":float(s["p"]),
          "p_minus_first_order_prediction":float(s["p"])-(1.-am),
          "q":float(s["q"]),
          "q_minus_first_order_prediction":float(s["q"])-(-am),
          "mapper_ratio_to_present":float(s["mapper_ratio_to_present"]),
          "M_over_Minf_est":Mrel,
          "M_minus_first_order_prediction":Mrel-predM,
          "lapse_mapper_product_minus_1":Nrel*Mrel-1.,
        })

    out["candidates"][name]={
      "r":r,
      "kappa":kappa,
      "A":A,
      "N_inf":Ninf,
      "Z_inf":Zinf,
      "Rstar":Rstar,
      "M_inf_est_from_latest_sample":Minf_est,
      "late_rows":rows,
      "latest":{
        "ln_a":rows[-1]["ln_a"],
        "A_m":rows[-1]["A_m"],
        "N_over_Ninf_minus_1":rows[-1]["N_over_Ninf"]-1.,
        "M_over_Minf_minus_1":rows[-1]["M_over_Minf_est"]-1.,
        "lapse_mapper_product_minus_1":
          rows[-1]["lapse_mapper_product_minus_1"],
        "p_minus_1":rows[-1]["p"]-1.,
        "q":rows[-1]["q"],
      }
    }

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("STRUCTURAL_LAPSE_MAPPER_ASYMPTOTIC_AUDIT")
for name,row in out["candidates"].items():
    print(name,json.dumps(row["latest"],sort_keys=True))
