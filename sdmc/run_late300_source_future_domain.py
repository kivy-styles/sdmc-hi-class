#!/usr/bin/env python3
from pathlib import Path
import json,math
from scipy.optimize import brentq

# Accepted late300 source coordinates used to reconstruct the frozen action.
AF=0.02048146490100771
ZC=3.927876388467848
WIDTH=0.33114133956842123
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.05181542627513409
OX=0.7014079815036496
LAMBDA_E=18.40625
ZT=16.19317622240633
DNT=0.5
A_LATE=0.024822032167576252
TAUA=0.25
B_LATE=0.01264704344701022
TAUB=1.5

F_MATURE=1.02388542769
D_MATURE=2.0

def S_of_N(N):
    Nc=-math.log1p(ZC)
    return 0.5*(1+math.tanh((N-Nc)/(2*WIDTH)))

def F_of_N(N):
    return math.exp(AF*S_of_N(N))

def D_of_N(N):
    return DFLOOR+D0*S_of_N(N)**POWER

def delta_H_z(z):
    return A_LATE*z*math.exp(-z/TAUA)-B_LATE*z*z*math.exp(-z/TAUB)

def g_z(z):
    return 1.0+delta_H_z(z)

# Future z in (-1,0), a=1/(1+z). Find first zero of the late profile factor.
grid=[-i/10000 for i in range(10001)]
brackets=[]
prevz,prev=g_z(grid[0]),None
lastz=grid[0]; lastg=g_z(lastz)
for z in grid[1:]:
    gz=g_z(z)
    if lastg*gz<0:
        brackets.append((lastz,z))
    lastz,lastg=z,gz
roots=[brentq(g_z,a,b) for a,b in brackets]
root=max(roots) if roots else None # closest-to-present crossing in (-1,0)
if root is not None:
    aroot=1/(1+root)
    Nroot=math.log(aroot)
else:
    aroot=Nroot=None

# Analytic source asymptotes.
F_source_inf=math.exp(AF)
D_source_inf=DFLOOR+D0
F_gap=F_MATURE-F_source_inf
D_gap=D_MATURE-D_source_inf

# Tracker handoff W -> 0 as N->+inf, so it cannot repair the late-profile zero.
# Sample the analytic source ahead of today without running or modifying CLASS.
Ns=[0.0,0.25,0.5,1.0,1.5,2.0,2.5]
samples=[]
for N in Ns:
    a=math.exp(N); z=1/a-1
    samples.append({
      "N":N,"a":a,"z":z,
      "S":S_of_N(N),"F_source":F_of_N(N),"D_source":D_of_N(N),
      "delta_H":delta_H_z(z),"one_plus_delta_H":g_z(z)
    })

out={
  "status":"analytic future-domain audit of the late300 parameterized source construction",
  "important_scope":{
    "accepted_frozen_action_future_prediction":False,
    "what_is_tested":"the analytic sdmc_full / independent-kinetic source functions used upstream to reconstruct late300",
    "why":"to determine whether the source construction itself can supply the manuscript mature future without a new canonicalization regime"
  },
  "source_parameters":{
    "AF":AF,"zc":ZC,"width":WIDTH,"D0":D0,"power":POWER,"Dfloor":DFLOOR,
    "Omega_X0":OX,"lambda_e":LAMBDA_E,"z_t":ZT,"DeltaN":DNT,
    "A_late":A_LATE,"tau_A":TAUA,"B_late":B_LATE,"tau_B":TAUB
  },
  "source_analytic_limits":{
    "F_source_N_to_inf":F_source_inf,
    "D_source_N_to_inf":D_source_inf,
    "alphaM_source_N_to_inf":0.0,
    "alphaB_source_N_to_inf":0.0
  },
  "manuscript_mature_targets":{
    "F_inf":F_MATURE,"D_inf":D_MATURE,
    "F_gap_target_minus_source_limit":F_gap,
    "F_fractional_gap":F_gap/F_MATURE,
    "D_gap_target_minus_source_limit":D_gap,
    "D_fractional_gap":D_gap/D_MATURE
  },
  "late_background_profile":{
    "formula":"1+delta_H(z), delta_H=A z exp(-z/tau_A)-B z^2 exp(-z/tau_B)",
    "future_zero_z":root,
    "future_zero_a":aroot,
    "future_zero_N":Nroot,
    "g_at_z_minus_1_limit":g_z(-1.0),
    "interpretation":"The analytic late-profile multiplier reaches zero at finite future scale factor; beyond that point the sdmc_full source fails its own 1+delta_H>0 admissibility condition."
  },
  "samples":samples,
  "verdict":{
    "source_continues_to_manuscript_mature_branch":False,
    "reason":[
      "The source F sigmoid saturates below the manuscript F_inf.",
      "The source independent-kinetic D saturates near 0.394 rather than D=2.",
      "The sdmc_full late background factor reaches 1+delta_H=0 at finite future a."
    ],
    "scientific_interpretation":"late300's upstream source parameterization is a finite-epoch reconstruction/fit, not the microscopic mature SDMC completion.",
    "accepted_action_implication":"The successful frozen-action replay through today remains valid. This audit only shows that the mature canonical future cannot be obtained by naively continuing the upstream late300 source formulas.",
    "required_next_step":"A distinct, derived canonicalization/transverse-flow regime is required to connect late300 to the manuscript mature action."
  }
}
Path("output").mkdir(exist_ok=True)
Path("output/late300_source_future_domain.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("LATE300_SOURCE_FUTURE_DOMAIN",json.dumps(out,sort_keys=True))
