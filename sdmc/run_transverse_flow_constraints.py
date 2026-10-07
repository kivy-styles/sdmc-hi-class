#!/usr/bin/env python3
from pathlib import Path
import json, math

ROOT=Path("inputs")
OUT=Path("output/transverse_flow_constraints")
OUT.mkdir(parents=True,exist_ok=True)

def load_one(d,name):
    hits=list((ROOT/d).rglob(name))
    if len(hits)!=1:
        raise RuntimeError(f"expected one {name} under {d}, got {hits}")
    return json.loads(hits[0].read_text())

approach=load_one("approach","accepted_action_canonical_approach.json")
chron=load_one("chronology","accepted_action_handoff_chronology.json")

ep=approach["present_side_endpoint_trends"]
today={k:v["today"] for k,v in ep.items()}
target={k:v["mature_target"] for k,v in ep.items()}
slope={k:v["d_dN_today"] for k,v in ep.items()}

# Local linear scales are diagnostics only: gap divided by present derivative magnitude.
# They are explicitly NOT future forecasts.
local_scale={}
for k in ["F","alphaM","alphaB","D","cs2","gH2","q","wtot"]:
    gap=target[k]-today[k]
    ds=slope[k]
    local_scale[k]={
      "gap_target_minus_today":gap,
      "present_d_dN":ds,
      "present_sign_points_toward_target":bool(gap*ds>0),
      "abs_gap_over_abs_present_slope":None if ds==0 else abs(gap/ds)
    }

# Topological/continuity constraints. q and w are currently below their mature
# targets and still moving downward; any continuous route to mature coasting
# must turn around at least once in the future.
turning_required={
  "q":{
    "today":today["q"],"target":target["q"],"d_dN_today":slope["q"],
    "future_turning_point_required":bool(today["q"]<target["q"] and slope["q"]<0),
    "required_after_turn":"dq/dN must become positive for some later interval"
  },
  "wtot":{
    "today":today["wtot"],"target":target["wtot"],"d_dN_today":slope["wtot"],
    "future_turning_point_required":bool(today["wtot"]<target["wtot"] and slope["wtot"]<0),
    "required_after_turn":"dw_tot/dN must become positive for some later interval"
  }
}

events=chron["events"]
event_summary=[{
 "name":e["name"],"z":e["state"]["z"],"N":e["state"]["N"],
 "F":e["state"]["F"],"D":e["state"]["D"],"cs2":e["state"]["cs2"],
 "alphaM":e["state"]["alphaM"],"gH2":e["state"]["gH2"],
 "q":e["state"]["q"],"wtot":e["state"]["wtot"]
} for e in events]

out={
 "status":"necessary-condition synthesis for the missing accepted-action transverse canonicalization flow",
 "inputs":{
   "accepted_action":"unchanged frozen linear-G3 covariant action; free evolution through today",
   "future_extrapolation_used":False,
   "manuscript_endpoint":{
     "F_inf":1.02388542769,
     "alphaM_inf":0.0,
     "alphaB_inf":0.0,
     "physical_G3_braiding_inf":0.0,
     "D_inf":2.0,
     "cs2_inf":1.0,
     "q_inf":0.0,
     "wtot_inf":-1.0/3.0,
     "G3_inf":0.0,
     "kinetic_shape":"linear f(Z), f_ZZ=0 under mature minimal-derivative/robust-luminality selection"
   }
 },
 "supported_handoff_chronology":event_summary,
 "present_local_directions":local_scale,
 "continuity_constraints":turning_required,
 "necessary_conditions":{
   "geometry_and_braiding":[
     "F must continue to F_inf while alpha_M tends to zero.",
     "Exact No-Slip must remain satisfied and the physical G3 braiding contribution must vanish.",
     "The flow must end on G3=0 and constant G4=F_inf/2 without spoiling the already verified late300 replay."
   ],
   "kinetic_sector":[
     "D must rise from its present value to 2 while remaining strictly positive.",
     "c_s^2 must rise from its present value to 1 without crossing a gradient instability; exact mature robust luminality requires the connected mature kinetic shape to become linear in Z.",
     "The higher-kinetic curvature must disappear in the mature regime rather than merely become observationally small on one trajectory."
   ],
   "background_kinematics":[
     "Because q and w_tot crossed coasting at z~0.707 and are now below the mature values while still moving downward, any continuous completion to q=0 and w_tot=-1/3 requires at least one future kinematic turning point.",
     "The future flow therefore cannot be a monotonic continuation of the present q and w_tot slopes."
   ],
   "dynamical_attractor":[
     "The final canonical branch must preserve the manuscript's attractive matter and scalar modes (eigenvalues -1 and -2).",
     "Finite-handoff freedom may be multi-dimensional, but the mature flow must collapse onto the one-dimensional matter-loaded infrared manifold described in the manuscript."
   ],
   "statistical_and_provenance":[
     "The microscopic law must not be replaced by an arbitrary fitted interpolation and then counted as a prediction.",
     "Any new coupling that affects accepted-trajectory observables must be either independently derived or included in the model-family prior/complexity accounting.",
     "The rejected null-G2 stabilizer is not part of this accepted-action flow."
   ]
 },
 "diagnostic_conclusion":{
   "already_correct_direction":["F","alphaM","alphaB","D","cs2","gH2"],
   "currently_wrong_direction":["q","wtot"],
   "largest_present_fractional_action_deficit":"D",
   "interpretation":"The accepted action already shuts down Planck-mass running and physical braiding, and its kinetic sector has turned toward the canonical targets after the sound-speed trough. The unresolved canonicalization is dominated by the large kinetic-shape deficit plus the required reversal of the present acceleration overshoot.",
   "what_is_not_derived":"No microscopic beta function or future trajectory is inferred."
 }
}
(OUT/"transverse_flow_constraints.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("TRANSVERSE_FLOW_CONSTRAINTS",json.dumps(out,sort_keys=True))
