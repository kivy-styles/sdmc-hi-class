#!/usr/bin/env python3
"""Corrigendum: classify old N=57.89/50.13 vs chronological N~3.41;
reconstruct same late300 age, bridge and Jacobian without cross-branch mixing.

This is a chronology/normalization identity audit, not independent scalar action
evolution or an early-galaxy likelihood. The accepted late300 candidate is read-only.
"""
import json, math
from pathlib import Path
from scipy.integrate import quad

V=json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]
H0=V["H0"]; h=H0/100
Om=(V["omega_b"]+V["omega_cdm"])/h**2
Or=4.17998772e-5/h**2
Ox=1-Om-Or
Y=977.7922216807891 # 1 (km/s/Mpc)^-1 in Gyr
C_GYR_MPC=306.60139378555056
Xi=math.sqrt(8*math.pi/3)
n_legacy_r=Xi*20
n_legacy_m=math.sqrt(2*math.pi)*20
n_accepted_reference=3.41350846
age_reference_Gyr=13.6129780
Mstar_reference=.02946748423
zstar=1090

def E(z):
    x=-math.log1p(z)
    rm=Om*(1+z)**3; rr=Or*(1+z)**4
    b=rm+rr+Ox
    wb=rr/(3*(rm+rr))
    W=.5*(1-math.tanh((x+math.log1p(V["z_t"]))/.5))
    ft=W*3*(1+wb)/(V["lambda_e"]**2)
    d=V["A_late"]*z*math.exp(-z/.25)-V["B_late"]*z*z*math.exp(-z/1.5)
    return math.sqrt(b)*(1+d)/math.sqrt(1-ft)

def cosmic_age(z):
    y=math.log1p(z)
    val,err=quad(lambda x:1/E(math.expm1(x)),y,40,
                 epsabs=8e-13,epsrel=2e-11,limit=400)
    return Y/H0*val

age0=cosmic_age(0.)
age_star=cosmic_age(zstar)
R0= n_accepted_reference*age_reference_Gyr*C_GYR_MPC
Nchrono=R0/(age0*C_GYR_MPC)
p0=Y/(H0*age0)
pstar=Y/(H0*E(zstar)*age_star)
Rstar=Nchrono*age_star*C_GYR_MPC
Mstar=Rstar/(R0/(1+zstar))
reference_Rstar=Mstar_reference*R0/(1+zstar)
pbridge_ref=1.6814565
# Previous incompatible reconstruction mixed an earlier present p0 with later
# late300 H0 (and set R0=N0/p0*c/H0); do not interpret its negative radius
# as an independent failure of the chronological closure.
R0_prev=(n_accepted_reference/1.03421566)*299792.458/H0
mixed_R0_diff=R0-R0_prev
snap={}
for z in (0,7,20,1090,3466,1000000):
    t=cosmic_age(z)
    snap[str(z)]={"age_Gyr":t,"N_const_chrono":Nchrono,
                 "p=1/(H t)":Y/(H0*E(z)*t),
                 "M=(t/t0)*(1+z)":(t/age0)*(1+z)}
out={"status":"MANUSCRIPT_CHRONOLOGICAL_LAPSE_RECONSTRUCTION_CORRECTS_CROSS_BRANCH_NEGATIVE_RADIUS",
 "sources":{"old_N":"Manuscript A Part III Section 17 old direct-active-density canonical scalar tracker",
            "new_N":"Manuscript A Part IX Section 138, Structural Clock Lagrangian Sections 4-6, Manuscript B Part I Sections 12-14,28"},
 "legacy_direct_source_N_r_lambda20":n_legacy_r,
 "legacy_direct_source_N_m_lambda20":n_legacy_m,
 "legacy_values_status":"conditional old chi_act=1 tracker, not mandatory in accepted Manuscript B",
 "revised_prethermal_sequence":{"Planck":1,"kinetic_crossover":2/3,"transfer_midpoint_illustrative":1.6155,"classical_capture_reference":n_accepted_reference},
 "accepted_reference":{"N_late266_chrono":n_accepted_reference,"t0_Gyr":age_reference_Gyr,"p0_original_late266":1.03421566,"Mstar":Mstar_reference},
 "late300_recalibrated":{"t0_Gyr":age0,"N_const_classical":Nchrono,"p0":p0,"pstar":pstar,
 "tstar_Gyr":age_star,"R0_Mpc":R0,"Rstar_Mpc":Rstar,"Mstar":Mstar,
 "Mstar_relative_discrepancy":Mstar/Mstar_reference-1,
 "pstar_relative_discrepancy":pstar/pbridge_ref-1,
 "older_incompatible_R0_Mpc":R0_prev,
 "mixed_R0_deficit_Mpc":mixed_R0_diff,
 "snapshots":snap},
 "qualification":"Chronological closure is a tested postulate/trajectory and not a uniquely derived full all-epoch kinetic solution; do not mix p0 and H0 from separate branch calibrations.",
 "new_candidate_promoted":False}
assert math.isclose(n_legacy_r,57.88810036466141,rel_tol=1e-12)
assert math.isclose(n_legacy_m,50.13256549262,rel_tol=1e-12)
assert age0>13.5 and age0<13.7
assert Mstar>0 and abs(Mstar/Mstar_reference-1)<.005
assert abs(pstar/pbridge_ref-1)<.005
assert Rstar>0 and mixed_R0_diff>50
o=Path("output/lapse_chronology_corrigendum")
o.mkdir(exist_ok=True,parents=True)
(o/"chronology.json").write_text(json.dumps(out,indent=2)+"\n")
print("SDMC_LAPSE_CHRONOLOGY",json.dumps({
 "legacy_radiation":n_legacy_r,"legacy_matter":n_legacy_m,
 "late300_t0_Gyr":age0,"N_const_classical":Nchrono,
 "p0":p0,"pstar":pstar,"Mstar":Mstar,"Rstar_Mpc":Rstar,
 "relative_Mstar_error":Mstar/Mstar_reference-1,
 "normalization_mismatch_Mpc":mixed_R0_diff,
 "status":out["status"],"promote_candidate":False},sort_keys=True))
