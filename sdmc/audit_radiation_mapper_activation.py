#!/usr/bin/env python3
"""SDMC radiation-epoch mapper/lapse/activation identity and transfer identifiability audit.

Works only with exact manuscript identities and frozen late300 parameters.
This does not modify a Boltzmann solver or establish a new scalar-radiation coupling.
Python standard library only.
"""
import json
import math
from pathlib import Path

p = json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]
Xi = math.sqrt(8*math.pi/3)
Nrad = 3.41350846
p_rad = 2.
omega_rad_endpoint = 0.01
Mstar = 0.02946748423
zstar = 1090.
zeq = 3466.
lambda_e = p["lambda_e"]

def chi_density_route(omega_x, N=Nrad, structural_index=p_rad):
    return omega_x*(N/(structural_index*Xi))**2

def fmass(z):
    n = -math.log1p(z)
    nc = -math.log1p(p["z_c"])
    x = (n-nc)/p["width"]
    u = (math.exp(x)/(1+math.exp(x)) if x < 0 else 1./(1.+math.exp(-x)))
    F = math.exp(p["A_F"]*u)
    alphaM = p["A_F"]/p["width"] * u*(1-u)
    return {"F":F,"F_minus_1":F-1,"alphaM":alphaM,"gate_U":u}

chi_rad = chi_density_route(omega_rad_endpoint)
# Radiation tracker p=2 implies M proportional a. Matter tracker p=3/2 implies M proportional sqrt a.
Meq = Mstar*math.sqrt((1+zstar)/(1+zeq))
def M_rad(z):
    assert z>=zeq
    return Meq*(1+zeq)/(1+z)

# Background action split: rhoX=chi rho_v, rho_v proportional sigma^-2,
# so chi'/chi=OmegaX'/OmegaX + 2 N'/N - 2 p'/p for Xi constant.
# Interacting conservation: 2p=Achi+3(1+wX)+q_J.
# q_J is not thereby *predicted*, since OmegaX',N',p' require dynamics.
w_rad=1/3
Achi=0.
q_J = 2*p_rad-3*(1+w_rad)-Achi

# Independent F(z) gate in frozen late300: check whether it can serve as a sizable
# radiation-era late Planck-mass activation at recombination/equality.
f_gate={str(z):fmass(z) for z in [zstar,zeq,1e6]}

# Heavy relative-mode envelope: u proportional a^-3/2 only after the UV lock forms.
late_damping={}
for zi in [zstar,zeq]:
    late_damping[str(zi)] = (21/(1+zi))**1.5

# Superhorizon curvature conservation under adiabatic δP_nad=0 is conditional,
# but does not generate independent CDM enhancement:
# dot(zeta) = -H δPnad/(rho+P) + gradients = 0 in this limit.
out = {
 "status":"DERIVED_RELATIONS_BUT_NO_RADIATION_CDM_TRANSFER_PREDICTION",
 "definitions":{"M":"sigma/a = R/(a R0)", "p":"1+dlnM/dlna",
                "chi":"rho_X/rho_v", "Xi_v":Xi,
                "N":"structural lapse (not ln a)"},
 "radiation_endpoint":{"p":p_rad,"N":Nrad,"Omega_X":omega_rad_endpoint,
                       "chi_density_route":chi_rad, "dlnM_dln_a":p_rad-1,
                       "dlnchi_dln_a":Achi,"q_J":q_J},
 "nonidentifiability":{"same_Np_different_OmegaX":[
        {"Omega_X":ox,"chi":chi_density_route(ox)} for ox in [0.005,0.01,0.02]],
    "Achi_identity":"dlnOmegaX/dlna + 2 dlnN/dlna - 2 dlnp/dlna",
    "q_J_identity":"2p-3(1+wX)-Achi",
    "radiation_to_CDM_transfer_coefficient":"not specified by background M,p,N",
    "adiabatic_superhorizon_delta_zeta":"0 only if deltaP_nonadiabatic=0 and gradient negligible"},
 "geometry":{"M_star":Mstar,"R0_over_chiJ_star":1/Mstar,
    "M_equality_tracker_approx":Meq,
    "M_z1e6_radiation_tracker_approx":M_rad(1e6),
    "M_star_squared_inverse":1/Mstar**2,
    "same_epoch_R_phys_over_chiJ_star":1/(1+zstar)},
 "late300_F_gate":f_gate,
 "heavy_mode_envelope_u_z20_over_ui":late_damping,
 "conclusion":"M determines p and geometric ratio, not chi or new matter perturbation transfer. Baseline late300 unchanged."
}
assert math.isclose(chi_rad,0.003477149557758102,rel_tol=1e-12)
assert math.isclose(q_J,0.,abs_tol=1e-12)
assert 0<Meq<Mstar<1
assert f_gate[str(zstar)]["alphaM"]<1e-7
assert chi_density_route(0.02)==2*chi_rad
d=Path("output/radiation_mapper_activation")
d.mkdir(parents=True,exist_ok=True)
(d/"diagnostics.json").write_text(json.dumps(out,indent=2)+"\n")
print("SDMC_RADIATION_MAPPER_AUDIT",json.dumps({
  "status":out["status"],"chi_rad":chi_rad,
  "M_star":Mstar,"q_J_if_exact_tracker":q_J,
  "F_minus_1_recombination":f_gate[str(zstar)]["F_minus_1"],
  "alphaM_recombination":f_gate[str(zstar)]["alphaM"],
  "promote_new_candidate":False
},sort_keys=True))
