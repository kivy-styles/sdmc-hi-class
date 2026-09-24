#!/usr/bin/env python3
"""
Matter-loaded mature-coasting bridge for the SDMC structural lapse.

Once the mature structural source has reached
    rho_X ~ a^-2,  w_X = -1/3,  p = d ln sigma/d ln a = 1,
a flat constant-F fluid-level background obeys

    H^2 = H_X^2 a^-2 [1 + mu_m/a + mu_r/a^2],

where mu_m and mu_r are the matter/scalar and radiation/scalar amplitudes
referenced at a=1.  Because p=1 implies sigma proportional to a,

    N/N_inf = sqrt(1 + mu_m/a + mu_r/a^2),

and

    q = 1/2 (mu_m/a + 2 mu_r/a^2)
        / (1 + mu_m/a + mu_r/a^2)  > 0.

Thus residual matter/radiation load the structural speed above its pure-scalar
endpoint and make exact mature coasting approach q=0 from the decelerating side.

The present universe is not yet on this mature manifold (w_X ~ -0.98 rather
than -1/3), so applying the relation at a=1 is explicitly a diagnostic
extrapolation, not an exact fit.  It nevertheless provides a useful numerical
reference inside the independently derived N_inf=2.860--2.894 bracket.
"""
from pathlib import Path
import json, math

OUT=Path("output/matter_loaded_coasting_bridge.json")

H0=69.45160505326464
h=H0/100.
omega_b=0.02208511574370519
omega_cdm=0.12298509428428875
Neff=3.046
Tcmb=2.7255

OMEGA_X0=0.7237401509603739
CHI0=0.9410895692368776
P0=1.034215658562453
N0=3.4135084646534866
AF=0.023604633340554453
FINF=math.exp(AF)
XI=math.sqrt(8.*math.pi/3.)

# Photon density at Tcmb=2.7255 K in the standard CLASS convention.
omega_gamma=2.472e-5*(Tcmb/2.7255)**4
omega_r=omega_gamma*(1.+0.22710731766*Neff)

Omega_m=(omega_b+omega_cdm)/(h*h)
Omega_r=omega_r/(h*h)

mu_m=Omega_m/OMEGA_X0
mu_r=Omega_r/OMEGA_X0

Ninf_loading=N0/math.sqrt(1.+mu_m+mu_r)
chiinf_loading=FINF*(Ninf_loading/XI)**2

Ninf_bounded=XI/math.sqrt(FINF)
Ninf_matched=XI

def bridge(a,Ninf=Ninf_loading):
    den=1.+mu_m/a+mu_r/(a*a)
    return {
      "a":a,
      "N_over_Ninf":math.sqrt(den),
      "N":Ninf*math.sqrt(den),
      "q":.5*(mu_m/a+2.*mu_r/(a*a))/den,
      "dlnN_dlna":-.5*(mu_m/a+2.*mu_r/(a*a))/den,
    }

out={
  "status":(
    "fluid-level mature-manifold diagnostic only; the accepted present action "
    "is still in the activation-driven w_X~-1 handoff and is not assumed to "
    "satisfy the mature p=1,w=-1/3 conditions exactly at a=1"
  ),
  "inputs":{
    "H0":H0,"h":h,"omega_b":omega_b,"omega_cdm":omega_cdm,
    "Omega_m0":Omega_m,"Omega_r0":Omega_r,"Omega_X0":OMEGA_X0,
    "chi0":CHI0,"p0":P0,"N0":N0,"F_inf":FINF,"Xi_v":XI,
  },
  "mature_bridge":{
    "mu_m":mu_m,"mu_r":mu_r,
    "formula_N":"N/Ninf=sqrt(1+mu_m/a+mu_r/a^2)",
    "formula_q":"q=0.5*(mu_m/a+2*mu_r/a^2)/(1+mu_m/a+mu_r/a^2)",
    "Ninf_inferred_if_present_were_on_mature_manifold":Ninf_loading,
    "chiinf_corresponding_at_fixed_Finf":chiinf_loading,
  },
  "independent_structural_bracket":{
    "bounded_chi_Ninf":Ninf_bounded,
    "matched_chi_F_Ninf":Ninf_matched,
    "loading_reference_fraction_from_lower":
      (Ninf_loading-Ninf_bounded)/(Ninf_matched-Ninf_bounded),
  },
  "present_loading_comparison":{
    "lapse_excess_squared_bounded":(N0/Ninf_bounded)**2-1.,
    "lapse_excess_squared_matched":(N0/Ninf_matched)**2-1.,
    "physical_matter_to_active_ratio":mu_m,
    "physical_matter_plus_radiation_to_active_ratio":mu_m+mu_r,
  },
  "samples":[bridge(a) for a in [1.,2.,3.,5.,10.,100.,1000.]]
}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

print("MATTER_LOADED_COASTING_BRIDGE")
print("Omega_m0",Omega_m)
print("Omega_r0",Omega_r)
print("mu_m",mu_m,"mu_r",mu_r)
print("Ninf_loading_reference",Ninf_loading)
print("chiinf_loading_reference",chiinf_loading)
for row in out["samples"]:
    print("BRIDGE",json.dumps(row,sort_keys=True))
