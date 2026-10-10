#!/usr/bin/env python3
"""SDMC nonlinear Jacobian / high-z galaxy / SMBH microphysics audit.

Independent conditional consistency calculations. NOT a full JWST catalogue,
SMBH population likelihood, free-action UV solution, or derived Wilson coefficients.
Preserves the accepted late300 candidate. Python + scipy only.

Manuscript B: spatial map sections 53--56 and nonlinear sections 26--36.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
from scipy.integrate import quad

P = json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]
H0 = P["H0"]
h = H0/100
om = (P["omega_b"] + P["omega_cdm"])/h**2
orad = 4.17998772e-5/h**2
ox = 1-om-orad
c = 299792.458
nt = -math.log1p(P["z_t"])
zhalo = 5.5
a = 1/(1+zhalo)

def E(z):
    n = -math.log1p(z)
    rm,rr = om*(1+z)**3,orad*(1+z)**4
    wb = rr/(3*(rm+rr))
    W = .5*(1-math.tanh((n-nt)/.5))
    ft = W*3*(1+wb)/P["lambda_e"]**2
    d = (P["A_late"]*z*math.exp(-z/.25)
         -P["B_late"]*z*z*math.exp(-z/1.5))
    return math.sqrt(rm+rr+ox)*(1+d)/math.sqrt(1-ft)

hz = H0*E(zhalo)
Hlen = hz/c

# Linearized 3D conformal metric gamma_ij=M^2 a^2 exp(2u) delta_ij:
# Ricci_ij = (k_i k_j+delta_ij k^2) u for Fourier plane wave.
# Ricci_ij Ricci_ij = 6 k^4 u^2; R^2=16 k^4 u^2.
# Dimension=3 identity Riem^2 = 4 Ricci^2 - R^2 = 8 k^4 u^2.
riem_prefactor = 4*6-16
assert riem_prefactor == 8
# Selector B W etaK Riem^2/2, kinetic trace Z_u=9 B f_ell^2:
# alpha_K_spatial = 8 W etaK / (9 f_ell^2 M^4).
# This is symbolic, not a value: W,etaK,f_ell,M still needed.

# Quasistatic response to cubic curvature source:
# R(k,z)=m^2/[m^2+c_ell^2 (k/a)^2+alpha_K(k/a)^4],
# with m and H expressed as inverse lengths, alpha_K [Mpc^2].
def response(k, mass_H, c_ell, alpha):
    m=mass_H*Hlen
    q=k/a
    return m*m/(m*m+c_ell*c_ell*q*q+alpha*q**4)

rows=[]
for mass_H in [5,10,300]:
    m=mass_H*Hlen
    for k in [.553,1.191,2.566,3.0,5.528,11.91]:
        q=k/a
        rows.append(dict(m_over_H=mass_H,k_Mpc_inverse=k,
            response_cell_1_alpha_0=response(k,mass_H,1.,0.),
            cell_max_for_50percent=m/q,
            cell_max_for_90percent=m/(3*q),
            alphaK_max_Mpc2_for_90percent_if_cell_0=m*m/(9*q**4)))
cutoffs=[]
for mass_H in [5,10,300]:
    m=mass_H*Hlen
    q=3./a
    for cell in [0,.001,.01,.1,1.]:
        alpha=(m*m-cell*cell*q*q)/q**4
        cutoffs.append(dict(m_over_H=mass_H,c_ell=cell,
                           alphaK_at_F50_k3_Mpc2=alpha,
                           positive_alphaK_allowed=alpha>=0))

g3_reported=.02048146490100771
g3_95=.03481602388850428
g3_68=.04422747975131273
g3_corrected=[]
for mass_H in [10,300]:
    for k in [.553,1.191]:
        R=response(k,mass_H,1.,0.)
        g3_corrected.append(dict(m_over_H=mass_H,k_Mpc_inverse=k,
            response=R,g3_needed_to_keep_original_95_target_if_local_source_linear=g3_95/R))
# rare-tail heuristic, not an independently normalized mass function.
rare=[]
for nu in [4,5,6,7]:
    for A in [1.,1.1555,1.2122596034,1.36,1.457]:
        rare.append(dict(nu=nu,collapse_amplification=A,
                 abundance_ratio_approx=math.exp(.5*nu**2*(1-A**-2))))

# Convert H0 to Gyr^-1 with 1/H0 = 9.778131/h Gyr.
def ageMyr(z):
    return 1000*9.778131/h*quad(lambda x:1/E(math.expm1(x)),
                               math.log1p(z),math.log1p(1e7),
                               epsrel=1e-10,limit=500)[0]
agez={str(z):ageMyr(z) for z in [5.5,7,10,11.71,20,23.26,30,34.82]}
smbh=[]
for zseed in [20,23.26,30,34.82]:
    dt=agez["7"]-agez[str(zseed)]
    for duty in [1.,.5,.2]:
        fac=math.exp(duty*dt/45.)
        smbh.append(dict(seed_redshift=zseed,target_redshift=7.,
            duty_cycle=duty,elapsed_Myr=dt,growth_factor_ideal_Eddington=fac,
            minimum_seed_Msun_for_1e9_at_z7=1e9/fac))

# Structural sector was suppressed by W -> 0 following thermal handoff,
# so the manuscript's prethermal selector cannot remain an effective halo k^4
# operator without a newly derived reactivation W_halo(z,environment).
out={
 "audit_status":"CONDITIONAL_BOUNDS_ONLY_NO_UNIQUE_UV_CLOSURE",
 "scientific_gates":{
   "g3_microphysically_derived":False,
   "alphaK_numeric_microphysically_derived":False,
   "alphaK_symbolic_derived_with_normalization_convention":True,
   "halo_epoch_selector_reactivation_derived":False,
   "exact_JWST_catalog_likelihood":False,
   "SMBH_seed_accretion_population_fit":False,
   "promote_candidate":False},
 "coefficient_identity":{
   "flat_trace_Riemann_squared_Fourier_prefactor":riem_prefactor,
   "assumed_kinetic_Zu":"9 B f_ell^2 (ell=3u)",
   "alphaK_spatial":"8 W etaK/(9 f_ell^2 M^4)",
   "selector_after_thermal_handoff":"W -> 0 per Manuscript B §55.2; hence alphaK -> 0 without reactivation",
   "not_to_confuse_with":"Horndeski EFT alpha_K kinetic alpha parameter"},
 "late300":{
   "H0":H0,"Omega_m0":om,"Omega_r0":orad,
   "H_z5p5_km_s_Mpc":hz,"H_over_c_Mpc_inv":Hlen,
   "a_z5p5":a,"linear_observables_not_recomputed":True},
 "conditional_g3":{
   "late300_diagnostic_g3_equals_AF":g3_reported,
   "JWST_FRESCO_epsilon20_95_equiv_g3":g3_95,
   "JWST_FRESCO_epsilon20_68_equiv_g3":g3_68,
   "note":"the g3 targets come from a prior likelihood-equivalent transfer, not a UV prediction; earlier identifications g3=AF are unproven"},
 "screening":{
   "response_formula":"m^2/[m^2 + c_ell^2(k/a)^2 + alphaK_spatial(k/a)^4]",
   "gradient_bounds":rows,"alphaK_for_assumed_kcut3":cutoffs,
   "g3_reweight_if_unscreened_c_ell1":g3_corrected},
 "rare_tail_diagnostics":rare,
 "cosmic_ages_Myr":agez,
 "SMBH_ideal_Eddington_examples":smbh,
 "caveats":[
  "g3 multiplies a new u*kappa^3 operator absent from original spatial quadratic action and accepted Horndeski IR action",
  "alphaK depends on UV curvature selector etaK, trace kinetic normalization, mapper M and W; none of these are uniquely fixed",
  "W prethermal curvature selector vanishes by the matter era in the written action; galaxy filter needs an independently derived matter-era reactivation",
  "for c_ell=1 and m/H=5,10,300, halo-scale gradients suppress the quasi-static response; very low c_ell is possible but not derived",
  "the EdS top-hat nonlinear force is a phenomenological closure, not a variation of a complete coupled action",
  "SMBH seed requirements assume ideal exponential Eddington accretion and no interruptions, feedback, mergers or radiative-efficiency variation"
 ]
}
outdir=Path("output/early_galaxy_smbh_uv")
outdir.mkdir(parents=True,exist_ok=True)
(outdir/"microphysics_bounds.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
print("SDMC_JWST_SMBH_UV_AUDIT",json.dumps({
 "status":out["audit_status"],"H_z5p5":hz,
 "g3_95_equivalent":g3_95,
 "alphaK_F50_m10_c0":next(x["alphaK_at_F50_k3_Mpc2"] for x in cutoffs if x["m_over_H"]==10 and x["c_ell"]==0),
 "R_k1p191_m10_cell1":response(1.191,10,1.,0.),
 "delta_t_z20_to_z23p26_Myr":agez["20"]-agez["23.26"],
 "promote_candidate":False},sort_keys=True))
assert all(x["response_cell_1_alpha_0"]<.01 for x in rows if x["k_Mpc_inverse"]==1.191)
assert all(x["positive_alphaK_allowed"] for x in cutoffs if x["c_ell"]==0)
assert agez["20"] > agez["23.26"] > agez["30"]
