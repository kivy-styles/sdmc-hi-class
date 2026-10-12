#!/usr/bin/env python3
"""Auditable accepted Late300 baryon loading and H0/acoustic-scale reinterpretation.

THREE DISTINCT TESTS:
(A) accepted physical Jordan baryons/photons, no structural rescaling (valid);
(B) FORCED baryons-only S^-2 instead of a^-3, photons conserved (unphysical
    without a baryon source; acoustic-only stress test, NOT model fit);
(C) FORCED photons and baryons both S^-2 (unphysical common thermal mapping;
    baryon loading constant and photon T(z) disagrees with observations).

Numerical toy baseline H(z) and historical S/S0=t/t0 are those used in earlier
SDMC radius calculations. This does not reconstruct microscopic F(S), vary z_*,
or run Planck spectra. No new H0 fit or candidate is claimed.
"""
import json, math, argparse
from pathlib import Path
from functools import lru_cache
from scipy.integrate import quad
from scipy.optimize import brentq

C=299792.458
OMEGA_GAMMA_PHYSICAL=2.469e-5
OMEGA_RAD_PHYSICAL=4.17998772e-5
PLANCK_100THETA_STAR=1.0411
Z_STAR=1089.2775

def main():
  a=argparse.ArgumentParser()
  a.add_argument("--out",type=Path,default=Path("output/structural_rb_h0_audit.json"))
  opts=a.parse_args()
  pars=json.loads((Path(__file__).parent/"late300_candidate.json").read_text())["parameters"]
  H_ACCEPTED=pars["H0"]
  WB=pars["omega_b"];WC=pars["omega_cdm"]
  W=(WB+WC); LAMBDA=pars["lambda_e"]; ZT=pars["z_t"]
  A=pars["A_late"];B=pars["B_late"]

  def E(z,H0=H_ACCEPTED):
    h=H0/100.;om=W/h**2;orr=OMEGA_RAD_PHYSICAL/h**2
    ox=1-om-orr
    xr=(orr*(1+z))/(om+orr*(1+z))
    n=-math.log1p(z);nt=-math.log1p(ZT)
    f=.5*(1-math.tanh((n-nt)/.5))*(3+xr)/LAMBDA**2
    delta=A*z*math.exp(-z/.25)-B*z*z*math.exp(-z/1.5)
    return math.sqrt((om*(1+z)**3+orr*(1+z)**4+ox)*(1+delta)**2/(1-f))

  @lru_cache(maxsize=256)
  def ageH0(z):
    # dimensionless cosmic age t(z)*H0 on the accepted H(z), integrated in a
    aa=1/(1+z)
    return quad(lambda a:1/(a*E(1/a-1)),0,aa,
                epsabs=1e-11,epsrel=1e-10,limit=300)[0]

  age0=ageH0(0.)
  def mapper(z):
    # S/S0 = t/t0, i.e. historical constant structural lapse approximation
    return (ageH0(z)/age0)*(1+z)
  def rho_struct_over_m(z):
    return (age0/ageH0(z))**2/(1+z)**3
  def Rb(z,form="physical"):
    ordinary=.75*(WB/OMEGA_GAMMA_PHYSICAL)/(1+z)
    if form=="physical":return ordinary
    if form=="forced_baryon":return ordinary*rho_struct_over_m(float(z))
    if form=="forced_both":return .75*(WB/OMEGA_GAMMA_PHYSICAL)
    raise ValueError(form)
  def cs(z,form):return C/math.sqrt(3*(1+Rb(z,form)))
  def acoustic(H0,form):
    DM=quad(lambda z:C/(H0*E(z,H0)),0,Z_STAR,
            epsrel=3e-8,limit=500)[0]
    rs=quad(lambda z:cs(z,form)/(H0*E(z,H0)),
            Z_STAR,1e8,epsrel=8e-8,limit=500)[0]
    return {"rs_Mpc":rs,"DM_Mpc":DM,"theta100":100*rs/DM}
  def theta_only_H0(form):
    return brentq(lambda H:acoustic(H,form)["theta100"]-PLANCK_100THETA_STAR,
           50,125,xtol=5e-8)

  z0=1090
  Zpoints=(0,20,30,1090,10000,1e7)
  history=[]
  for z in Zpoints:
    q=rho_struct_over_m(float(z))
    p=1/(E(z)*ageH0(float(z)))
    history.append({"z":z,"S_over_S0_chronological":ageH0(float(z))/age0,
     "structural_M_mapper":mapper(float(z)),
     "p_dlnS_dlna":p,"rho_struct_over_Jordan_m":q,
     "Rb_physical":Rb(z),"Rb_forced_baryons_only":Rb(z,"forced_baryon"),
     "baryon_number_source_Q_over_Hrho_if_forced":3-2*p,
     "photon_energy_source_Q_over_Hrho_if_forced":4-2*p})

  theta_accepted=acoustic(H_ACCEPTED,"physical")
  theta_forced=acoustic(H_ACCEPTED,"forced_baryon")
  H0_theta_physical=theta_only_H0("physical")
  H0_theta_forced=theta_only_H0("forced_baryon")
  tstar=ageH0(1090.)/age0
  rho_struct_over_rad_early=rho_struct_over_m(1e7)*W/(OMEGA_RAD_PHYSICAL*(1e7+1))
  coeff=9.74e53/(2.176434e-8*2.71678534e61)
  physical_rb1090=Rb(1090.)
  forced_rb1090=Rb(1090.,"forced_baryon")
  assert .614<physical_rb1090<.616
  assert 1.05<rho_struct_over_m(1090)<1.06
  assert .646<rho_struct_over_m(20)<.648
  assert 1.13<rho_struct_over_rad_early<1.15
  assert abs(theta_accepted["theta100"]-1.04164647)<2e-5
  assert 69.3<H0_theta_physical<69.8
  assert 91<H0_theta_forced<93
  result={
    "status":"PASSED; structural density is an algebraic reparameterization at z=0, not independently predictive",
    "accepted_Late300_H0_full_model_km_s_Mpc":H_ACCEPTED,
    "accepted_chronological_structural_history_not_action_derived":True,
    "M0_over_mpS0_rounded_PartII":coeff,
    "mass_relation_identity":"M0/(mP*S0)=(4*pi/3)*rho_m0/(rhoP/S0^2)",
    "Rb_standard_recomb_z1090":physical_rb1090,
    "cs_over_c_standard_recomb":1/math.sqrt(3*(1+physical_rb1090)),
    "Rb_forced_baryon_Sminus2_z1090":forced_rb1090,
    "Rb_forced_baryon_and_photon_Sminus2_all_redshifts":Rb(1090,"forced_both"),
    "cs_over_c_forced_both":1/math.sqrt(3*(1+Rb(1090,"forced_both"))),
    "T_recomb_K_physical":2.7255*1091,
    "T_recomb_K_forced_photons_Sminus2":2.7255/math.sqrt(tstar),
    "rho_struct_over_radiation_at_z1e7_if_matter_replaced":rho_struct_over_rad_early,
    "BBN_H_boost_if_extra_rhostr_source_against_radiation":
       math.sqrt(1+rho_struct_over_rad_early),
    "accepted_tracker_BBN_H_boost":1/math.sqrt(1-4/LAMBDA**2),
    "acoustic_at_fixed_accepted_H0":{
        "accepted_physical":theta_accepted,
        "forced_baryon_Sminus2_no_H_change":theta_forced
    },
    "theta_only_diagnostic_H0": {
       "Planck_target_100theta":PLANCK_100THETA_STAR,
       "with_accepted_baryon_loading":H0_theta_physical,
       "with_forced_baryon_loading":H0_theta_forced,
       "CAUTION":"Not a full CMB H0 inference: fixed recombination z*, target H(z) for density violation, unchanged perturbations."
    },
    "redshift_history":history,
    "verdict":"No new independent density constraint/H0 shift. Preserve separately conserved Jordan baryons and radiation unless a new covariant interaction is supplied. Forced histories fail conservation and high-z background."
  }
  opts.out.parent.mkdir(parents=True,exist_ok=True)
  opts.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
  print("PASS STRUCTURAL_MASS_IDENTITY C=",coeff)
  print("PASS RECOMB Rb",physical_rb1090,forced_rb1090)
  print("PASS ACOUSTIC H0 THETA_ONLY",H0_theta_physical,H0_theta_forced)
  print("PASS BBN_STRESS",rho_struct_over_rad_early,result["BBN_H_boost_if_extra_rhostr_source_against_radiation"])
  print("NO NEW ACCEPTED H0 OR LATE300 CANDIDATE")
if __name__=="__main__":main()
