#!/usr/bin/env python3
"""Era-spanning mass–stretch closure for accepted SDMC Late300.

This is a source-identity and physical-consistency investigation, NOT
a new solution of the field equations, new CMB/BBN likelihood, or fitter.

Model and epistemic separation:
  (1) M_G(0) from Manuscript Part II, present observed matter normalization;
  (2) hypothesize M_G(S) = k*m_P*S for ALL epochs, same k from (1);
  (3) calculate rho_G = (3*k/4pi)*rhoP/S² using the HISTORICAL
      chronological S/S0=t/t0, not an independently solved S equation;
  (4) compare this rho_G with separately conserved physical matter, radiation,
      and the full accepted Jordan-frame scalar-tensor Friedmann source;
  (5) show the interaction or source residual required by replacing physical
      matter, and show unchanged photon-baryon loading for a coordinate map;
  (6) DO NOT relabel any toy Friedmann inversion as a measured/derived H0.

Ordinary matter: rho_m0 a^-3; radiation: rho_r0 a^-4;
rho_X^J(z) = F(z)*rho_crit0*E(z)^2 - rho_m(z)-rho_r(z).
Normalized structural reference rho_P=m_P/l_P³; volume 4*pi/3 R³.
"""
import argparse
import json
import math
from pathlib import Path
from functools import lru_cache
from scipy.constants import G, pi
from scipy.integrate import quad
import numpy as np

MPC_M=3.085677581491367e22
MP=2.176434e-8
LP=1.616255e-35
S0=2.71678534e61
MG0_PARTII=9.74e53
C_KMS=299792.458
OMEGA_GAMMA_PHYS=2.469e-5
OMEGA_R_PHYS=4.17998772e-5
ZSTAR=1089.2775

class Audit:
    def __init__(self):
        self.p=json.loads((Path(__file__).parent/"late300_candidate.json").read_text())["parameters"]
        self.H0=self.p["H0"]
        self.h=self.H0/100
        self.omega_m=self.p["omega_b"]+self.p["omega_cdm"]
        self.om=self.omega_m/self.h**2
        self.orr=OMEGA_R_PHYS/self.h**2
        self.ox=1-self.om-self.orr
        self.H0_si=self.H0*1000/MPC_M
        self.rhoc0=3*self.H0_si**2/(8*pi*G)
        self.rhoP=MP/LP**3
        self.k=MG0_PARTII/(MP*S0)
        self.k_matter_today=(4*pi/3)*(self.rhoc0*self.om/self.rhoP)*S0**2
        self.R0=LP*S0
        self.rhoG0=MG0_PARTII/(4*pi/3*self.R0**3)
        self.rhom0=self.rhoc0*self.om
        self.rhor0=self.rhoc0*self.orr

    def E(self,z):
        a=1/(1+z); n=-math.log1p(z); nt=-math.log1p(self.p["z_t"])
        xr=self.orr*(1+z)/(self.om+self.orr*(1+z))
        f=.5*(1-math.tanh((n-nt)/.5))*(3+xr)/self.p["lambda_e"]**2
        delta=(self.p["A_late"]*z*math.exp(-z/.25)
               -self.p["B_late"]*z*z*math.exp(-z/1.5))
        return math.sqrt(
            (self.om*a**-3+self.orr*a**-4+self.ox)*(1+delta)**2/(1-f))

    @lru_cache(maxsize=1000)
    def age_H0(self,z):
        # Numerically integrate dimensionless t H0 over all Late300 eras.
        amax=1/(1+z)
        return quad(lambda a: 1/(a*self.E(1/a-1)),0,amax,
            epsabs=1e-11,epsrel=5e-10,limit=400)[0]

    def F(self,z):
        n=-math.log1p(z); nc=-math.log1p(self.p["z_c"])
        u=.5*(1+math.tanh((n-nc)/(2*self.p["width"])))
        return math.exp(self.p["A_F"]*u)

    def state(self,z):
        q=1+z
        a=1/q
        time_ratio=self.age_H0(z)/self.age_H0(0.)
        S=S0*time_ratio
        p=1/(self.E(z)*self.age_H0(z))
        fr=self.rhoG0*time_ratio**-2
        matter=self.rhom0*q**3
        radiation=self.rhor0*q**4
        photon=(self.rhoc0*OMEGA_GAMMA_PHYS/self.h**2)*q**4
        baryons=(self.rhoc0*self.p["omega_b"]/self.h**2)*q**3
        homogeneous_total=matter+radiation
        F=self.F(z)
        jordan_total=F*self.rhoc0*self.E(z)**2
        scalar=jordan_total-homogeneous_total
        x_rad=radiation/(matter+radiation)
        peff=(3+x_rad)/2
        OmG=fr/(self.rhoc0*self.E(z)**2)
        OmG_J=OmG/F
        p_exact_to_match=peff
        k_needed_mr=self.k*homogeneous_total/fr
        k_needed_matter=self.k*matter/fr
        source_Q_for_m=3-2*p
        source_Q_for_rad=4-2*p
        # Computed direct matter-domain stretch from physical matter density:
        Sm=(self.rhoP/matter)**(1/3)
        Km=S/Sm
        # Structural density ratio source decomposition; not fractions of a
        # total if sum >1 (at z=20 baryons+CDM > proposed structural rho).
        mass_eq_coeff=fr/matter
        rb=.75*baryons/photon
        grav_H_struct=math.sqrt((8*pi*G*fr)/(3*F))
        grav_H_sdmc=self.E(z)*self.H0_si
        # if rho_G replaces Jordan ordinary matter, with rho_r and scalar
        # fixed at the accepted background, it is NOT a full field solution:
        h_rhoG_replace_m=math.sqrt(max(0,(jordan_total-matter+fr)/jordan_total))
        # if rho_G replaces all Jordan sources:
        h_rhoG_replace_all=grav_H_struct/grav_H_sdmc
        return {
            "z":float(z),"S_over_S0":time_ratio,"S":S,"p_chronological":p,
            "p_needed_if_G_equals_matter_plus_radiation":p_exact_to_match,
            "p_mismatch":p-p_exact_to_match,
            "F":F,
            "R_struct_Mpc":LP*S/MPC_M,
            "rho_G_kg_m3":fr,"rho_m_kg_m3":matter,
            "rho_r_kg_m3":radiation,"rho_XJ_kg_m3":scalar,
            "rho_total_J_kg_m3":jordan_total,
            "rhoG_over_m":mass_eq_coeff,
            "rhoG_over_m_plus_r":fr/homogeneous_total,
            "rhoG_over_full_Jordan":fr/jordan_total,
            "rhoG_over_r":fr/radiation,
            "k_required_for_full_m_plus_r":k_needed_mr,
            "k_required_for_matter_only":k_needed_matter,
            "Omega_G_standard":OmG,"Omega_G_Jordan":OmG_J,
            "H_structural_only_over_H_Late300":h_rhoG_replace_all,
            "H_trial_replace_matter_over_H_Late300":h_rhoG_replace_m,
            "Rb_conserved":rb,
            "Rb_structurally_reexpressed":.75*(baryons/fr)/(photon/fr),
            "Sm_independent_cubic":Sm,"Km_matter_mapper":Km,
            "Qm_required_if_matter_tracks_Sminus2_over_Hrho":source_Q_for_m,
            "Qr_required_if_radiation_tracks_Sminus2_over_Hrho":source_Q_for_rad,
            "w_of_standalone_conserved_structural_fluid":-1+2*p/3
        }

    def result(self):
        zpoints=[0,.5,3.927876388467848,7,10,16.19317622240633,20,30,
           100,1090,3400,10000,1e7,1e9]
        rows=[self.state(z) for z in zpoints]
        # Analytic simplification using units: rhoP=(c²/G) lP^-2
        # -> rhoG = (3*k/(4pi))c²/(G R²)
        identity=[]
        for r in rows:
            z=r["z"];E=self.E(z);p=r["p_chronological"]
            S=S0*r["S_over_S0"]
            lapse=LP*S/(299792458.*(self.age_H0(z)/self.H0_si))
            # Historically the numerical structural lapse implied by S/t.
            lhs=r["Omega_G_standard"]
            rhs=2*self.k*p*p/lapse**2*(G*MP/(LP*299792458.**2))
            identity.append(abs(lhs-rhs))
        maxid=max(identity)
        assert maxid<1e-11, maxid
        assert abs(rows[0]["rhoG_over_m"]-(self.rhoG0/self.rhom0))<1e-11
        assert .63<self.state(20)["rhoG_over_m"]<.67
        assert 1.03<self.state(1090)["rhoG_over_m"]<1.1
        assert 1.1<self.state(1e7)["rhoG_over_r"]<1.2
        assert .28<self.state(0)["rhoG_over_full_Jordan"]<.31
        assert .78<self.state(1090)["rhoG_over_m_plus_r"]<.83
        assert .1<abs(self.state(3400)["p_mismatch"]) or abs(self.state(3400)["p_mismatch"])>.02
        # All epochs must preserve the particle-number physical baryon loading
        assert max(abs(row["Rb_conserved"]-row["Rb_structurally_reexpressed"])
                   for row in rows)<1e-12
        assert abs(self.state(1090)["Rb_conserved"]-.61486221617)<1e-8
        # Full source should not be claimed without stress-energy split.
        assert self.state(20)["k_required_for_full_m_plus_r"] > 2.4
        assert self.state(1090)["k_required_for_full_m_plus_r"] > 2
        assert self.state(1e7)["k_required_for_full_m_plus_r"] < 1.5
        # Show direct mass-conservation conflict if M_G itself is ordinary mass.
        ratio_MG_recomb=self.state(1090)["S_over_S0"]
        assert ratio_MG_recomb<1e-4
        coeffG0=self.k
        two_mass_sources={
          "proposed_epoch_mass_MG_over_MG0":ratio_MG_recomb,
          "ordinary_observable_patch_conserved_mass_ratio":1.,
          "postulated_MG_at_S1_mP":coeffG0,
          "boundary_condition_if_one_planck_mass_requires_k_at_S1":1.
        }
        out={
         "scope":"Algebraic structural-mass closure, historical chronological S approximation, conserved Jordan fluid and covariant-source diagnostics. No new field-action fit.",
         "accepted_Late300_H0":self.H0,
         "Planck_length_m":LP,"Planck_mass_kg":MP,"Planck_density_kg_m3":self.rhoP,
         "S0_chronological_reference":S0,
         "M_G_present_partII_kg":MG0_PARTII,
         "coefficient_k_from_partII":self.k,
         "coefficient_k_needed_to_match_Late300_matter_today":self.k_matter_today,
         "present_structural_mass_normalization_relative_to_accepted_matter_density":self.rhoG0/self.rhom0,
         "structural_density_coefficient_relative_rhoP":3*self.k/(4*pi),
         "present_structural_density_kg_m3":self.rhoG0,
         "present_matter_density_kg_m3":self.rhom0,
         "max_analytic_OmegaG_vs_2k_p2_over_N2_error":maxid,
         "formulae":{
           "mass":"M_G(S)=k*mP*S (hypothesis at epochs other than z=0)",
           "density":"rho_G=(3*k/(4*pi))*rhoP/S**2",
           "omegaG":"rho_G/rho_crit_ref=2*k*[c/(H*Rstruct)]**2=2*k*p**2/lapse**2",
           "structural_era_exponent":"n_G=2*p",
           "K_m_if_rhoG_identified_with_m":"Km=(3*k/(4*pi))**(1/3)*S**(1/3)",
           "Km_actual_independent_matter":"Km=S/(rhoP/rhom)**(1/3)",
           "dlog_Km_dlog_a_if_conserved_matter":"p-1",
           "density_to_critical_Jordan":"rhoG/(F*rhoCrit)",
           "structural_as_separately_conserved_fluid":"w_G=-1+2*p/3 (conditional, no interaction)",
           "required_era_coefficient":"k_req(z)=k*[rho_m(z)+rho_r(z)]/rho_G(z)",
           "Rb_exact":"3/4*(rho_b/rho_G)/(rho_gamma/rho_G)=3/4*(omega_b/omega_gamma)/(1+z)",
           "mass_source_vs_patches":"M_G(S)/M_G0=S/S0, versus same mass in conserved comoving patch."
         },
         "same_mass_consistency":two_mass_sources,
         "BBN_reference":{
           "z":1e7,
           "early_rhoG_to_radiation":self.state(1e7)["rhoG_over_r"],
           "early_H_ratio_if_rhoG_replaces_all":self.state(1e7)["H_structural_only_over_H_Late300"],
           "early_H_ratio_if_rhoG_replaces_m_but_radiation_also_added":self.state(1e7)["H_trial_replace_matter_over_H_Late300"],
           "accepted_early_H_speedup":1/math.sqrt(1-4/self.p["lambda_e"]**2),
           "warning":"Some experiments double count existing radiation; none is an admissible full BBN model."
         },
         "rows":rows,
         "conclusion":"Current k is derived from present matter normalization; universal fixed-k structural density is not the separately conserved baryon/radiation/total Jordan source under historical S; it cannot give a new Rb or H0 without an independent interaction/action."
        }
        return out

def plot(out,folder):
    import matplotlib;matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder.mkdir(parents=True,exist_ok=True)
    a=Audit()
    zs=np.logspace(-3,7,120)
    zs=np.concatenate(([0.],zs))
    r=[a.state(float(z)) for z in zs]
    x=np.log10(1+zs)
    fig,ax=plt.subplots(figsize=(9,5.25))
    for k,label in [("rhoG_over_m","Structural / conserved matter"),
                    ("rhoG_over_m_plus_r","Structural / (matter + radiation)"),
                    ("rhoG_over_full_Jordan","Structural / full Jordan source")]:
      ax.plot(x,[v[k] for v in r],label=label,linewidth=2)
    ax.axhline(1,linewidth=1,linestyle="--",color="black")
    for z,label in [(20,"z=20"),(1090,"recomb."),(3400,"equality")]:
      ax.axvline(np.log10(1+z),color="grey",alpha=.35,linestyle=":")
    ax.set(xlabel="log₁₀(1+z)",ylabel="Density ratio",
       title="One fixed structural mass law cannot represent every conserved fluid")
    ax.set_ylim(0,5);ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.tight_layout();fig.savefig(folder/"density_closure.png",dpi=160);plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,5.25))
    ax.plot(x,[v["p_chronological"] for v in r],label="Late300 chronological p",linewidth=2)
    ax.plot(x,[v["p_needed_if_G_equals_matter_plus_radiation"] for v in r],
      label="p required to reproduce matter + radiation",linewidth=2,linestyle="--")
    ax.set(xlabel="log₁₀(1+z)",ylabel="Structural scaling p",
           title="Era-dependent exponent: required versus chronological structural history")
    ax.grid(alpha=.2);ax.legend();fig.tight_layout()
    fig.savefig(folder/"epoch_structural_clock.png",dpi=160);plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,5.25))
    ax.plot(x,[v["k_required_for_full_m_plus_r"] for v in r],
     linewidth=2,label="k required to match matter+radiation")
    ax.axhline(out["coefficient_k_from_partII"],linestyle="--",color="black",
     label="constant k from Part II present mass")
    ax.set(xlabel="log₁₀(1+z)",ylabel="Structural mass coefficient k",
     title="Testing epoch invariance of M_G/(mP S)")
    ax.grid(alpha=.2);ax.legend();fig.tight_layout()
    fig.savefig(folder/"mass_coefficient_closure.png",dpi=160);plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("output/structural_mass_full_closure.json"))
    p.add_argument("--plots",type=Path,default=Path("output/structural_mass_full_closure_plots"))
    v=p.parse_args();a=Audit();data=a.result()
    v.out.parent.mkdir(parents=True,exist_ok=True)
    v.out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    plot(data,v.plots)
    print("SUCCESS UNIVERSAL-MASS ERA CLOSURE REPORT")
    print("COEFF",data["coefficient_k_from_partII"],"NORMALIZED_TO_M",data["coefficient_k_needed_to_match_Late300_matter_today"])
    print("OMEGA_IDENTITY_MAX_ERR",data["max_analytic_OmegaG_vs_2k_p2_over_N2_error"])
    for d in data["rows"]:
        print("MASS_ERA_ROW", json.dumps({k:round(d[k],8) for k in [
        "z","p_chronological","p_needed_if_G_equals_matter_plus_radiation",
        "rhoG_over_m","rhoG_over_m_plus_r","rhoG_over_full_Jordan",
        "k_required_for_full_m_plus_r","Omega_G_Jordan",
        "H_structural_only_over_H_Late300","H_trial_replace_matter_over_H_Late300",
        "Km_matter_mapper","Rb_conserved"]},sort_keys=True),flush=True)
    print("NO NEW H0 NO VALID CMB OR BBN REOPTIMIZATION")

if __name__=="__main__":main()
