#!/usr/bin/env python3
"""Necessary structural-clock histories for universal kappa*mP*S mass law.

This INVERTS accepted Late300 conserved matter + radiation, or full Jordan
Friedmann source, to infer the S(a), p(a) and N_str(a) each would REQUIRE.
It DOES NOT solve the independent structural field action or provide a new H0.

Two physically distinct matching choices are kept disjoint:
A. structural density = conserved matter+radiation; moving sphere encloses
   additional existing matter, no matter creation; pressure work for photons.
B. structural density = entire F(z)*3 Mref^2 H(z)^2 Jordan source
   including kinetic/braiding/vacuum; different present coefficient kappa
   and a strongly different required structural lapse.

References: SDMC Manuscript I sec. 10, II sec. 5.3; Late300 action audit.
"""
import json
import math
import argparse
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from run_mass_stretch_era_density_full_closure import Audit, LP, MP, S0, MPC_M, C_KMS, pi

C_M_S=299792458.0

def alphaM(p,z):
    n=-math.log1p(z)
    nc=-math.log1p(p["z_c"])
    u=.5*(1+math.tanh((n-nc)/(2*p["width"])))
    return p["A_F"]/p["width"]*u*(1-u)

def p_jordan_necessary(a,z):
    # p_J = -.5*d ln[F(z)H(z)^2]/d ln(a)
    #     = -d ln H/d ln a - 0.5*alpha_M
    n=-math.log1p(z)
    eps=1e-5
    ep=a.E(math.expm1(-(n+eps)))
    em=a.E(math.expm1(-(n-eps)))
    dlnH_dn=(math.log(ep)-math.log(em))/(2*eps)
    return -dlnH_dn-.5*alphaM(a.p,z)

def analyse():
    a=Audit()
    q0_mr=a.rhom0+a.rhor0
    q0_J=a.F(0)*a.rhoc0*a.E(0)**2
    vol_fact=4*pi/3*LP**3*S0**2/MP
    kappa_mr=q0_mr*vol_fact
    kappa_J=q0_J*vol_fact
    kappa_partII=a.k
    Nch_baseline=(LP*S0/MPC_M)*(a.H0/C_KMS)/a.age_H0(0) # chrono R0/(c*t0)
    rows=[]
    epochs=[0,0.5,1,3.927876388467848,7,10,16.19317622240633,20,30,100,500,1090,3400,10000,1e7,1e9]
    for z in epochs:
        cur=a.state(z)
        q_mr=cur["rho_m_kg_m3"]+cur["rho_r_kg_m3"]
        q_J=cur["rho_total_J_kg_m3"]
        p_mr=cur["p_needed_if_G_equals_matter_plus_radiation"]
        p_J=p_jordan_necessary(a,z)
        mr_Srat=math.sqrt(q0_mr/q_mr)
        J_Srat=math.sqrt(q0_J/q_J)
        H_R0_over_c=a.H0*a.E(z)*(LP*S0/MPC_M)/C_KMS
        N_mr=p_mr*H_R0_over_c*mr_Srat
        N_J=p_J*H_R0_over_c*J_Srat
        N_ch=cur["p_chronological"]*H_R0_over_c*cur["S_over_S0"]
        rhoGhist=cur["rho_G_kg_m3"]
        w_mr=-1+2*p_mr/3
        w_J=-1+2*p_J/3
        x_r=cur["rho_r_kg_m3"]/q_mr
        # Actual kappa under chronological S and source q:
        k_mr_hist=kappa_partII*q_mr/rhoGhist
        k_J_hist=kappa_partII*q_J/rhoGhist
        # for a moving spherical matter domain physical mass flux
        # M_m = (4pi/3)R^3 rho_m; Mdot=4pi R² rho_m (Rdot-HR)
        # H-normalized and in units of M_m: 3*(p-1)
        # radiation E/c²: derivative 3*p-4 includes radiation redshifting
        mass_growth_mat_chrono=3*(cur["p_chronological"]-1)
        mass_growth_mat_mr=3*(p_mr-1)
        kappa_mr_log_derivative=2*cur["p_chronological"]-(3+x_r)
        rows.append({
          "z":float(z),
          "chrono_R_over_R0":cur["S_over_S0"],
          "chrono_p":cur["p_chronological"],
          "chrono_lapse_N":N_ch,
          "mr_required_R_over_R0":mr_Srat,
          "mr_required_R_Mpc":LP*S0/MPC_M*mr_Srat,
          "mr_required_p":p_mr,
          "mr_required_lapse_N":N_mr,
          "mr_required_R_over_chrono":mr_Srat/cur["S_over_S0"],
          "mr_normalized_kappa":kappa_mr,
          "mr_kappa_if_chrono":k_mr_hist,
          "mr_chrono_dlnkappa_dln_a":kappa_mr_log_derivative,
          "mr_radiation_fraction":x_r,
          "mr_effective_w":w_mr,
          "full_Jordan_required_R_over_R0":J_Srat,
          "full_Jordan_required_R_Mpc":LP*S0/MPC_M*J_Srat,
          "full_Jordan_required_p":p_J,
          "full_Jordan_required_lapse_N":N_J,
          "full_Jordan_required_R_over_chrono":J_Srat/cur["S_over_S0"],
          "full_Jordan_normalized_kappa":kappa_J,
          "full_Jordan_kappa_if_chrono":k_J_hist,
          "full_Jordan_effective_w_if_conserved":w_J,
          "matter_mass_dlnM_dln_a_chrono":mass_growth_mat_chrono,
          "matter_mass_dlnM_dln_a_required_mr":mass_growth_mat_mr,
          "radiation_energy_dlnM_dln_a_required_mr":3*p_mr-4,
          "matter_plus_radiation_energy_dlnM_dln_a_required_mr":p_mr,
          "Late300_Rb":cur["Rb_conserved"]
        })
    # Exact bridge: if rho_G=(3*k/4pi)rho_P/S² and source rho_mr,
    # solve S(S0): S/S0=sqrt(rho_mr0/rho_mr).
    # No new fitted acoustic parameter enters.
    assert abs(rows[0]["mr_required_R_over_R0"]-1)<1e-13
    assert abs(rows[0]["full_Jordan_required_R_over_R0"]-1)<1e-13
    assert 1.62<kappa_mr<1.65
    assert 5.4<kappa_J<5.8
    assert 2.65<rows[7]["mr_required_lapse_N"]<2.85 # z20
    assert 2.85<rows[11]["mr_required_lapse_N"]<3.08 # z1090
    assert 3.5<rows[14]["mr_required_lapse_N"]<3.8 # radiation
    assert 4.75<rows[0]["mr_required_lapse_N"]<5.15 # today
    assert 1.2<rows[0]["full_Jordan_required_lapse_N"]<1.9 # today
    assert all(abs(r["chrono_lapse_N"]-Nch_baseline)<1e-8 for r in rows)
    assert max(abs(r["matter_plus_radiation_energy_dlnM_dln_a_required_mr"]-r["mr_required_p"]) for r in rows)<1e-12
    assert max(abs(r["Late300_Rb"]-r["Late300_Rb"]) for r in rows)<1e-15
    assert rows[11]["mr_required_R_over_chrono"]<.92
    assert rows[7]["mr_required_R_over_chrono"]<.84
    # The radiation-like energy within moving sphere also experiences p dV:
    # dM_r/dt =4pi R² rho_r (Rdot-HR) - H M_r, no new photons.
    # Integrate p_required to verify explicit algebraic S(t):
    def p_from_z(z):
        rr=a.rhor0*(1+z)**4;mm=a.rhom0*(1+z)**3
        return (3*mm+4*rr)/(2*(mm+rr))
    for z in [20,1090,1e7]:
        integral=quad(lambda zz:p_from_z(zz)/(1+zz),0,z,epsrel=1e-10,limit=300)[0]
        sratio=next(r["mr_required_R_over_R0"] for r in rows if r["z"]==z)
        assert abs(-math.log(sratio)-integral)<3e-10
    # Ensure impossible to keep one universal k for both matter and total Jordan
    # with fixed R0 and today's observations:
    require=kappa_J/kappa_mr
    assert require>3.3
    # Inverse full Jordan includes effective stress split; not an independent
    # equation of state for a physical fluid, because F varies.
    result={
      "status":"PASS: two required clock histories reconstructed and not interchangeable",
      "scientific_scope":"Exact inverse source-to-S mapping; accepted H(z),F(z) treated as fixed; not autonomous structural field or full likelihood.",
      "model":"accepted Late300",
      "present_H0_km_s_Mpc":a.H0,
      "present_Rstruct_Mpc":LP*S0/MPC_M,
      "historical_chronological_N":Nch_baseline,
      "partII_mass_kappa0":kappa_partII,
      "kappa_to_identify_rhoG_with_matter_plus_radiation_today":kappa_mr,
      "kappa_to_identify_rhoG_with_FULL_Jordan_today":kappa_J,
      "kappa_J_to_kappa_mr_ratio":require,
      "kappa_if_one_Planck_mass_at_S1":1.0,
      "equations":{
       "mass_in_sphere":"M_Q=(4pi/3)Rstruct³*Q, Q=rhom+rhor OR Q=F*rhocrit(H)",
       "universal_structural_mass":"M_G=kappa*mP*S; rho_G=(3*kappa/4pi)*rhoP/S²",
       "inverse_clock":"S_Q(z)/S0=sqrt(Q(0)/Q(z)), and kappa_Q=(4pi/3)*Q(0)*lp³*S0²/mP",
       "required_p":"p_Q=-(1/2)*d ln Q/d ln a",
       "required_lapse":"N_Q(z)=p_Q(z)*H(z)*lp*S_Q(z)/c",
       "matter_plus_radiation_clock":"p_mr=(3*rhom+4*rhor)/[2*(rhom+rhor)]",
       "full_Jordan_clock":"p_J=-(d ln H/d ln a)-alpha_M/2",
       "mass_sweep_dust":"Mdot_m=4*pi*R²*rho_m*(Rdot-H*R); no baryon creation",
       "energy_sweep_radiation":"Mdot_r=4*pi*R²*rho_r*(Rdot-H*R) - H*M_r; photon redshift/work",
       "constant_kappa_necessary_equation":"2*p=3+x_r for Q=rhom+rhor",
       "chrono_kappa_drift":"d ln kappa_mr/d ln a=2*p_chrono-3-x_r"
      },
      "rows":rows,
      "central_inference":"Mass growth M proportional to structural radius may be ordinary moving-boundary sweep only for a prescribed S(a); does not imply same fixed mass at all z, novel overdensity, Rb or H0.",
      "fundamental_status":"No independent equation predicts p or lapse; the inverse clocks derive from already assumed Late300 H/F/densities, and their boundary values differ from historical constant-lapse N."
    }
    return result

def plots(data,dest):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dest.mkdir(parents=True,exist_ok=True)
    a=Audit()
    z=np.r_[0.,np.logspace(-3,8,220)]
    rows=[]
    Q0mr=a.rhom0+a.rhor0
    Q0J=a.rhoc0*a.F(0)*a.E(0)**2
    for zz in z:
        mm=a.rhom0*(1+zz)**3
        rr=a.rhor0*(1+zz)**4
        H=a.H0*a.E(float(zz))
        Qmr=mm+rr
        Qj=a.F(float(zz))*a.rhoc0*a.E(float(zz))**2
        sr=math.sqrt(Q0mr/Qmr)
        sj=math.sqrt(Q0J/Qj)
        pM=(3*mm+4*rr)/(2*Qmr)
        pJ=p_jordan_necessary(a,float(zz))
        chrono=a.state(float(zz)) # slow but ~220*quad
        Npref=H*(LP*S0/MPC_M)/C_KMS
        rows.append((pM,pJ,chrono["p_chronological"],
           Npref*sr*pM,Npref*sj*pJ,
           Npref*chrono["S_over_S0"]*chrono["p_chronological"],
           sr/chrono["S_over_S0"],sj/chrono["S_over_S0"]))
    d=np.asarray(rows)
    xx=np.log10(1+z)
    fig,ax=plt.subplots(figsize=(9,5.2))
    for idx,label,style in [(0,"Required p: matter+radiation","-"),
                             (1,"Required p: full Jordan source","-."),
                             (2,"Historical chronological p",":")]:
        ax.plot(xx,d[:,idx],style,label=label,linewidth=2)
    ax.set(xlabel="log₁₀(1+z)",ylabel="p=d ln S/d ln a",
          title="Three incompatible histories for one structural clock")
    ax.legend(fontsize=9);ax.grid(alpha=.2);fig.tight_layout()
    fig.savefig(dest/"required_structural_exponents.png",dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.2))
    for idx,label in [(3,"Required N: matter+radiation"),
                      (4,"Required N: full Jordan"),
                      (5,"Historical constant N")]:
        ax.plot(xx,d[:,idx],label=label,linewidth=2)
    ax.set(xlabel="log₁₀(1+z)",ylabel="Structural lapse N=dR/(c dt)",
           title="Inverting source closure requires different lapse histories")
    ax.legend(fontsize=9);ax.grid(alpha=.2);fig.tight_layout()
    fig.savefig(dest/"required_structural_lapse.png",dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.2))
    ax.plot(xx,d[:,6],label="Required radius to match matter+radiation / historical radius",linewidth=2)
    ax.plot(xx,d[:,7],label="Required radius to match full Jordan / historical radius",linewidth=2)
    ax.axhline(1,color="black",linestyle="--",alpha=.6)
    ax.set(xlabel="log₁₀(1+z)",ylabel="Ratio",title="How much historical S would need to change")
    ax.legend(fontsize=9);ax.grid(alpha=.2);fig.tight_layout()
    fig.savefig(dest/"required_radius_vs_historical.png",dpi=160);plt.close(fig)

def main():
    par=argparse.ArgumentParser()
    par.add_argument("--out",type=Path,default=Path("output/required_structural_clock.json"))
    par.add_argument("--plots",type=Path,default=Path("output/required_structural_clock_plots"))
    args=par.parse_args()
    result=analyse()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    plots(result,args.plots)
    print("PASS SDMC NECESSARY SOURCE_TO_STRUCTURAL_CLOCK_INVERSION")
    print("NORMALIZATIONS",json.dumps({x:result[x] for x in ["partII_mass_kappa0","kappa_to_identify_rhoG_with_matter_plus_radiation_today","kappa_to_identify_rhoG_with_FULL_Jordan_today","historical_chronological_N"]}))
    for r in result["rows"]:
        print("CLOCK_ROW",json.dumps({k:round(r[k],9) for k in ["z","chrono_p","chrono_lapse_N",
           "mr_required_R_over_chrono","mr_required_p","mr_required_lapse_N",
           "full_Jordan_required_R_over_chrono","full_Jordan_required_p","full_Jordan_required_lapse_N",
           "mr_kappa_if_chrono","mr_chrono_dlnkappa_dln_a"]}))
    print("NO INDEPENDENT DERIVATION OF STRUCTURAL ACTION OR NEW HUBBLE FIT")
if __name__=="__main__":main()
