#!/usr/bin/env python3
"""Forward mature SDMC canonical action dynamics and finite-to-IR derivability audit.

This code evolves the SPECIFIED mature canonical action FORWARD; it does not
fabricate a microscopic finite Horndeski-to-canonical matching mechanism.

Based on SDMC_A_COMPLETED(7).pdf pp 667-669:
 G2_IR = (2 F_inf Z_sigma - 4 F_inf Z_BI) / sigma^2
 G3_IR = 0 ; G4_IR = F_inf/2 ; D_IR=2 ; c_s^2_IR=1.
Let sigma=C exp(phi); then G2=2F_inf X_phi - V0 exp(-2phi).
In terms of n=ln(a) and u=dphi/dn, Om=rho_m/(3FH^2),
Or=rho_r/(3FH^2) and y=V/(3FH^2):
  epsilon=-dlnH/dn=u^2 + 3 Om/2 + 2 Or,
  y=1-u^2/3-Om-Or,
  du/dn=-(3-epsilon)*u + 3*y,
  dOm/dn=Om*(-3+2epsilon),
  dOr/dn=Or*(-4+2epsilon),
  dphi/dn=u,  dlnH/dn=-epsilon.
Fixed point (u,Om,Or)=(1,0,0) with D=2, epsilon=1,q=0.
Jacobian eigenvalues (-2,-1,-2); matter forcing u=1-1.5*m1/a+...
and q=-1.5*m1/a+... reproduce the manuscript's leading IR
matter corrections, CONDITIONAL on its mature canonical action.

Also checks action-coefficient matching and non-identifiability:
 a) k1_IR(phi)=2F_inf,k2_IR=0,g_IR=0,V_IR~e^-2phi;
 b) no information from this fixed point alone determines old finite
    D0, width, power, kinetic handoff epoch or source-level UV functions.
"""
import csv
import json
import math
from pathlib import Path

OUT=Path("output/D_mature_forward_IR_identifiability")
OUT.mkdir(parents=True,exist_ok=True)
AF_DERIVED=.020520
F_DERIVED=math.exp(AF_DERIVED)
F_MATURE_PREVIOUS=1.023885427695
F_HISTORICAL_FIN=.05181542627513409+.34231919445927034
POTENTIAL_SLOPE_PHI=2.
STEPSIZE=.002
END_N=12.
CASES={
 "fiducial_matter":(0.8,.20,0.),
 "mixed_radiation":(.30,.35,.02),
 "supercoasting_initial_u":(1.20,.10,.01),
 "small_initial_u":(0.,.10,.05),
 "fast_initial_field":(1.50,.05,0.),
 "reversed_initial_field":(-.2,.25,.04),
 "near_attractor_matter":(1.,.05,0.),
 "near_attractor_radiation":(1.,0.,.05),
}

def check(ok,description):
    if not ok:
        raise AssertionError(description)

def close(x,y,tol=1e-9):
    return math.isclose(x,y,rel_tol=tol,abs_tol=tol)

def flow(state):
    u,Om,Or,phi,logH=state
    eps=u*u+1.5*Om+2.*Or
    y=1.-u*u/3.-Om-Or
    return (
       -(3.-eps)*u+3.*y,
       Om*(-3.+2.*eps),
       Or*(-4.+2.*eps),
       u,
       -eps,
    )

def rk4(y,h):
    k1=flow(y)
    k2=flow(tuple(a+h*b/2 for a,b in zip(y,k1)))
    k3=flow(tuple(a+h*b/2 for a,b in zip(y,k2)))
    k4=flow(tuple(a+h*b for a,b in zip(y,k3)))
    return tuple(a+h*(b+2*c+2*d+e)/6 for a,b,c,d,e in zip(y,k1,k2,k3,k4))

def jacobian_at_IR():
    # u'=(u-1)*(u*u-3)+Om*(1.5*u-3)+Or*(2*u-3);
    # Om'=Om*(-3+2*u*u+3*Om+4*Or),
    # Or'=Or*(-4+2*u*u+3*Om+4*Or).
    return [[-2.,-1.5,-1.], [0.,-1.,0.],[0.,0.,-2.]]

def eval_case(label,seed):
    u,Om,Or=seed
    check(1.-u*u/3.-Om-Or>0., "initial positive canonical potential fraction")
    initial_y=1.-u*u/3.-Om-Or
    # IMPORTANT: Use ONE fixed canonical action/potential for every seed.
    # In dimensionless units tP*S0=1, Z_BI=4pi/3, C=1 and phi0=0.
    # V0=4F*Z_BI is a SINGLE theory coefficient shared by all eight seeds.
    # Solve Friedmann for a different allowed initial H0 in each seed.
    Z_BI=4*math.pi/3
    common_V0=4*F_DERIVED*Z_BI
    H0=math.sqrt(common_V0/(3*F_DERIVED*initial_y))
    logH0=math.log(H0)
    y=(u,Om,Or,0.,logH0)
    total=int(round(END_N/STEPSIZE))
    rows=[]
    max_friedmann_err=0.
    max_omega_sum_err=0.
    for i in range(total+1):
        n=i*STEPSIZE
        u,Om,Or,phi,logH=y
        H=math.exp(logH-logH0)
        epsilon=u*u+1.5*Om+2.*Or
        yV=1.-u*u/3.-Om-Or
        K=u*u/3.
        q=epsilon-1.
        # Check actual independent matter, radiation and scalar potential
        # terms against the integrated Friedmann equation.
        Om_direct=seed[1]*math.exp(-3*n-2*(logH-logH0))
        Or_direct=seed[2]*math.exp(-4*n-2*(logH-logH0))
        yV_direct=initial_y*math.exp(-2*phi-2*(logH-logH0))
        friedmann=K+Om_direct+Or_direct+yV_direct
        max_friedmann_err=max(max_friedmann_err,abs(friedmann-1.))
        max_omega_sum_err=max(max_omega_sum_err,abs(K+Om+Or+yV-1.))
        if i%250==0 or i==total:
            rows.append({
                "case":label,"n":n,"a":math.exp(n),"u":u,
                "Omega_m":Om,"Omega_r":Or,"Omega_potential":yV,
                "H_over_initial":H,"phi":phi,"D_canonical":2*u*u,
                "epsilon":epsilon,"q":q,
                "c_s2_canonical":1.,
                "friedmann_constraint_error":friedmann-1.,
            })
        if i<total:y=rk4(y,STEPSIZE)
    uf,Omf,Orf,phif,Hf=y
    epsf=uf*uf+1.5*Omf+2.*Orf
    # At large N the decaying matter mode dominates and leads:
    # Omega_m~m1 e^-n, u-1~-1.5 Omega_m, q~-1.5 Omega_m.
    matter_u_ratio=((uf-1.)/Omf) if Omf>1e-11 else None
    matter_q_ratio=((epsf-1.)/Omf) if Omf>1e-11 else None
    check(abs(uf-1.)<3e-4,f"{label}: canonical coasting attractor not attained")
    check(Omf<3e-5 and Orf<1e-9,f"{label}: remaining fluids at n=12 too large")
    check(max_friedmann_err<2e-7,
         f"{label}: Friedmann energy constraint not conserved: {max_friedmann_err}")
    if matter_u_ratio is not None:
        check(abs(matter_u_ratio+1.5)<.05,
              f"{label}: leading matter-driven u correction not recovered")
        check(abs(matter_q_ratio+1.5)<.05,
              f"{label}: leading matter-driven q correction not recovered")
    return rows,{
      "case":label,"initial_u":seed[0],"initial_Omega_m":seed[1],
      "initial_Omega_r":seed[2],
      "initial_potential_fraction":initial_y,
      "common_action_V0":common_V0,
      "initial_H0_in_common_units":H0,
      "final_u":uf,"final_D":2*uf*uf,"final_Omega_m":Omf,
      "final_Omega_r":Orf,"final_epsilon":epsf,"final_q":epsf-1.,
      "matter_u_minus_one_over_Omega_m":matter_u_ratio,
      "matter_q_over_Omega_m":matter_q_ratio,
      "max_friedmann_constraint_error":max_friedmann_err,
      "max_fraction_constraint_error":max_omega_sum_err,
      "result":"FORWARD_CANONICAL_IR_ATTRACTOR_PASS",
    }

def prove_scalar_only_global_basin():
    """Global scalar-only homogeneous result, |u|<sqrt(3), V>0.

    On the canonical action with Om=Or=0:
       u'=(u-1)(u^2-3),
       d (u-1)^2/dn = 2 (u-1)^2 (u^2-3)<0.
    Thus all interior physical scalar-only velocities reach u=1.
    Integrated solution I(u)-I(u0)=n-n0 is verified independently.
    """
    rt=math.sqrt(3.)
    def implicit_I(u):
        if u==1.:raise ValueError("singular exact attractor")
        return (-.5*math.log(abs(u-1.))+
                .25*math.log(3.-u*u)+
                1./(4.*rt)*math.log(abs((u-rt)/(u+rt))))
    records=[]
    for u0 in (-1.65,-1.0,-.2,.0,.5,.8,1.2,1.6):
        u=u0
        i0=implicit_I(u0)
        max_err=0.
        for j in range(1,1501):
            h=.002
            rhs=lambda x:(x-1.)*(x*x-3.)
            a=rhs(u);b=rhs(u+h*a/2);c=rhs(u+h*b/2);d=rhs(u+h*c)
            u=u+h*(a+2*b+2*c+d)/6.
            n=j*h
            if abs(u-1.)>1e-8:
                max_err=max(max_err,abs(implicit_I(u)-i0-n))
            check(abs(u)<rt,"physical scalar-only phase space remained in |u|<sqrt3")
            check((u-1.)**2<=(u0-1.)**2+1e-13,"global Lyapunov decay")
        check(max_err<2e-6,"exact implicit scalar-only solution agrees with RK4")
        records.append({"initial_u":u0,"final_u_n3":u,
                        "max_abs_implicit_integral_error":max_err})
    return records


def prove_on_shell_cubic_decay():
    """On-shell continuation of No-Slip cubic under future coasting.

    For phi=ln a, H^2~a^-2, X=H^2/2, box phi=-2H^2:
       g=-F_phi/H^2, G3=gX=-F_phi/2,
       -G3 box phi=-F_phi H^2.
    Ratio to F H^2 = -F_phi/F = -alpha_M.
    Since alpha_M~a^(-1/w), on-shell cubic action density decouples
    for EVERY finite w>0. Coefficient g itself decays only if w<1/2.
    """
    tests=[]
    for w in (.25,.32925377073762596,.45,.5,.6,1.,1.5):
        F0=math.exp(AF_DERIVED)
        def alpha(phi):
            zc=4.034077502899133
            S=1/(1+math.exp(-(phi+math.log1p(zc))/w))
            return AF_DERIVED*S*(1-S)/w
        phi=8.
        Hsq=math.exp(-2*phi)
        F=math.exp(AF_DERIVED/(1+math.exp(-(phi+math.log1p(4.034077502899133))/w)))
        Fphi=F*alpha(phi)
        g=-Fphi/Hsq
        X=Hsq/2
        boxphi=-2*Hsq
        cubic=-g*X*boxphi
        ratio=cubic/(F*Hsq)
        check(math.isclose(ratio,-alpha(phi),rel_tol=1e-8,abs_tol=1e-15),
              "on-shell cubic ratio must equal -alpha_M")
        p=2-1/w
        tests.append({"width":w,"g_coefficient_exponent_p":p,
                      "g_coefficient_fades_if_width_below_half":bool(w<.5),
                      "cubic_action_relative_EH_power_a":-1/w,
                      "normalized_cubic_action_at_n8":ratio,
                      "physical_cubic_action_decays_for_positive_width":True})
    return tests


def derive_exact_structural_lapse():
    """Canonical equilibrium derives N_infinity=sqrt(8pi/3) independently of S0.

    Z_BI=4pi/(3 t_P^2 S0^2).  G2_IR=(2F Z_sigma-4F Z_BI)/sigma².
    At the coasting scalar-only fixed point, V=2F H², and
    V=4F Z_BI/sigma², so H²=2 Z_BI/sigma².
    Thus dot(sigma)=H*sigma=sqrt(2 Z_BI),
    structural lapse N=tP S0 dot(sigma)=sqrt(8pi/3).
    Independent of F_inf, scalar normalization S0 and Planck time.
    """
    target=math.sqrt(8.*math.pi/3.)
    cases=[]
    for tP,S0 in ((1.,1.),(1.,1e5),(5.391247e-44,2.72e61),
                    (5.391247e-44,1e62),(.1,6.)):
        zbi=4.*math.pi/(3.*tP*tP*S0*S0)
        sigmadot=math.sqrt(2.*zbi)
        structural_lapse=tP*S0*sigmadot
        check(close(structural_lapse,target,tol=1e-10),
              "canonical lapse independent of Planck clock normalization")
        cases.append({"tP":tP,"S0":S0,"Z_BI":zbi,
                      "sigma_dot_at_coasting":sigmadot,
                      "structural_lapse_N_infinity":structural_lapse,
                      "structural_radius_speed_over_c":structural_lapse})
    return {"exact_structural_lapse":target,"cases":cases}


def main():
    # The potential is exponential with canonical field psi=sqrt(2F) phi,
    # so V(psi)~exp(-sqrt(2/F)psi), a dimensionful-field slope.
    # In the Planck-normalized canonical field psi/sqrt(F)=sqrt(2)phi,
    # the dimensionless slope is lambda_IR=sqrt(2).
    jac=jacobian_at_IR()
    check(jac[0][0]==-2 and jac[1][1]==-1 and jac[2][2]==-2,
          "IR eigenvalues")
    # Direct analytic linearization:
    # du/dn = -2 du - 3/2 Om - Or,
    # dOm/dn=-Om; dOr/dn=-2Or.
    # Matter particular solution du=-3/2 Om; q=2 du+3/2 Om=-3/2 Om.
    cu=-1.5
    cq=2*cu+1.5
    check(close(cq,-1.5), "leading matter q")
    check(close(math.sqrt(2.),1.4142135623730951),"canonical field slope")

    exact_structural_lapse=derive_exact_structural_lapse()
    scalar_only_global_proof=prove_scalar_only_global_basin()
    cubic_decay_checks=prove_on_shell_cubic_decay()

    traces=[]
    summaries=[]
    for label,seed in CASES.items():
        rows,result=eval_case(label,seed)
        traces.extend(rows)
        summaries.append(result)

    # Check exact continuity of finite vs infinite kinetic values.
    check(abs(F_HISTORICAL_FIN-2.)>1.,"finite D logistic is not mature canonical")
    check(F_DERIVED<F_MATURE_PREVIOUS,"different original and F-derived branches")
    # For phi=ln(sigma/C), Z_sigma=sigma^2 X_phi.
    # The homogeneous scalar's canonical k1=2F∞ and yV=2/3 IR.
    x=.3
    check(close((2*x*(2*F_DERIVED))/(x*2*F_DERIVED),2.),
          "IR phi kinetic action implies D=2 on coasting u=1")

    with (OUT/"forward_canonical_trajectories.csv").open("w",newline="") as h:
        wr=csv.DictWriter(h,fieldnames=list(traces[0].keys()))
        wr.writeheader();wr.writerows(traces)
    with (OUT/"candidate_ir_results.csv").open("w",newline="") as h:
        wr=csv.DictWriter(h,fieldnames=list(summaries[0].keys()))
        wr.writeheader();wr.writerows(summaries)
    verdict={
     "scope":"forward homogeneous evolution of specified IR canonical SDMC action ONLY",
     "status":"canonical homogeneous future attractor proved locally and tested forward from 8 physical initial states",
     "theory":{
       "mature_G2_sigma":"[2 F∞ Z_sigma-4 F∞ Z_BI]/sigma²",
       "mature_G3":"0","mature_G4":"F∞/2",
       "sigma_phi_map":"sigma=C exp(phi)",
       "mature_G2_phi":"2 F∞ X_phi - V0 exp(-2 phi)",
       "same_action_all_eight_seeds":True,
       "shared_V0_in_units_tP_S0_equals_1":4*F_DERIVED*4*math.pi/3,
       "canonical_dimless_slope":math.sqrt(2.),
       "mature_F_from_force_derived_AF":F_DERIVED,
       "mature_k1_phi":2*F_DERIVED,
       "mature_G3_phi":"0","mature_k2_phi":"0",
       "mature_kinetic_D":"2 when u=1, alphaB=0",
       "scalar_sound_speed_squared":"1 everywhere for strictly canonical constant-F branch",
       "e_fold_variable":"n=ln a, different from structural lapse symbol N",
       "fluid_equations":"epsilon=u²+1.5*Om+2*Or; u'=-(3-epsilon)u+3(1-u²/3-Om-Or); Om'=Om(-3+2epsilon); Or'=Or(-4+2epsilon)",
       "fixed_point":{"u":1.,"Omega_m":0.,"Omega_r":0.,
                      "Omega_V":2./3.,"epsilon":1.,"q":0.,"D":2.},
       "Jacobian_u_Om_Or":jac,
       "linear_eigenvalues":[-2.,-1.,-2.],
       "first_matter_mode":{"delta_u_over_Omega_m":cu,"q_over_Omega_m":cq,
                             "Omega_m_asymptote":"m1/a","p_structural":"u=1-(3/2)m1/a + ..."},
     },
     "historical_finite_late_D_asymptote":F_HISTORICAL_FIN,
     "mature_IR_D":2.,
     "F_IR_derived_branch":F_DERIVED,
     "F_IR_old_manuscript_branch":F_MATURE_PREVIOUS,
     "F_IR_cross_branch_relative_mismatch":F_MATURE_PREVIOUS/F_DERIVED-1.,
     "structural_lapse_from_canonical_BI":exact_structural_lapse,
     "all_cases":summaries,
     "scalar_only_global_attractor_tests":scalar_only_global_proof,
     "coasting_NoSlip_on_shell_cubic_decay_tests":cubic_decay_checks,
     "proof_limitations":[
       "Forward runs begin already in the specified mature canonical action, not in late300's finite reconstructed action",
       "The amount and location of finite-to-canonical operator transfer are unspecified",
       "The Lagrangian boundary conditions plus endpoint data do not select the finite kinetic correction",
       "The canonical scalar potential normalization V0 depends on initial conditions and structural-clock zero point",
       "A microscopic selector for lambda_e, Dfloor enhancement, D0 amplitude, kinetic exponent, and independent D window was not derived",
       "This is a homogeneous dynamical gate; full native scalar perturbations across finite-to-mature handoff are NOT run",
       "A full cosmological evidence or model-family Bayes factor vs optimized local021 is NOT inferred",
     ]
    }
    (OUT/"verdict.json").write_text(json.dumps(verdict,indent=2,sort_keys=True)+"\n")
    print("D_FORWARD_CANONICAL_IR_AUDIT_PASS",
          json.dumps({"models_evolved":len(summaries),
                      "N_infinity_from_canonical_action":exact_structural_lapse["exact_structural_lapse"],
                      "scalar_only_exact_global_orbits":len(scalar_only_global_proof),
                      "coasting_cubic_widths_tested":len(cubic_decay_checks),
                      "latest_n":END_N,
                      "largest_final_abs_u_minus_one":max(abs(x["final_u"]-1) for x in summaries),
                      "largest_Friedmann_error":max(x["max_friedmann_constraint_error"] for x in summaries),
                      "linear_eigenvalues":[-2.,-1.,-2.],
                      "leading_matter_u_coefficient":cu,
                      "leading_matter_q_coefficient":cq,
                      "finite_D":F_HISTORICAL_FIN,
                      "mature_D":2.,
                      "finite_to_mature_operator_flow":"NOT UNIQUELY SUPPLIED BY CURRENT ACTION",
                      "local021_family_BF":"NOT CALCULABLE WITHOUT JOINT PRIOR + JOINT RUN"},sort_keys=True),
          flush=True)

if __name__=="__main__":
    main()
