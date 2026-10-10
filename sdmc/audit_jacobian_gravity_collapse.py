#!/usr/bin/env python3
"""SDMC structural Jacobian/lapse -> expansion and matter-growth diagnostic.

Frozen late300 expansion, sourced from sdmc/apply_sdmc_full_background_patch.py.
Unlike a covariant evolution/Einstein-Boltzmann run, this is an inverse-background
identity and conditional pressureless, subhorizon growth test. Does not promote
new candidates, infer an independent scalar force, or modify late300.
"""
from pathlib import Path
import csv,json,math
from scipy.integrate import quad,solve_ivp

P=json.loads(Path("sdmc/late300_candidate.json").read_text())["parameters"]
C=299792.458
H0=P["H0"]; h=H0/100
Om=(P["omega_b"]+P["omega_cdm"])/h**2
Or=4.17998772e-5/h**2
Ox=1-Om-Or
dN=.5; tauA=.25; tauB=1.5
nT=-math.log1p(P["z_t"])
nC=-math.log1p(P["z_c"])
Xi0=3.41350846/1.03421566    # independent manuscript N0,p0 reference
N0=3.41350846
R0=Xi0*C/H0
zstar=1090.
Mstar=0.02946748423
xstar=-math.log1p(zstar)
Rstar=Mstar*R0/(1+zstar)
GYR_PER_INVERSE_KMSMPC=977.7922216807891
MPC_PER_GYR=C/(GYR_PER_INVERSE_KMSMPC) # WRONG? C*1Gyr /1Mpc (see below)
# dt=GYR_PER_INVERSE_KMSMPC/H0*(dx/E),
# c dt in Mpc = c/H0 *(dx/E).

def state(x):
    z=math.expm1(-x)
    r=Om*math.exp(-3*x); s=Or*math.exp(-4*x)
    b=r+s+Ox
    xr=s/(r+s)
    th=math.tanh((x-nT)/dN)
    W=.5*(1-th)
    Wx=-.5*(1-th*th)/dN
    f0=(3+xr)/P["lambda_e"]**2
    f=W*f0
    fx=Wx*f0-W*xr*(1-xr)/P["lambda_e"]**2
    d=P["A_late"]*z*math.exp(-z/tauA)-P["B_late"]*z*z*math.exp(-z/tauB)
    dz=P["A_late"]*math.exp(-z/tauA)*(1-z/tauA)-P["B_late"]*math.exp(-z/tauB)*(2*z-z*z/tauB)
    dxdelta=-(1+z)*dz
    E=math.sqrt(b)*(1+d)/math.sqrt(1-f)
    lnHx=.5*(-3*r-4*s)/b+dxdelta/(1+d)+.5*fx/(1-f)
    u=1/(1+math.exp(-(x-nC)/P["width"])) if x-nC> -200*P["width"] else 0.
    F=math.exp(P["A_F"]*u)
    return {"E":E,"q":-1-lnHx,"lnHx":lnHx,"Omega_m":r/(E*E),"F":F,
            "mu_no_slip_proxy":1/F,"Wtracker":W,"tracker_fraction":f}

def get_age_gyr(z):
    x=-math.log1p(z)
    return GYR_PER_INVERSE_KMSMPC/H0*quad(lambda t:1/state(t)["E"],-35,x,
        epsrel=2e-10,limit=500)[0]

# Reconstruct the structural radius with a fixed N: this MUST be positive.
It=quad(lambda t:1/state(t)["E"],xstar,0.,epsrel=2e-11)[0]
Rstar_if_Nconstant=R0-C*N0/H0*It
Nmean_needed=(R0-Rstar)*H0/(C*It)

def shape_one(x):
    a=math.exp(x)
    return 4*a*(1-a)
def shape_two(x):
    a=math.exp(x)
    return 6.75*a*(1-a)**2

lapse_families={}
for name, shape in [("symmetric_a_window",shape_one),("skewed_a_window",shape_two)]:
    Ig=quad(lambda x:shape(x)/state(x)["E"],xstar,0.,epsrel=2e-11)[0]
    amp=(N0*It-(R0-Rstar)*H0/C)/(N0*Ig)
    def Ns(x):return N0*(1-amp*shape(x))
    def Rs(x):return R0-C/H0*quad(lambda v:Ns(v)/state(v)["E"],x,0.,
                                  epsrel=2e-11,limit=250)[0]
    arr=[]
    for z in [0.,1.,3.,7.,10.,20.,50.,1090.]:
        x=-math.log1p(z)
        H=H0*state(x)["E"]
        R=Rs(x); N=Ns(x)
        idx=C*N/(H*R)
        M=R/(math.exp(x)*R0)
        arr.append({"z":z,"lapse_N":N,"R_Mpc":R,"M":M,"structural_exponent_p":idx,
                    "q_from_H":state(x)["q"],"H_km_s_Mpc":H})
    assert abs(arr[-1]["R_Mpc"]-Rstar)<2e-7
    assert abs(arr[0]["structural_exponent_p"]-1.03421566)<1e-6
    assert all(x["R_Mpc"]>0 for x in arr)
    lapse_families[name]={"amplitude_fitted_to_Mstar":amp,"R_values":arr}
assert Rstar_if_Nconstant<0

# Growth: δ''+(1-q)δ'-(3/2) Ωm μ δ=0. This is subhorizon pressureless,
# with same approximate initial growing-mode δ=δ'=1 at z_i=50.
# μ=1/F is a conditional No-Slip proxy, NOT the full Horndeski μ(k,z).
initial_x=-math.log1p(50.)
def growth(mode):
    def rhs(x,y):
        s=state(x)
        if mode=="lcdm":
            r=Om*math.exp(-3*x); ra=Or*math.exp(-4*x)
            b=r+ra+Ox
            q=-1+.5*(3*r+4*ra)/b; omega=r/b; mu=1.
        else:
            q=s["q"]; omega=s["Omega_m"]
            mu={"sdmc_unit_mu":1., "sdmc_noslip_proxy":1./s["F"],
                "sdmc_illustrative_mu1p1":1.1}[mode]
        return [y[1],-(1-q)*y[1]+1.5*omega*mu*y[0]]
    sol=solve_ivp(rhs,(initial_x,0),[1.,1.],dense_output=True,
                  rtol=2e-10,atol=5e-13,max_step=.08)
    if not sol.success:raise RuntimeError(sol.message)
    return {str(z):float(sol.sol(-math.log1p(z))[0]) for z in [50,30,20,16.19,10,7,3,1,0]}
growths={mode:growth(mode) for mode in
 ["lcdm","sdmc_unit_mu","sdmc_noslip_proxy","sdmc_illustrative_mu1p1"]}

# q can be recovered from N and p, but that is an identity once H was
# supplied; it does not independently derive the expansion history.
# q=p-1+(dlnp/dln a)-(dlnN/dln a).
# Verify via finite differences at z<=20 to avoid endpoint subtractive error.
jacobian_q_residual=[]
for name,f in [("symmetric_a_window",shape_one),("skewed_a_window",shape_two)]:
    amp=lapse_families[name]["amplitude_fitted_to_Mstar"]
    def Nm(x):return N0*(1-amp*f(x))
    def Rm(x):return R0-C/H0*quad(lambda t:Nm(t)/state(t)["E"],x,0.,
                                  epsrel=2e-11)[0]
    def px(x):return C*Nm(x)/(H0*state(x)["E"]*Rm(x))
    for z in [0.,3.,7.,20.]:
        xx=-math.log1p(z);hh=2e-4
        pval=px(xx)
        dlnp=(math.log(px(xx+hh))-math.log(px(xx-hh)))/(2*hh)
        dlnN=(math.log(Nm(xx+hh))-math.log(Nm(xx-hh)))/(2*hh)
        rec=pval-1+dlnp-dlnN
        actual=state(xx)["q"]
        jacobian_q_residual.append({"family":name,"z":z,"q_reconstructed":rec,
                                     "q_from_background":actual,"abs_residual":abs(rec-actual)})
assert max(a["abs_residual"] for a in jacobian_q_residual)<1e-5

snap=[]
for z in [0,.5,1,3,7,10,20,50,1090,3466]:
    s=state(-math.log1p(z))
    snap.append({"z":z,"H_km_s_Mpc":s["E"]*H0,"q":s["q"],
                 "Omega_m":s["Omega_m"],"F":s["F"],
                 "mu_no_slip_proxy":s["mu_no_slip_proxy"]})

diag={
"status":"JACOBIAN_RECONSTRUCTION_VALID_NONUNIQUE_LAPSE_GROWTH_CONDITIONAL",
"source":"late300 expansion_smg frozen; R0 inferred from N0,p0; recombination mapper Mstar from Manuscript B",
"scientific_disclaimer":"No independent scalar growth prediction, covariant nonlinear halo or JWST/SMBH likelihood. Lapse windows are inverse reconstructions fitted to two manuscript geometry endpoints.",
"baseline_late300_unchanged":True,
"geometry":{"R0_Mpc":R0,"Mstar":Mstar,"matched_Rstar_Mpc":Rstar,
           "recombination_manuscript_approx_Rstar_Mpc":.3848,
           "constant_N_Rstar_Mpc":Rstar_if_Nconstant,
           "constant_N_negative_unphysical":Rstar_if_Nconstant<0,
           "interval_average_N_required":Nmean_needed,"endpoint_N0":N0,
           "constant_N_insufficiency_Mpc":Rstar_if_Nconstant-Rstar,
           "constant_N_assumption":"diagnostic only; manuscript already permits evolving N across transitions",
           "reconstruction_same_late300_background_different_lapse":lapse_families},
"reconstruction_q":{"identity":"q=p-1+(d ln p/d ln a)-(d ln N/d ln a)",
                    "numerical_residuals":jacobian_q_residual,
                    "status":"tautological once H(a) supplied, not a new gravity force law"},
"background_snapshots":snap,
"conditional_growth":{"initial":"z=50; delta=1, ddelta/dln a=1; same normalisation all runs",
                      "equation":"delta_xx+(1-q)delta_x-1.5*Omega_m*mu*delta=0",
                      "mu_no_slip_proxy":"1/F(z) only; not full scale-dependent Horndeski response",
                      "growth_normalized_to_z50":growths,
                      "relative_sdmc_muF_to_sdmc_mu1":{
                           str(z):growths["sdmc_noslip_proxy"][str(z)]/growths["sdmc_unit_mu"][str(z)]
                           for z in [20,10,7,0]},
                      "relative_sdmc_muF_to_lcdm":{
                           str(z):growths["sdmc_noslip_proxy"][str(z)]/growths["lcdm"][str(z)]
                           for z in [20,10,7,0]},
                      "illustrative_mu1p1_not_candidate":True,
                      "scope":"pressureless subhorizon linear growth proxy; not radiation/acoustic physics, nonlinear collapse or JWST population"},
"conclusion":"Background Jacobian+N reproduces expansion q but does not determine N independently nor mu(k,z). No derived galaxy/SMBH accelerator. No candidate promoted."
}
out=Path("output/jacobian_gravity_collapse");out.mkdir(parents=True,exist_ok=True)
(out/"diagnostic.json").write_text(json.dumps(diag,indent=2,sort_keys=True)+"\n")
with (out/"background.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(snap[0]));writer.writeheader();writer.writerows(snap)
print("SDMC_JACOBIAN_GRAVITY_AUDIT",json.dumps({
"status":diag["status"],"R0_Mpc":R0,"Rstar_target_Mpc":Rstar,
"Rstar_constN_Mpc":Rstar_if_Nconstant,
"Nmean_required":Nmean_needed,
"q_z1090":state(-math.log1p(1090))["q"],"q_z20":state(-math.log1p(20))["q"],
"q_z7":state(-math.log1p(7))["q"],"q_z0":state(0.)["q"],
"growth_delta_ratio_mu1overF_to_unit_z7":diag["conditional_growth"]["relative_sdmc_muF_to_sdmc_mu1"]["7"],
"max_q_residual":max(a["abs_residual"] for a in jacobian_q_residual),
"promote_candidate":False},sort_keys=True))
