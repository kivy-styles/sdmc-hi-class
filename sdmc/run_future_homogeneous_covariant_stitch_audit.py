#!/usr/bin/env python3
"""
Future homogeneous covariant stitch audit for the accepted SDMC action.

This script does NOT modify the successful z>=0 action or its likelihood
history.  It reconstructs the accepted structural-field action at sigma<=1,
attaches perturbation-ready future tails for sigma>1, and then evolves the
exact homogeneous Horndeski background equations forward with matter and
radiation dilution.  The F and g tails preserve the complete C3 jet at the
present boundary because hi_class uses third field derivatives of G3 and G4
in its perturbation gravity functions.

Action:
    G2 = k1(sigma) Z + k2(sigma) Z^2 - V(sigma)
    G3 = g(sigma) Z
    G4 = F(sigma)/2
    Z  = dot(sigma)^2/2

The future normal form is
    k1 -> kappa1 sigma^-2
    k2 -> kappa2 sigma^-2
    V  -> U0     sigma^-2
    F  -> F_inf
    g  -> 0

with the mature constant-Z coasting relations
    kappa1 + 2 kappa2 Z_inf = 2 F_inf
    U0 = 2 kappa1 Z_inf + 3 kappa2 Z_inf^2.

Two structural endpoints are audited:
  bounded_chi : chi_inf=1,          N_inf=2.860445...
  matched_chi : chi_inf=F_inf,      N_inf=2.894405...

The exact homogeneous equations are reduced to a 2x2 system for dot(v) and
dot(H), while H itself is obtained from the Friedmann constraint.  The script
reports both a term-balance No-Slip residual and the physical dimensionless
combination alpha_B+2 alpha_M.  The exact Bellini-Sawicki D and c_s^2 background health functions are checked
along the released future trajectory.  A future-capable hi_class perturbation
propagation remains a separate final implementation gate.
"""
from pathlib import Path
import json, math, re
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

BG=Path("output/linear_cov_free_00_background.dat")
OUT=Path("output/future_homogeneous_covariant_stitch_audit.json")

AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
MPC_TIME_GYR=MPC_M/C_MS/SEC_PER_GYR

XI=math.sqrt(8.*math.pi/3.)
N0_STRUCT=3.4135084646534866
CHI0=0.9410895692368776
FINF=math.exp(AF)

# Fixed reproducible candidate seeds obtained from the preceding design scan.
CANDIDATES=[
  dict(name="bounded_chi",chi_inf=1.0,
       r=3.436844444274898,muQ=3.7502361106872546),
  dict(name="matched_chi_F",chi_inf=FINF,
       r=3.756118836402892,muQ=3.753897593021392),
]

def read(path):
    lines=path.read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and
         re.search(r"(?:^|\s)1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

def coeffs_match(qder,qinf,mu):
    # Exponential-polynomial tail matching the complete supplied endpoint jet.
    # Three derivatives give C2 matching; four give perturbation-ready C3.
    A=[]
    for k in range(len(qder)):
        s=0.
        for j in range(k):
            s += math.comb(k,j)*((-mu)**(k-j))*A[j]
        A.append(qder[k]-(qinf if k==0 else 0.)-s)
    return A

def p_der(A,x,k):
    s=0.
    n=len(A)-1
    for j in range(k,n+1):
        s += A[j]*x**(j-k)/math.factorial(j-k)
    return s

def tail_der(A,qinf,mu,x,k):
    if k==0:
        return qinf+math.exp(-mu*x)*p_der(A,x,0)
    s=0.
    for j in range(k+1):
        s += math.comb(k,j)*((-mu)**(k-j))*p_der(A,x,j)
    return math.exp(-mu*x)*s

d=read(BG)
need=["z","proper time [Gyr]","H [1/Mpc]","(.)rho_smg","(.)p_smg",
      "(.)rho_b","(.)rho_cdm","(.)rho_g","(.)rho_ur"]
for k in need:
    if k not in d: raise RuntimeError(f"missing {k}")

z=np.asarray(d["z"]); N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d["(.)rho_smg"]); pre=np.asarray(d["(.)p_smg"])
tg=np.asarray(d["proper time [Gyr]"])
rb=np.asarray(d["(.)rho_b"]); rc=np.asarray(d["(.)rho_cdm"])
rg=np.asarray(d["(.)rho_g"]); rur=np.asarray(d["(.)rho_ur"])

o=np.argsort(N)
N,H,rho,pre,tg,z,rb,rc,rg,rur=[x[o] for x in
                                (N,H,rho,pre,tg,z,rb,rc,rg,rur)]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre,tg,z,rb,rc,rg,rur=[x[keep] for x in
                                (N,H,rho,pre,tg,z,rb,rc,rg,rur)]

lnH=CubicSpline(N,np.log(H))
h=lnH(N,1)
dotH=H*H*h
X=.5*H*H

Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH)
tt=np.tanh(xx); uu=1.-tt*tt
ST=.5*(1.+tt); S1=uu/(4.*WIDTH); S2=-tt*uu/(4.*WIDTH*WIDTH)
F=np.exp(AF*ST); F1=F*(AF*S1); F2=F*((AF*S1)**2+AF*S2)

g=-F1/(H*H)
g1=CubicSpline(N,g)(N,1)
Dtarget=DFLOOR+D0*ST**POWER
alphaM=F1/F; alphaB=-2.*alphaM
alphaK=Dtarget-1.5*alphaB*alphaB

R=3.*rho+2.*g1*X*X+6.*F1*H*H+3.*(F-1.)*H*H
P=3.*pre+2.*g1*X*X-2.*X*F2 \
  -3.*(F-1.)*H*H-2.*(F-1.)*dotH \
  -2.*F1*(H*H+dotH)
Cbg=(R+P)/(2.*X)
k2=(F*alphaK-Cbg+4.*g1*X-6.*g*H*H)/(4.*X)
k1=Cbg-2.*k2*X
V=.5*(R-P)-k2*X*X

t=tg*SEC_PER_GYR*C_MS/MPC_M
i0=int(np.argmin(np.abs(N)))
t0=float(t[i0])
sigma=t/t0
A=H*t
Apsi=A+h*A*A

# phi -> psi -> sigma field redefinition.
gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4+2.*g*A**2*Apsi
gs=gt/sigma**3
k1s=k1t/sigma**2
k2s=(k2t-2.*gt)/sigma**4
Vs=V

spg=CubicSpline(sigma,gs)
spk1=CubicSpline(sigma,k1s)
spk2=CubicSpline(sigma,k2s)
spV=CubicSpline(sigma,Vs)
spF=CubicSpline(sigma,F)

def qbar_derivs_from_sigma(sp):
    Q=float(sp(1.)); Q1=float(sp(1.,1)); Q2=float(sp(1.,2))
    return [Q,2.*Q+Q1,4.*Q+5.*Q1+Q2]

def xderivs_direct_from_sigma(sp):
    Q=float(sp(1.)); Q1=float(sp(1.,1)); Q2=float(sp(1.,2)); Q3=float(sp(1.,3))
    # x=ln(sigma): d/dx=sigma d/dsigma.
    return [Q,Q1,Q1+Q2,Q1+3.*Q2+Q3]

bars={
  "k1":qbar_derivs_from_sigma(spk1),
  "k2":qbar_derivs_from_sigma(spk2),
  "V":qbar_derivs_from_sigma(spV),
  "F":xderivs_direct_from_sigma(spF),
}
g_xder=[
  float(spg(1.)),
  float(spg(1.,1)),
  float(spg(1.,1)+spg(1.,2)),
  float(spg(1.,1)+3.*spg(1.,2)+spg(1.,3)),
]

v0=1./t0
Z0=.5*v0*v0
H0=float(H[i0])
q0_bg=-1.-float(h[i0])
p0_bg=1./(H0*t0)

# Action-unit matter/radiation densities: CLASS background output is
# multiplied by 8piG/3, whereas the manuscript action equations use 8piG.
rho_m0=3.*float(rb[i0]+rc[i0])
rho_r0=3.*float(rg[i0]+rur[i0])

def asym_Z_from_mu(mu):
    AFc=coeffs_match(bars["F"],FINF,mu)
    mug=mu+1.
    Agc=coeffs_match(g_xder,0.,mug)
    # With the C3-preserving cubic-polynomial tails, the leading
    # x^3 exp(-mu x) terms give the asymptotic No-Slip ratio.
    return mu*AFc[3]/(2.*Agc[3])

def build_candidate(spec):
    chiinf=spec["chi_inf"]
    Ninf=XI*math.sqrt(chiinf/FINF)
    Zinf=Z0*(Ninf/N0_STRUCT)**2
    r=spec["r"]; muQ=spec["muQ"]
    mu_rates={
      "k1":float(spec.get("mu_k1",muQ)),
      "k2":float(spec.get("mu_k2",muQ)),
      "V":float(spec.get("mu_V",muQ)),
    }

    kappa1=2.*FINF*(1.-r)
    kappa2=r*FINF/Zinf
    U0=FINF*Zinf*(4.-r)

    # The C3 asymptotic equation has multiple positive roots.  Use the
    # fastest branch near mu_F~5.6; the slower branches leave a longer
    # modified-gravity transient while reaching the same fixed point.
    rootfun=lambda mu: asym_Z_from_mu(mu)-Zinf
    muF=brentq(rootfun,5.0,6.2)
    mug=muF+1.

    AFc=coeffs_match(bars["F"],FINF,muF)
    Agc=coeffs_match(g_xder,0.,mug)
    asym={"k1":kappa1,"k2":kappa2,"V":U0}
    Aq={q:coeffs_match(bars[q],asym[q],mu_rates[q]) for q in asym}

    def action(sig):
        x=max(0.,math.log(sig))
        out={}
        for q in ["k1","k2","V"]:
            qb=[tail_der(Aq[q],asym[q],mu_rates[q],x,j) for j in range(3)]
            qx=[]
            for k in range(3):
                qx.append(math.exp(-2.*x)*sum(
                    math.comb(k,j)*((-2.)**(k-j))*qb[j]
                    for j in range(k+1)))
            out[q]=qx[0]
            out[q+"s"]=qx[1]/sig
            out[q+"ss"]=(qx[2]-qx[1])/(sig*sig)

        fx=[tail_der(AFc,FINF,muF,x,j) for j in range(3)]
        out["F"]=fx[0]; out["Fs"]=fx[1]/sig
        out["Fss"]=(fx[2]-fx[1])/(sig*sig)

        gx=[tail_der(Agc,0.,mug,x,j) for j in range(3)]
        out["g"]=gx[0]; out["gs"]=gx[1]/sig
        out["gss"]=(gx[2]-gx[1])/(sig*sig)
        return out

    return {
      **spec,"N_inf":Ninf,"Z_inf":Zinf,"muF":muF,"mug":mug,
      "mu_rates":mu_rates,
      "kappa1":kappa1,"kappa2":kappa2,"U0":U0,"action":action
    }

def rhs_diag(Ne,y,model):
    sig=float(y[0]); v=float(y[1])
    if sig<1.:
        raise RuntimeError("future stitch received sigma<1")
    q=model["action"](sig)
    k1v=q["k1"]; k2v=q["k2"]; Vv=q["V"]; Fv=q["F"]
    gv=q["g"]; gsv=q["gs"]; gssv=q["gss"]
    k1sv=q["k1s"]; k2sv=q["k2s"]; Vsv=q["Vs"]
    Fsv=q["Fs"]; Fssv=q["Fss"]

    rm=rho_m0*math.exp(-3.*Ne)
    rr=rho_r0*math.exp(-4.*Ne)
    Z=.5*v*v

    # Friedmann constraint: 3F H^2 - 3 A_H H - R_E = 0.
    RE=rm+rr+.5*k1v*v*v+(.75*k2v-.5*gsv)*v**4+Vv
    AH=gv*v**3-Fsv*v
    disc=9.*AH*AH+12.*Fv*RE
    if Fv<=0. or disc<=0.:
        raise RuntimeError(f"invalid Friedmann root F={Fv} disc={disc}")
    Hn=(3.*AH+math.sqrt(disc))/(6.*Fv)

    B=.25*k2v-gsv/6.
    Bs=.25*k2sv-gssv/6.
    Khom=k1v+12.*B*v*v+6.*Hn*gv*v
    Q=Fsv-gv*v*v

    R0=(rm+4.*rr/3.
        +k1v*v*v+(k2v-gsv)*v**4
        +3.*Hn*gv*v**3-Hn*Fsv*v+Fssv*v*v)

    C0=(Vsv+.5*k1sv*v*v+3.*Bs*v**4
        +3.*Hn*k1v*v+12.*Hn*B*v**3+2.*Hn*gsv*v**3
        +9.*Hn*Hn*gv*v*v-6.*Hn*Hn*Fsv)

    det=2.*Fv*Khom+3.*Q*Q
    if not np.isfinite(det) or det==0.:
        raise RuntimeError(f"singular homogeneous matrix det={det}")

    dotv=-(2.*Fv*C0+3.*Q*R0)/det
    dotH=(Q*C0-Khom*R0)/det

    ns=Z*gv+.5*Fsv
    nsden=abs(Z*gv)+.5*abs(Fsv)+1e-300
    # For G4=F/2 and G3=g(sigma)Z,
    # alpha_B+2 alpha_M = 2 v (Z g + F_,sigma/2)/(H F).
    noslip_alpha_combo=2.*v*ns/(Hn*Fv)

    # Exact Bellini-Sawicki alpha functions for
    # G2=k1 Z+k2 Z^2-V, G3=g Z, G4=F/2.
    alphaM=v*Fsv/(Hn*Fv)
    alphaB=v*(-Fsv+2.*Z*gv)/(Hn*Fv)
    alphaK=2.*Z*Khom/(Hn*Hn*Fv)
    Dfull=alphaK+1.5*alphaB*alphaB

    # Exact d alpha_B/d ln a, needed by the hi_class c_s^2 numerator.
    Cb=(-Fsv+2.*Z*gv)/Fv
    Cb_s=((-Fssv+2.*Z*gsv)*Fv-(-Fsv+2.*Z*gv)*Fsv)/(Fv*Fv)
    Cb_Z=2.*gv/Fv
    Abo=v/Hn
    Abo_N=Abo*(dotv/(Hn*v)-dotH/(Hn*Hn))
    Cb_N=(v/Hn)*Cb_s+(v*dotv/Hn)*Cb_Z
    alphaB_N=Abo_N*Cb+Abo*Cb_N

    # Keep two action-level density conventions separate.
    #
    # The positive SDMC structural source used in chi and p_X is the KGB
    # matter-side density BEFORE the non-minimal -3 H dot(F) term is folded
    # into the effective RHS:
    #
    #   rho_X^struct = k1 Z + 3 k2 Z^2 + V
    #                  + 6 H v Z g - 2 Z^2 g_,sigma
    #                = 3 H^2 (F + F') - rho_m-rho_r.
    #
    # This is Part-V Eq.(43)/(70) and Part-IX Eq.(458).  The metric-RHS
    # effective source in Eq.(455)/(457) is instead
    #
    #   rho_X^eff = rho_X^struct - 3 H dot(F),
    #
    # and can change sign during the handoff without invalidating the
    # positive structural source.
    rhoX_structural=(k1v*Z+3.*k2v*Z*Z+Vv
                     +6.*Hn*v*Z*gv-2.*Z*Z*gsv)
    rhoX_metric_rhs=rhoX_structural-3.*Hn*Fsv*v

    # Exact hi_class effective rho_smg and p_smg for this subclass.
    # These are the bookkeeping variables used by gravity_functions_smg.c
    # in the Bellini-Sawicki c_s^2 numerator. They are deliberately distinct
    # from the manuscript's positive action-level rho_X.
    rho_smg_class=((k1v*Z+3.*k2v*Z*Z+Vv-2.*gsv*Z*Z)/3.
                   -Hn*v*(Fsv-2.*Z*gv)-(Fv-1.)*Hn*Hn)
    p_smg_class=(k1v*Z+k2v*Z*Z-Vv-2.*gsv*Z*Z+2.*Fssv*Z
                 +3.*(Fv-1.)*Hn*Hn+2.*(Fv-1.)*dotH
                 +2.*Fsv*Hn*v+(Fsv-2.*Z*gv)*dotv)/3.

    csproxy=np.nan
    den=k1v+6.*k2v*Z
    if abs(den)>1e-300:
        csproxy=(k1v+2.*k2v*Z)/den

    return np.array([v/Hn,dotv/Hn]), {
      "H":Hn,
      "p":v/(Hn*sig),
      "q":-1.-dotH/(Hn*Hn),
      "Khom":Khom,
      "det":det,
      "N_struct":N0_STRUCT*v/v0,
      "N_ratio":v/v0,
      "Z_ratio":Z/Z0,
      "F":Fv,
      "noslip_rel":ns/nsden,
      "noslip_alpha_combo":noslip_alpha_combo,
      "alphaM":alphaM,
      "alphaB":alphaB,
      "alphaK":alphaK,
      "D_full":Dfull,
      "alphaB_N":alphaB_N,
      "rhoX_structural":rhoX_structural,
      "rhoX_metric_rhs":rhoX_metric_rhs,
      "rho_smg_class":rho_smg_class,
      "p_smg_class":p_smg_class,
      "rho_m_action":rm,
      "rho_r_action":rr,
      "cs2_kessence_proxy":csproxy,
    }

def evolve(model,Nmax=10.,npts=2001):
    def fun(Ne,y):
        return rhs_diag(Ne,y,model)[0]

    sol=solve_ivp(fun,(0.,Nmax),[1.,v0],rtol=3e-9,
                  atol=[1e-11,1e-14],max_step=.01,dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)

    NN=np.linspace(0.,Nmax,npts)
    yy=sol.sol(NN)
    diag=[
      rhs_diag(float(ne),[float(sig),float(v)],model)[1]
      for ne,(sig,v) in zip(NN,yy.T)
    ]
    HH=np.array([x["H"] for x in diag])

    # Structural activation history.  Use the positive structural source,
    # not the metric-RHS effective source containing -3 H dot(F).
    # Constants cancel in chi/chi0=(rho_X/rho_X0) sigma^2.
    rhoXa=np.array([x["rhoX_structural"] for x in diag])
    chi_hist=CHI0*(rhoXa/rhoXa[0])*yy[0]*yy[0]
    F_hist=np.array([x["F"] for x in diag])
    gamma_hist=chi_hist/F_hist
    mapper_hist=yy[0]/np.exp(NN)
    for j,x in enumerate(diag):
        x["chi_action"]=float(chi_hist[j])
        x["Gamma_action"]=float(gamma_hist[j])
        x["mapper_ratio_to_present"]=float(mapper_hist[j])

    # Exact hi_class/Bellini-Sawicki scalar sound speed in two equivalent
    # representations.  The first mirrors gravity_functions_smg.c; the second
    # is the compact Horndeski form.  Their agreement is an internal audit.
    bra=np.array([x["alphaB"] for x in diag])
    run=np.array([x["alphaM"] for x in diag])
    DD=np.array([x["D_full"] for x in diag])
    dbradN=np.array([x["alphaB_N"] for x in diag])

    cs2_source=np.empty_like(NN)
    cs2_compact=np.empty_like(NN)
    cs2num_source=np.empty_like(NN)
    cs2num_compact=np.empty_like(NN)

    for j,x in enumerate(diag):
        Fv=x["F"]; Hn=x["H"]; b=bra[j]; m=run[j]
        xp=x["rho_smg_class"]+x["p_smg_class"]
        mp=(x["rho_m_action"]+4.*x["rho_r_action"]/3.)/3.

        num1=((2.-b)*(b+2.*m)/2.
              +1.5*(2.-b)*xp/(Hn*Hn)
              -1.5*(2.-2.*Fv+b*Fv)*mp/(Hn*Hn*Fv)
              +dbradN[j])

        hln=-1.-x["q"]
        matter_action=x["rho_m_action"]+4.*x["rho_r_action"]/3.
        num2=((2.-b)*(-hln+.5*b+m)
              -matter_action/(Hn*Hn*Fv)
              +dbradN[j])

        cs2num_source[j]=num1
        cs2num_compact[j]=num2
        cs2_source[j]=num1/DD[j]
        cs2_compact[j]=num2/DD[j]

        x["cs2_hiclass_source"]=float(cs2_source[j])
        x["cs2_hiclass_compact"]=float(cs2_compact[j])
        x["cs2num_hiclass_source"]=float(num1)
        x["cs2num_hiclass_compact"]=float(num2)

    cs2_identity_err=float(np.max(np.abs(cs2_source-cs2_compact)))

    # hi_class Horndeski quasi-static effective Newton and slip diagnostics.
    # With alpha_T=alpha_H=0, beta_1=alpha_B+2 alpha_M.  On exact No-Slip
    # beta_1=0, so G_eff=1/F and slip=1.
    Geff=np.empty_like(NN)
    slip=np.empty_like(NN)
    for j,x in enumerate(diag):
        Fv=x["F"]; Hn=x["H"]; b=bra[j]; m=run[j]
        xp=x["rho_smg_class"]+x["p_smg_class"]
        mp=(x["rho_m_action"]+4.*x["rho_r_action"]/3.)/3.
        beta1=b+2.*m
        beta2=(2.*beta1
               -3.*(2.-2.*Fv+b*Fv)*mp/(Hn*Hn*Fv)
               -3.*(-2.+b)*xp/(Hn*Hn)
               +2.*dbradN[j])
        bb=b*beta1
        if abs(bb)<1e-30:
            ge=1./Fv
        else:
            ge=(1.-bb/(bb-beta2))/Fv
        snum=2.*m*beta1
        sden=snum+beta2
        sl=1. if abs(snum)<1e-30 or abs(sden)<1e-30 else 1.-snum/sden
        Geff[j]=ge
        slip[j]=sl
        x["G_eff_horndeski"]=float(ge)
        x["slip_horndeski"]=float(sl)
        x["G_eff_F_minus_1"]=float(ge*Fv-1.)

    dt=np.zeros_like(NN)
    for i in range(1,len(NN)):
        dN=NN[i]-NN[i-1]
        dt[i]=dt[i-1]+.5*dN*(1./HH[i]+1./HH[i-1])*MPC_TIME_GYR

    qq=np.array([x["q"] for x in diag])
    crossings=[]
    for i in range(len(NN)-1):
        if qq[i]==0. or qq[i]*qq[i+1]<0.:
            f=-qq[i]/(qq[i+1]-qq[i])
            ne=NN[i]+f*(NN[i+1]-NN[i])
            tt=dt[i]+f*(dt[i+1]-dt[i])
            crossings.append({
              "ln_a":float(ne),
              "a":float(math.exp(ne)),
              "delta_t_Gyr":float(tt)
            })

    imax=int(np.argmax(qq))
    ichi=int(np.argmax(chi_hist))
    ichimin=int(np.argmin(chi_hist))
    igamma=int(np.argmax(gamma_hist))
    igammamin=int(np.argmin(gamma_hist))
    nos=np.array([abs(x["noslip_rel"]) for x in diag])
    imns=int(np.argmax(nos))
    nosa=np.array([abs(x["noslip_alpha_combo"]) for x in diag])
    imnsa=int(np.argmax(nosa))
    csproxy=np.array([x["cs2_kessence_proxy"] for x in diag])
    Kh=np.array([x["Khom"] for x in diag])
    iD=int(np.argmin(DD))
    ics=int(np.argmin(cs2_compact))
    icsmax=int(np.argmax(cs2_compact))

    samples=[]
    for nt in [0.,.5,1.,2.,3.,5.,8.,10.]:
        j=int(np.argmin(abs(NN-nt)))
        row={
          k:float(v) for k,v in diag[j].items()
          if isinstance(v,(int,float,np.floating)) and np.isfinite(v)
        }
        row.update(
          ln_a=float(NN[j]),
          a=float(math.exp(NN[j])),
          sigma=float(yy[0,j]),
          delta_t_Gyr=float(dt[j]),
          D_horndeski=float(DD[j]),
          cs2_horndeski=float(cs2_compact[j])
        )
        samples.append(row)

    return {
      "present_eom_check":{
        "H_rel_vs_free":float(diag[0]["H"]/H0-1.),
        "p_minus_free":float(diag[0]["p"]-p0_bg),
        "q_minus_free":float(diag[0]["q"]-q0_bg),
        "alphaB_plus_2alphaM":float(diag[0]["alphaB"]+2.*diag[0]["alphaM"]),
        "D_horndeski":float(DD[0]),
        "cs2_horndeski":float(cs2_compact[0]),
        "cs2_source_compact_identity_max_abs":cs2_identity_err,
        "G_eff_horndeski":float(Geff[0]),
        "slip_horndeski":float(slip[0]),
      },
      "asymptotic_target":{
        "N_inf":model["N_inf"],
        "N_inf_over_N0":model["N_inf"]/N0_STRUCT,
        "Z_inf_over_Z0":model["Z_inf"]/Z0,
      },
      "tail_parameters":{
        "r":model["r"],"muQ":model["muQ"],
        "mu_rates":model.get("mu_rates",{"k1":model["muQ"],"k2":model["muQ"],"V":model["muQ"]}),
        "muF":model["muF"],
        "mug":model["mug"],"kappa1":model["kappa1"],
        "kappa2":model["kappa2"],"U0":model["U0"],
      },
      "health_diagnostics":{
        "min_Khom":float(np.min(Kh)),
        "min_sigma2_Khom":float(np.min(Kh*yy[0]*yy[0])),
        "min_D_horndeski":float(np.min(DD)),
        "min_D_horndeski_at_ln_a":float(NN[iD]),
        "max_D_horndeski":float(np.max(DD)),
        "min_cs2_horndeski":float(np.min(cs2_compact)),
        "min_cs2_horndeski_at_ln_a":float(NN[ics]),
        "max_cs2_horndeski":float(np.max(cs2_compact)),
        "max_cs2_horndeski_at_ln_a":float(NN[icsmax]),
        "cs2_source_compact_identity_max_abs":cs2_identity_err,
        "min_cs2_kessence_proxy":float(np.nanmin(csproxy)),
        "max_cs2_kessence_proxy":float(np.nanmax(csproxy)),
        "chi_action_max":float(chi_hist[ichi]),
        "chi_action_max_at_ln_a":float(NN[ichi]),
        "chi_action_min":float(chi_hist[ichimin]),
        "chi_action_final":float(chi_hist[-1]),
        "Gamma_action_max":float(gamma_hist[igamma]),
        "Gamma_action_max_at_ln_a":float(NN[igamma]),
        "Gamma_action_min":float(gamma_hist[igammamin]),
        "Gamma_action_final":float(gamma_hist[-1]),
        "mapper_ratio_final":float(mapper_hist[-1]),
        "max_abs_noslip_term_balance_rel":float(nos[imns]),
        "max_abs_noslip_term_balance_at_ln_a":float(NN[imns]),
        "max_abs_alphaB_plus_2alphaM":float(nosa[imnsa]),
        "max_abs_alphaB_plus_2alphaM_at_ln_a":float(NN[imnsa]),
        "max_abs_G_eff_F_minus_1":float(np.max(np.abs(Geff*np.array([x["F"] for x in diag])-1.))),
        "max_abs_slip_minus_1":float(np.max(np.abs(slip-1.))),
        "final_G_eff_horndeski":float(Geff[-1]),
        "final_slip_horndeski":float(slip[-1]),
      },
      "future_kinematics":{
        "q_max":float(qq[imax]),
        "q_max_ln_a":float(NN[imax]),
        "q_max_delta_t_Gyr":float(dt[imax]),
        "q_zero_crossings":crossings,
      },
      "final":{
        "ln_a":float(NN[-1]),
        "a":float(math.exp(NN[-1])),
        "sigma":float(yy[0,-1]),
        "delta_t_Gyr":float(dt[-1]),
        **{
          k:float(v) for k,v in diag[-1].items()
          if isinstance(v,(int,float,np.floating)) and np.isfinite(v)
        },
        "D_horndeski":float(DD[-1]),
        "cs2_horndeski":float(cs2_compact[-1]),
      },
      "samples":samples,
    }

# Accepted present alpha/c_s values, if present in the saved hi_class table.
accepted_present_health={}
for key in ["kineticity_smg","braiding_smg","M2_running_smg","c_s^2","kin (D)","M*^2_smg"]:
    if key in d:
        accepted_present_health[key]=float(np.asarray(d[key])[o][keep][i0])


def release_trajectory(model,Nmax=10.,npts=2001):
    sol=solve_ivp(lambda Ne,y:rhs_diag(Ne,y,model)[0],
                  (0.,Nmax),[1.,v0],rtol=3e-9,
                  atol=[1e-11,1e-14],max_step=.01,dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    NN=np.linspace(0.,Nmax,npts)
    yy=sol.sol(NN)
    diag=[rhs_diag(float(ne),[float(sig),float(v)],model)[1]
          for ne,(sig,v) in zip(NN,yy.T)]
    return NN,yy,diag

def refine_noslip_model(model,yy,blend_rate=500.):
    """
    Use the released Z(sigma) trajectory to reconstruct the exact on-trajectory
    No-Slip value g_NS=-F_,sigma/(2Z), while preserving the original future
    tail's complete C3 jet at sigma=1.

    The blend
      W=1-exp(-u)(1+u+u^2/2+u^3/6), u=blend_rate*ln(sigma)
    satisfies W^(n)(0)=0 for n=0,1,2,3 and tends to unity rapidly.  This
    matters because hi_class uses G3_phiphiphi and G4_phiphiphi in the
    perturbation gravity functions.
    """
    base=model["_base_action"]
    sig=np.asarray(yy[0]); vel=np.asarray(yy[1])
    x=np.log(sig); Ztraj=.5*vel*vel
    zsp=PchipInterpolator(x,Ztraj,extrapolate=False)

    xmax=max(12.,float(x[-1])+1.)
    xg=np.linspace(0.,xmax,10001)
    xclip=np.minimum(xg,x[-1])
    Zg=np.asarray(zsp(xclip))
    Zg[xg>x[-1]]=model["Z_inf"]

    gns=np.empty_like(xg)
    for i,(xx,zz) in enumerate(zip(xg,Zg)):
        aq=base(math.exp(xx))
        gns[i]=-aq["Fs"]/(2.*zz)
    gnsp=CubicSpline(xg,gns,bc_type="natural")

    def action(sigv):
        aq=base(sigv).copy()
        xx=max(0.,math.log(sigv))

        # Baseline g and its x=ln(sigma) derivatives.
        go=aq["g"]
        go1=sigv*aq["gs"]
        go2=sigv*sigv*aq["gss"]+sigv*aq["gs"]

        gn=float(gnsp(xx))
        gn1=float(gnsp(xx,1))
        gn2=float(gnsp(xx,2))
        de=gn-go; de1=gn1-go1; de2=gn2-go2

        u=blend_rate*xx
        ee=math.exp(-u)
        W=1.-ee*(1.+u+.5*u*u+u*u*u/6.)
        W1=blend_rate*ee*u*u*u/6.
        W2=blend_rate*blend_rate*ee*(3.*u*u-u*u*u)/6.

        gx=go+W*de
        gx1=go1+W1*de+W*de1
        gx2=go2+W2*de+2.*W1*de1+W*de2

        aq["g"]=gx
        aq["gs"]=gx1/sigv
        aq["gss"]=(gx2-gx1)/(sigv*sigv)
        return aq

    new=dict(model)
    new["action"]=action
    return new

def build_refined_noslip_candidate(spec,iterations=2,blend_rate=500.):
    model=build_candidate(spec)
    model["_base_action"]=model["action"]
    for _ in range(iterations):
        _,yy,_=release_trajectory(model)
        model=refine_noslip_model(model,yy,blend_rate=blend_rate)
    return model

if __name__ == "__main__":
    results={}
    refined_results={}
    for spec in CANDIDATES:
        model=build_candidate(spec)
        results[spec["name"]]=evolve(model)
        refined=build_refined_noslip_candidate(spec,iterations=2,blend_rate=500.)
        refined_results[spec["name"]]=evolve(refined)

    out={
      "status":(
        "action-native homogeneous future audit; z>=0 action is unchanged. "
        "The exact homogeneous Bellini-Sawicki D and c_s^2 functions are also "
        "evaluated along the future trajectory. Full perturbation propagation in "
        "a future-capable hi_class background remains outstanding. The physical "
        "No-Slip diagnostic is alpha_B+2 alpha_M; the separate term-balance "
        "ratio can look large when both terms are individually tiny."
      ),
      "accepted_present":{
        "H0_1Mpc":H0,"p0":p0_bg,"q0":q0_bg,
        "N0_struct":N0_STRUCT,"Z0":Z0,"F0":float(F[i0]),
        "F_inf":FINF,"Xi_v":XI,
        "hi_class_saved_health":accepted_present_health,
      },
      "candidates":results,
      "noslip_refined_candidates":refined_results,
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

    print("FUTURE_HOMOGENEOUS_COVARIANT_STITCH_AUDIT")
    for name,r in results.items():
        print("CANDIDATE",name,json.dumps({
          "present_check":r["present_eom_check"],
          "target":r["asymptotic_target"],
          "tail":r["tail_parameters"],
          "health":r["health_diagnostics"],
          "future":r["future_kinematics"],
          "final":r["final"],
        },sort_keys=True))

    for name,r in refined_results.items():
        print("NOSLIP_REFINED",name,json.dumps({
          "present_check":r["present_eom_check"],
          "target":r["asymptotic_target"],
          "tail":r["tail_parameters"],
          "health":r["health_diagnostics"],
          "future":r["future_kinematics"],
          "final":r["final"],
        },sort_keys=True))
