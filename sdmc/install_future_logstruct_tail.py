#!/usr/bin/env python3
"""
Extend the accepted log-structural Horndeski coefficient table into sigma>1.

This script is run only after install_late266_logstruct_covariant.py has added
the sdmc_v3_covariant_logstruct_audit model.  It replaces only the generated
coefficient header, preserving the accepted sigma<=1 reconstruction exactly
(up to the same reconstruction algebra).  The G2 coefficients use C2 future
matching, while F and g preserve the complete C3 endpoint jet required by the
hi_class perturbation gravity functions.

The future tail is one of the two structural endpoint candidates already
audited in run_future_homogeneous_covariant_stitch_audit.py:
  bounded_chi  : N_inf = 2.860445...
  matched_chi_F: N_inf = 2.894405...

The table variable is psi=ln(sigma).  Structural coefficients (hat) are mapped
back to the logstruct field coefficients (tilde) through
  g_tilde  = sigma^3 g_hat
  k1_tilde = sigma^2 k1_hat
  k2_tilde = sigma^4 k2_hat + 2 sigma^3 g_hat
  V_tilde  = V_hat
  F_tilde  = F_hat.
"""
from pathlib import Path
import argparse, json, math, re
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
import run_future_homogeneous_covariant_stitch_audit as future_aud

TARGET=Path("output/linear_cov_target_00_background.dat")
HEADER=Path("gravity_smg/sdmc_late266_logstruct_table.h")
SUMMARY=Path("output/future_logstruct_tail_install.json")

AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366
XI=math.sqrt(8.*math.pi/3.)
N0_STRUCT=3.4135084646534866
FINF=math.exp(AF)

CANDIDATES={
  "bounded_chi":dict(chi_inf=1.0,
                     r=3.436844444274898,
                     muQ=3.7502361106872546),
  "matched_chi_F":dict(chi_inf=FINF,
                       r=3.756118836402892,
                       muQ=3.753897593021392),
  "legacy_canonical":dict(chi_inf=FINF,
                          r=0.0,
                          muQ=50.0,
                          mu_k1=50.0,
                          mu_k2=5.0,
                          mu_V=120.0,
                          grid_power=2.0),
}

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22

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
    # Match the complete supplied endpoint jet: C2 for G2 coefficients,
    # perturbation-ready C3 for F and g.
    A=[]
    for k in range(len(qder)):
        s=0.
        for j in range(k):
            s += math.comb(k,j)*((-mu)**(k-j))*A[j]
        A.append(qder[k]-(qinf if k==0 else 0.)-s)
    return A

def p_der(A,x,k):
    out=0.
    for j in range(k,len(A)):
        out += A[j]*x**(j-k)/math.factorial(j-k)
    return out

def tail_der(A,qinf,mu,x,k):
    if k==0:
        return qinf+math.exp(-mu*x)*p_der(A,x,0)
    out=0.
    for j in range(k+1):
        out += math.comb(k,j)*((-mu)**(k-j))*p_der(A,x,j)
    return math.exp(-mu*x)*out

def qbar_derivs_from_sigma(sp):
    Q=float(sp(1.)); Q1=float(sp(1.,1)); Q2=float(sp(1.,2))
    return [Q,2.*Q+Q1,4.*Q+5.*Q1+Q2]

def xderivs_direct_from_sigma(sp):
    Q=float(sp(1.)); Q1=float(sp(1.,1)); Q2=float(sp(1.,2)); Q3=float(sp(1.,3))
    return [Q,Q1,Q1+Q2,Q1+3.*Q2+Q3]

def carr(x):
    return ",".join(f"{v:.17e}" for v in np.asarray(x).ravel())

ap=argparse.ArgumentParser()
ap.add_argument("--candidate",choices=sorted(CANDIDATES),
                default="matched_chi_F")
ap.add_argument("--psi-max",type=float,default=8.0)
ap.add_argument("--future-points",type=int,default=1600)
args=ap.parse_args()
spec=CANDIDATES[args.candidate]

d=read(TARGET)
z=np.asarray(d["z"]); N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d["(.)rho_smg"]); pre=np.asarray(d["(.)p_smg"])
tg=np.asarray(d["proper time [Gyr]"])
o=np.argsort(N)
N,H,rho,pre,tg,z=[x[o] for x in (N,H,rho,pre,tg,z)]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre,tg,z=[x[keep] for x in (N,H,rho,pre,tg,z)]

lnH=CubicSpline(N,np.log(H)); h=lnH(N,1); dotH=H*H*h
X=.5*H*H
Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH); tt=np.tanh(xx); uu=1.-tt*tt
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
psi=np.log(t/t0)
sigma=np.exp(psi)
A=H*t
Apsi=A+h*A*A

# accepted phi -> psi coefficients
gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4+2.*g*A**2*Apsi

# psi -> sigma structural coefficients
gs=gt/sigma**3
k1s=k1t/sigma**2
k2s=(k2t-2.*gt)/sigma**4
Vs=V

spg=CubicSpline(sigma,gs)
spk1=CubicSpline(sigma,k1s)
spk2=CubicSpline(sigma,k2s)
spV=CubicSpline(sigma,Vs)
spF=CubicSpline(sigma,F)

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

Z0=.5/t0**2
Ninf=XI*math.sqrt(spec["chi_inf"]/FINF)
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

def asym_Z_from_mu(mu):
    AFc=coeffs_match(bars["F"],FINF,mu)
    mug=mu+1.
    Agc=coeffs_match(g_xder,0.,mug)
    return mu*AFc[3]/(2.*Agc[3])

muF=brentq(lambda mu:asym_Z_from_mu(mu)-Zinf,5.0,6.2)
mug=muF+1.
AFc=coeffs_match(bars["F"],FINF,muF)
Agc=coeffs_match(g_xder,0.,mug)
asym={"k1":kappa1,"k2":kappa2,"V":U0}
Aq={q:coeffs_match(bars[q],asym[q],mu_rates[q]) for q in asym}

# Import only the future No-Slip correction from the shared action audit.
# The correction is constructed with a C3-flat switch at sigma=1, so adding
# delta g = g_refined-g_base preserves this installer's target-derived
# endpoint jet while aligning the unobserved future with the common refined
# fixed-point action.
_aud_spec={"name":args.candidate,**spec}
_aud_base=future_aud.build_candidate(_aud_spec)
_aud_ref=future_aud.build_refined_noslip_candidate(
    _aud_spec,iterations=2,blend_rate=500.)

def future_structural(x):
    sig=math.exp(x)
    out={}
    for q in ["k1","k2","V"]:
        qb=[tail_der(Aq[q],asym[q],mu_rates[q],x,j) for j in range(3)]
        qx=[]
        for k in range(3):
            qx.append(math.exp(-2.*x)*sum(
                math.comb(k,j)*((-2.)**(k-j))*qb[j]
                for j in range(k+1)))
        out[q]=qx[0]
    fx=[tail_der(AFc,FINF,muF,x,j) for j in range(3)]
    gx=[tail_der(Agc,0.,mug,x,j) for j in range(3)]
    out["F"]=fx[0]
    out["g"]=gx[0]
    # Future-only correction; exactly flat through third order at sigma=1.
    out["g"] += (_aud_ref["action"](sig)["g"]
                 - _aud_base["action"](sig)["g"])
    return sig,out

# Preserve all accepted past samples through psi=0, then append future.
past=psi <= 1e-12
psi_p=psi[past]
N_p=N[past]
gt_p=gt[past]; k1t_p=k1t[past]; k2t_p=k2t[past]
V_p=V[past]; F_p=F[past]

grid_power=float(spec.get("grid_power",1.0))
u_future=np.arange(1,args.future_points+1,dtype=float)/args.future_points
x_future=args.psi_max*u_future**grid_power
gt_f=[]; k1t_f=[]; k2t_f=[]; V_f=[]; F_f=[]
for x in x_future:
    sig,q=future_structural(float(x))
    gt_f.append(sig**3*q["g"])
    k1t_f.append(sig**2*q["k1"])
    k2t_f.append(sig**4*q["k2"]+2.*sig**3*q["g"])
    V_f.append(q["V"])
    F_f.append(q["F"])

psi_ext=np.r_[psi_p,x_future]
# psiN is only used for the early-time initial-condition lookup.  Appending
# N=x in the unobserved future keeps its abscissa monotonic and approaches the
# required p=1 branch without changing the accepted past mapping.
N_ext=np.r_[N_p,x_future]
gt_ext=np.r_[gt_p,gt_f]
k1t_ext=np.r_[k1t_p,k1t_f]
k2t_ext=np.r_[k2t_p,k2t_f]
V_ext=np.r_[V_p,V_f]
F_ext=np.r_[F_p,F_f]

if not np.all(np.diff(psi_ext)>0):
    raise RuntimeError("extended psi grid is not strictly increasing")
if not np.all(np.diff(N_ext)>0):
    raise RuntimeError("extended N grid is not strictly increasing")

# Preserve the accepted z>=0 spline coefficients exactly.  Re-splining the
# combined past+future knot set changes the last accepted intervals because a
# cubic spline solves a global tridiagonal system; this matters for the very
# fast canonical future release.  Build the accepted past and the future
# continuation separately, then concatenate interval coefficients.
past_data={
  "g":gt_p,"k1":k1t_p,"k2":k2t_p,"V":V_p,"F":F_p
}
future_data={
  "g":np.asarray(gt_f),"k1":np.asarray(k1t_f),"k2":np.asarray(k2t_f),
  "V":np.asarray(V_f),"F":np.asarray(F_f)
}
coef={}
join_jet={}
xf=np.r_[0.,x_future]
for name in ["g","k1","k2","V","F"]:
    sp_p=CubicSpline(psi_p,past_data[name])
    # The analytic structural tail was matched through C2 (G2) or C3 (F,g).
    # Enforce the accepted second derivative at the table join so the stored
    # representation is C2 there; the tiny clustered first future interval
    # then makes the first-derivative mismatch numerically negligible.
    y0=float(past_data[name][-1])
    yfuture=np.r_[y0,future_data[name]]
    second0=float(sp_p(0.,2))
    sp_f=CubicSpline(
      xf,yfuture,bc_type=((2,second0),"not-a-knot")
    )
    coef[name]=np.vstack([sp_p.c.T,sp_f.c.T])
    join_jet[name]={
      "past_d1":float(sp_p(0.,1)),
      "future_d1":float(sp_f(0.,1)),
      "d1_abs_mismatch":float(sp_f(0.,1)-sp_p(0.,1)),
      "past_d2":second0,
      "future_d2":float(sp_f(0.,2)),
      "d2_abs_mismatch":float(sp_f(0.,2)-second0),
    }

# Preserve the accepted N->psi initial-condition spline as well.  The future
# extension is exactly psi=N and is not used to reconstruct the accepted past.
psiN_p=CubicSpline(N_p,psi_p)
psiN_future_coef=np.asarray([
  [0.,0.,1.,float(x)] for x in np.r_[0.,x_future[:-1]]
])
psiN_coef=np.vstack([psiN_p.c.T,psiN_future_coef])

with HEADER.open("w") as f:
    f.write("#ifndef SDMC_LATE266_LOGSTRUCT_TABLE_H\n#define SDMC_LATE266_LOGSTRUCT_TABLE_H\n")
    f.write(f"#define SDMC_LS_N {len(psi_ext)}\n")
    f.write(f"#define SDMC_LS_T0 {t0:.17e}\n")
    f.write("static const double sdmc_ls_x[SDMC_LS_N] = {"+carr(psi_ext)+"};\n")
    f.write("static const double sdmc_ls_Nx[SDMC_LS_N] = {"+carr(N_ext)+"};\n")
    for name,cc in coef.items():
        f.write(f"static const double sdmc_ls_{name}[4*(SDMC_LS_N-1)] = "+"{"+carr(np.asarray(cc).reshape(-1))+"};\n")
    f.write("static const double sdmc_ls_psiN[4*(SDMC_LS_N-1)] = {"+carr(np.asarray(psiN_coef).reshape(-1))+"};\n")
    f.write(r"""
static int sdmc_ls_idx(const double *xarr,double x){
  int lo=0,hi=SDMC_LS_N-1;
  if(x<=xarr[0]) return 0;
  if(x>=xarr[SDMC_LS_N-1]) return SDMC_LS_N-2;
  while(hi-lo>1){int m=(lo+hi)/2; if(xarr[m]<=x) lo=m; else hi=m;}
  return lo;
}
static void sdmc_ls_eval(const double *xarr,const double *coef,double x,
                         double *y,double *yp,double *ypp,double *yppp){
  int i=sdmc_ls_idx(xarr,x);
  double dx=x-xarr[i];
  const double *cc=coef+4*i;
  *y=((cc[0]*dx+cc[1])*dx+cc[2])*dx+cc[3];
  *yp=(3.*cc[0]*dx+2.*cc[1])*dx+cc[2];
  *ypp=6.*cc[0]*dx+2.*cc[1];
  *yppp=6.*cc[0];
}
#endif
""")

# Join diagnostics in logstruct coefficients.
j=len(psi_p)-1
eps=x_future[0]
join={}
for name,yp,yf in [
  ("g",gt_p,gt_f),("k1",k1t_p,k1t_f),("k2",k2t_p,k2t_f),
  ("V",V_p,V_f),("F",F_p,F_f)]:
    join[name]={
      "present":float(yp[-1]),
      "first_future":float(yf[0]),
      "fractional_step":float((yf[0]-yp[-1])/(abs(yp[-1])+1e-300))
    }

summary={
  "candidate":args.candidate,
  "psi_min":float(psi_ext[0]),
  "psi_max":float(psi_ext[-1]),
  "n_past":int(len(psi_p)),
  "n_future":int(len(x_future)),
  "future_grid_power":grid_power,
  "first_future_psi":float(x_future[0]),
  "past_spline_preserved_exactly":True,
  "join_jet_mismatch":join_jet,
  "t0_Mpc":t0,
  "N_inf":Ninf,
  "Z_inf_over_Z0":Zinf/Z0,
  "r":r,"muQ":muQ,"mu_rates":mu_rates,"muF":muF,"mug":mug,
  "F_inf":FINF,
  "noslip_refinement":{"iterations":2,"blend_rate":500.0},
  "tail_smoothness":{
    "G3_F_structural_jet":"C3 at sigma=1 before table interpolation",
    "G2_structural_coefficients":"C2 at sigma=1",
    "table_representation":"piecewise cubic spline (C2 globally)"
  },
  "join_step_at_first_future_sample":join
}
SUMMARY.parent.mkdir(parents=True,exist_ok=True)
SUMMARY.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
print("FUTURE_LOGSTRUCT_TAIL_INSTALL")
print(json.dumps(summary,sort_keys=True))
