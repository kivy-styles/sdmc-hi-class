#!/usr/bin/env python3
from pathlib import Path
import re
import numpy as np
from scipy.interpolate import CubicSpline

TARGET=Path("output/linear_cov_target_00_background.dat")

def read(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

def replace_once(path,old,new):
    p=Path(path); s=p.read_text()
    if new in s: return
    if s.count(old)!=1:
        raise RuntimeError(f"{path}: anchor count={s.count(old)} for {old[:100]!r}")
    p.write_text(s.replace(old,new,1))

d=read(TARGET)
AF=0.023604633340554453
ZC=3.4328776987879466
WIDTH=0.36879721635160295
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.050275813910968366

z=np.asarray(d["z"]); N=-np.log1p(z)
H=np.asarray(d["H [1/Mpc]"])
rho=np.asarray(d["(.)rho_smg"]); pre=np.asarray(d["(.)p_smg"])
tg=np.asarray(d["proper time [Gyr]"])
o=np.argsort(N)
N,H,rho,pre,tg,z=[x[o] for x in (N,H,rho,pre,tg,z)]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre,tg,z=[x[keep] for x in (N,H,rho,pre,tg,z)]

lnH=CubicSpline(N,np.log(H)); h=lnH(N,1); Hpa=H*H*h
Xphi=.5*H*H
Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH); tt=np.tanh(xx); uu=1.-tt*tt
ST=.5*(1.+tt); S1=uu/(4.*WIDTH); S2=-tt*uu/(4.*WIDTH*WIDTH)
F=np.exp(AF*ST); F1=F*(AF*S1); F2=F*((AF*S1)**2+AF*S2)

g=-F1/(H*H)
g1=CubicSpline(N,g)(N,1)
Dtarget=DFLOOR+D0*ST**POWER
alphaM=F1/F; alphaB=-2.*alphaM
alphaK=Dtarget-1.5*alphaB*alphaB
R=3.*rho+2.*g1*Xphi*Xphi+6.*F1*H*H+3.*(F-1.)*H*H
P=3.*pre+2.*g1*Xphi*Xphi-2.*Xphi*F2 \
  -3.*(F-1.)*H*H-2.*(F-1.)*Hpa-2.*F1*(H*H+Hpa)
C=(R+P)/(2.*Xphi)
k2=(F*alphaK-C+4.*g1*Xphi-6.*g*H*H)/(4.*Xphi)
k1=C-2.*k2*Xphi
V=.5*(R-P)-k2*Xphi*Xphi

SEC_PER_GYR=1e9*365.25*86400.
C_MS=299792458.
MPC_M=3.085677581491367e22
t=tg*SEC_PER_GYR*C_MS/MPC_M
i0=int(np.argmin(np.abs(N)))
t0=float(t[i0])
psi=np.log(t/t0)

# phi=f(psi), A=dphi/dpsi=H t, A_psi=A+h A^2.
A=H*t
Apsi=A+h*A*A
Xpsi=.5/(t*t)

# Field-redefined Horndeski functions.
gt=g*A**3
k1t=k1*A**2
k2t=k2*A**4+2.*g*A**2*Apsi
Vt=V
Ft=F

# Exact algebraic closure checks before installing.
noslip=Xpsi*gt+0.5*(F1*A)
xmap=np.max(np.abs((A*A*Xpsi)/Xphi-1.))
if xmap>1e-10 or np.max(np.abs(noslip))>1e-10:
    raise RuntimeError(f"field-redefinition closure failed xmap={xmap} noslip={np.max(np.abs(noslip))}")

def carr(x):
    return ",".join(f"{v:.17e}" for v in np.asarray(x).ravel())

spl={}
for name,y in [("g",gt),("k1",k1t),("k2",k2t),("V",Vt),("F",Ft)]:
    spl[name]=CubicSpline(psi,y)
psiN=CubicSpline(N,psi)

hp=Path("gravity_smg/sdmc_late266_logstruct_table.h")
with hp.open("w") as f:
    f.write("#ifndef SDMC_LATE266_LOGSTRUCT_TABLE_H\n#define SDMC_LATE266_LOGSTRUCT_TABLE_H\n")
    f.write(f"#define SDMC_LS_N {len(psi)}\n")
    f.write(f"#define SDMC_LS_T0 {t0:.17e}\n")
    f.write("static const double sdmc_ls_x[SDMC_LS_N] = {"+carr(psi)+"};\n")
    f.write("static const double sdmc_ls_Nx[SDMC_LS_N] = {"+carr(N)+"};\n")
    for name,sp in spl.items():
        f.write(f"static const double sdmc_ls_{name}[4*(SDMC_LS_N-1)] = "+"{"+carr(np.asarray(sp.c.T).reshape(-1))+"};\n")
    f.write("static const double sdmc_ls_psiN[4*(SDMC_LS_N-1)] = {"+carr(np.asarray(psiN.c.T).reshape(-1))+"};\n")
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

replace_once("include/background.h",
             "sdmc_v3_covariant_linear_audit,",
             "sdmc_v3_covariant_linear_audit, sdmc_v3_covariant_logstruct_audit,")

replace_once("gravity_smg/gravity_models_smg.c",
             '#include "sdmc_late266_linear_covariant_table.h"\n',
             '#include "sdmc_late266_linear_covariant_table.h"\n#include "sdmc_late266_logstruct_table.h"\n')

parser_anchor='  if (strcmp(string1,"sdmc_v3_covariant_linear_audit") == 0) {'
parser=r'''  if (strcmp(string1,"sdmc_v3_covariant_logstruct_audit") == 0) {
    pba->gravity_model_smg = sdmc_v3_covariant_logstruct_audit;
    pba->field_evolution_smg = _TRUE_;
    pba->M2_evolution_smg = _FALSE_;
    flag2=_TRUE_;
    pba->parameters_size_smg = 1;
    class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);
    if (has_tuning_index_smg == _FALSE_) {
      pba->tuning_index_smg = 0;
      pba->tuning_dxdy_guess_smg = 1.;
    }
  }

  if (strcmp(string1,"sdmc_v3_covariant_linear_audit") == 0) {'''
replace_once("gravity_smg/gravity_models_smg.c",parser_anchor,parser)

ganchor='  else if (pba->gravity_model_smg == sdmc_v3_covariant_linear_audit) {'
gcode=r'''  else if (pba->gravity_model_smg == sdmc_v3_covariant_logstruct_audit) {
    double gt,gt1,gt2,gt3;
    double k1,k1p,k1pp,k1ppp;
    double k2,k2p,k2pp,k2ppp;
    double V,Vp,Vpp,Vppp;
    double F,Fp,Fpp,Fppp;
    sdmc_ls_eval(sdmc_ls_x,sdmc_ls_g,phi,&gt,&gt1,&gt2,&gt3);
    sdmc_ls_eval(sdmc_ls_x,sdmc_ls_k1,phi,&k1,&k1p,&k1pp,&k1ppp);
    sdmc_ls_eval(sdmc_ls_x,sdmc_ls_k2,phi,&k2,&k2p,&k2pp,&k2ppp);
    sdmc_ls_eval(sdmc_ls_x,sdmc_ls_V,phi,&V,&Vp,&Vpp,&Vppp);
    sdmc_ls_eval(sdmc_ls_x,sdmc_ls_F,phi,&F,&Fp,&Fpp,&Fppp);
    V += 3.*pba->H0*pba->H0*pba->parameters_smg[0];

    pgf->G2 = k1*X+k2*X*X-V;
    pgf->G2_X = k1+2.*k2*X;
    pgf->G2_XX = 2.*k2;
    pgf->G2_phi = k1p*X+k2p*X*X-Vp;
    pgf->G2_Xphi = k1p+2.*k2p*X;
    pgf->G2_XXphi = 2.*k2p;
    pgf->G2_phiphi = k1pp*X+k2pp*X*X-Vpp;
    pgf->G2_Xphiphi = k1pp+2.*k2pp*X;

    pgf->G3_X = gt;
    pgf->G3_XX = 0.;
    pgf->G3_XXX = 0.;
    pgf->G3_phi = gt1*X;
    pgf->G3_Xphi = gt1;
    pgf->G3_XXphi = 0.;
    pgf->G3_phiphi = gt2*X;
    pgf->G3_Xphiphi = gt2;
    pgf->G3_phiphiphi = gt3*X;

    pgf->DG4 = (F-1.)/2.;
    pgf->G4 = F/2.;
    pgf->G4_phi = Fp/2.;
    pgf->G4_phiphi = Fpp/2.;
    pgf->G4_phiphiphi = Fppp/2.;
  }

  else if (pba->gravity_model_smg == sdmc_v3_covariant_linear_audit) {'''
replace_once("gravity_smg/gravity_models_smg.c",ganchor,gcode)

p=Path("gravity_smg/gravity_models_smg.c")
s=p.read_text()
ic_anchor='case sdmc_v3_covariant_linear_audit:'
pos=s.index(ic_anchor,s.index("int gravity_models_initial_conditions_smg("))
ic=r'''case sdmc_v3_covariant_logstruct_audit:
    {
      double Na=log(a);
      double ps,ps1,ps2,ps3;
      sdmc_ls_eval(sdmc_ls_Nx,sdmc_ls_psiN,Na,&ps,&ps1,&ps2,&ps3);
      pvecback_integration[pba->index_bi_phi_smg]=ps;
      /* d psi / d tau = a / t, with t=t0 exp(psi). */
      pvecback_integration[pba->index_bi_phi_prime_smg]=a/(SDMC_LS_T0*exp(ps));
    }
    break;

  '''
s=s[:pos]+ic+s[pos:]
p.write_text(s)

out=[
 f"t0_Mpc={t0:.17e}",
 f"psi_min={psi.min():.17e}",
 f"psi_max={psi.max():.17e}",
 f"A_min={A.min():.17e}",
 f"A_max={A.max():.17e}",
 f"Xmap_max_rel={xmap:.17e}",
 f"noslip_redef_max_abs={np.max(np.abs(noslip)):.17e}",
 f"present_p={1./(H[i0]*t[i0]):.17e}",
 f"present_A={A[i0]:.17e}",
]
Path("output/logstruct_redefinition_summary.txt").write_text("\n".join(out)+"\n")
print("LOGSTRUCT_REDEFINITION")
print("\n".join(out))
