#!/usr/bin/env python3
from pathlib import Path
import numpy as np, re
from scipy.interpolate import CubicSpline

AF=0.020284758737310768
ZC=3.6509963892400266
WIDTH=0.35059972247108817
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.05064729745686054
H0=69.71482083084993/299792.458
XSTAR=H0*H0

def replace_once(path, old, new):
    p=Path(path); s=p.read_text()
    if new in s:
        return
    if s.count(old)!=1:
        raise RuntimeError(f"{path}: anchor count {s.count(old)} for {old[:80]!r}")
    p.write_text(s.replace(old,new,1))

path=Path("output/cov_exact_target_00_background.dat")
lines=path.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr))
names=[]
for i,m in enumerate(marks):
    end=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():end].strip())
a=np.loadtxt(path)
d={n:a[:,i] for i,n in enumerate(names)}

order=np.argsort(-np.log1p(d["z"]))
N=-np.log1p(d["z"][order])
H=d["H [1/Mpc]"][order]
rho=d["rho_tot_wo_smg_dbg"][order]
pres=d["p_tot_wo_smg_dbg"][order]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pres=N[keep],H[keep],rho[keep],pres[keep]

lnH=CubicSpline(N,np.log(H))
h=lnH(N,1)

Nc=-np.log1p(ZC)
x=(N-Nc)/(2*WIDTH)
t=np.tanh(x); u=1-t*t
S=.5*(1+t)
Sp=u/(4*WIDTH)
Spp=-t*u/(4*WIDTH*WIDTH)
F=np.exp(AF*S)
Fp=F*AF*Sp
Fpp=F*((AF*Sp)**2+AF*Spp)
g=Fp/2.
gp=Fpp/2.

D=DFLOOR+D0*S**POWER
alphaM=Fp/F
alphaB=-2*alphaM
alphaK=D-1.5*alphaB*alphaB

X=.5*H*H
lnx=np.log(X/XSTAR)

# Exact-logarithmic No-Slip ansatz:
# G4=F/2, G3=-G4_phi ln(X/Xstar), G2=K X + L X^2 - V.
# These reconstruction identities are evaluated on phi=N, X=H^2/2.
R=-4*g*(h-2.)-2*F*h-3*(rho+pres)/(H*H)
Q=F*alphaK
Y=(Q-R)/4.             # Y = L*X on target trajectory
K=1.5*R-.5*Q-2*gp*(1.+lnx)
L=Y/X
V=3*F*H*H+12*g*H*H-3*rho-K*X-3*L*X*X-2*X*gp*lnx
U=V/(H*H)
q=np.log(X)

for arr,name in [(K,"K"),(Y,"Y"),(U,"U"),(q,"q")]:
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"nonfinite {name}")

# Preserve the complete native target grid.  The earlier 1200-node
# uniform compression introduced percent-level errors in second derivatives
# near the low-z end even though the underlying all-node reconstruction closed
# the covariant equations to ~1e-5 or better.  Since the free field evolution
# is sensitive to those derivatives, use the original monotonic N nodes
# directly rather than resampling the reconstructed action.
nk=len(N)
grid=N.copy()
spl={}
for name,val in [("K",K),("Y",Y),("U",U),("QX",q)]:
    spl[name]=CubicSpline(grid,val)

def carr(x):
    return ",".join(f"{v:.17e}" for v in np.asarray(x).ravel())

head=Path("gravity_smg/sdmc_covariant_exact_data.h")
txt=[
"#ifndef __SDMC_COVARIANT_EXACT_DATA__",
"#define __SDMC_COVARIANT_EXACT_DATA__",
f"#define SDMC_COV_N {nk}",
f"#define SDMC_COV_AF {AF:.17e}",
f"#define SDMC_COV_ZC {ZC:.17e}",
f"#define SDMC_COV_WIDTH {WIDTH:.17e}",
f"#define SDMC_COV_XSTAR {XSTAR:.17e}",
"static const double sdmc_cov_x[SDMC_COV_N] = {"+carr(grid)+"};"
]
for name,s in spl.items():
    coeff=np.asarray(s.c.T).reshape(-1)
    txt.append(f"static const double sdmc_cov_{name}[4*(SDMC_COV_N-1)] = "+"{"+carr(coeff)+"};")
txt.append(r'''
static void sdmc_cov_eval(const double *coef,double x,double *f,double *fp,double *fpp){
  int lo=0,hi=SDMC_COV_N-1;
  if(x<=sdmc_cov_x[0]) lo=0;
  else if(x>=sdmc_cov_x[SDMC_COV_N-1]) lo=SDMC_COV_N-2;
  else {
    while(hi-lo>1){int mid=(lo+hi)/2; if(sdmc_cov_x[mid]<=x) lo=mid; else hi=mid;}
  }
  double dx=x-sdmc_cov_x[lo];
  const double *c=coef+4*lo;
  *f=((c[0]*dx+c[1])*dx+c[2])*dx+c[3];
  *fp=(3*c[0]*dx+2*c[1])*dx+c[2];
  *fpp=6*c[0]*dx+2*c[1];
}
#endif
''')
head.write_text("\n".join(txt))

Path("output/cov_exact_reconstruction_summary.txt").write_text(
    f"N=[{N[0]:.9g},{N[-1]:.9g}]\n"
    f"K=[{K.min():.9g},{K.max():.9g}]\n"
    f"LX=[{Y.min():.9g},{Y.max():.9g}]\n"
    f"V_H2=[{U.min():.9g},{U.max():.9g}]\n"
    f"D_independent={DFLOOR:.17g}+{D0:.17g}*S\\n"
    f"alphaK_z0={alphaK[-1]:.12g}\n"
)

replace_once(
    "include/background.h",
    "sdmc_v3_native_noslip,",
    "sdmc_v3_native_noslip, sdmc_v3_covariant_exact,"
)

replace_once(
    "gravity_smg/gravity_models_smg.c",
    '#include "gravity_models_smg.h"\n',
    '#include "gravity_models_smg.h"\n#include "sdmc_covariant_exact_data.h"\n'
)

parser_anchor='  if (strcmp(string1,"eft_alphas_power_law") == 0) {'
parser=r'''  if (strcmp(string1,"sdmc_v3_covariant_exact") == 0) {
    pba->gravity_model_smg = sdmc_v3_covariant_exact;
    pba->field_evolution_smg = _TRUE_;
    pba->M2_evolution_smg = _FALSE_;
    flag2=_TRUE_;
    pba->parameters_size_smg = 1;
    class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);
    if (has_tuning_index_smg == _FALSE_) {
      pba->tuning_index_smg = 0;
    }
    if (has_dxdy_guess_smg == _FALSE_) {
      pba->tuning_dxdy_guess_smg = 1.;
    }
  }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
replace_once("gravity_smg/gravity_models_smg.c",parser_anchor,parser)

ganchor='  else if(pba->gravity_model_smg == nkgb){'
gcode=r'''  else if(pba->gravity_model_smg == sdmc_v3_covariant_exact){
    double K,Kp,Kpp,Y,Yp,Ypp,U,Up,Upp,q,qp,qpp;
    sdmc_cov_eval(sdmc_cov_K,phi,&K,&Kp,&Kpp);
    sdmc_cov_eval(sdmc_cov_Y,phi,&Y,&Yp,&Ypp);
    sdmc_cov_eval(sdmc_cov_U,phi,&U,&Up,&Upp);
    sdmc_cov_eval(sdmc_cov_QX,phi,&q,&qp,&qpp);

    double Xt=exp(q);
    double L=Y/Xt;
    double Lp=(Yp-Y*qp)/Xt;
    double Lpp=(Ypp-2.*Yp*qp+Y*(qp*qp-qpp))/Xt;

    /* Constant potential offset is the one shooting parameter. */
    double V=2.*Xt*U + 3.*pba->H0*pba->H0*pba->parameters_smg[0];
    double Vp=2.*Xt*(Up+U*qp);
    double Vpp=2.*Xt*(Upp+2.*Up*qp+U*(qpp+qp*qp));

    double Nc=-log(1.+SDMC_COV_ZC);
    double xx=(phi-Nc)/(2.*SDMC_COV_WIDTH);
    double tt=tanh(xx);
    double uu=1.-tt*tt;
    double w=SDMC_COV_WIDTH;
    double S=.5*(1.+tt);
    double S1=uu/(4.*w);
    double S2=-tt*uu/(4.*w*w);
    double S3=-uu*(1.-3.*tt*tt)/(8.*w*w*w);
    double S4=tt*uu*(2.-3.*tt*tt)/(4.*w*w*w*w);
    double A=SDMC_COV_AF;
    double F=exp(A*S);
    double F1=F*(A*S1);
    double F2=F*(A*A*S1*S1+A*S2);
    double F3=F*(A*A*A*S1*S1*S1+3.*A*A*S1*S2+A*S3);
    double F4=F*(pow(A*S1,4)+6.*A*A*A*S1*S1*S2
                 +3.*A*A*S2*S2+4.*A*A*S1*S3+A*S4);

    double g=.5*F1, gp=.5*F2, gpp=.5*F3, gppp=.5*F4;
    double logx=log(fmax(X,1e-300)/SDMC_COV_XSTAR);

    pgf->G2 = K*X + L*X*X - V;
    pgf->G2_X = K + 2.*L*X;
    pgf->G2_XX = 2.*L;
    pgf->G2_phi = Kp*X + Lp*X*X - Vp;
    pgf->G2_Xphi = Kp + 2.*Lp*X;
    pgf->G2_XXphi = 2.*Lp;
    pgf->G2_phiphi = Kpp*X + Lpp*X*X - Vpp;
    pgf->G2_Xphiphi = Kpp + 2.*Lpp*X;

    pgf->G3_X = -g/X;
    pgf->G3_XX = g/(X*X);
    pgf->G3_XXX = -2.*g/(X*X*X);
    pgf->G3_phi = -gp*logx;
    pgf->G3_Xphi = -gp/X;
    pgf->G3_XXphi = gp/(X*X);
    pgf->G3_phiphi = -gpp*logx;
    pgf->G3_Xphiphi = -gpp/X;
    pgf->G3_phiphiphi = -gppp*logx;

    pgf->G4 = .5*F;
    pgf->DG4 = .5*(F-1.);
    pgf->G4_phi = .5*F1;
    pgf->G4_phiphi = .5*F2;
    pgf->G4_phiphiphi = .5*F3;
  }

  else if(pba->gravity_model_smg == nkgb){'''
replace_once("gravity_smg/gravity_models_smg.c",ganchor,gcode)

p=Path("gravity_smg/gravity_models_smg.c")
s=p.read_text()
start=s.index("int gravity_models_initial_conditions_smg(")
pos=s.index("case propto_omega:",start)
ic=r'''case sdmc_v3_covariant_exact:
    {
      double q,qp,qpp;
      double ph=log(a);
      sdmc_cov_eval(sdmc_cov_QX,ph,&q,&qp,&qpp);
      double Ht=sqrt(2.*exp(q));
      pvecback_integration[pba->index_bi_phi_smg]=ph;
      pvecback_integration[pba->index_bi_phi_prime_smg]=a*Ht;
    }
    break;

  '''
s=s[:pos]+ic+s[pos:]
p.write_text(s)

print("installed exact-logarithmic covariant reconstruction")
print(Path("output/cov_exact_reconstruction_summary.txt").read_text())
