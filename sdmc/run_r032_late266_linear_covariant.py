#!/usr/bin/env python3
from pathlib import Path
import re
import numpy as np
from scipy.interpolate import CubicSpline

TARGET_INI = Path("output/linear_cov_target.ini")

def ini_value(path, key):
    for raw in Path(path).read_text().splitlines():
        line = raw.split("#",1)[0].strip()
        if not line or "=" not in line:
            continue
        k,v = line.split("=",1)
        if k.strip() == key:
            return v.strip()
    raise RuntimeError(f"missing {key} in {path}")

def vals(key):
    return [float(x.strip()) for x in ini_value(TARGET_INI,key).split(",")]

H0_KMS = float(ini_value(TARGET_INI,"H0"))
KIN = vals("parameters_smg")
EXP = vals("expansion_smg")
if len(KIN) != 6:
    raise RuntimeError(f"expected 6 kinetic params, got {KIN}")
if len(EXP) != 8:
    raise RuntimeError(f"expected 8 expansion params, got {EXP}")
AF,ZC,WIDTH,D0,POWER,DFLOOR = KIN
OX,LAMBDA_E,ZT,DNT,A_LATE,TAUA,B_LATE,TAUB = EXP

def replace_once(path, old, new):
    p=Path(path); s=p.read_text()
    if new in s:
        return
    if s.count(old) != 1:
        raise RuntimeError(f"{path}: anchor count {s.count(old)} for {old[:90]!r}")
    p.write_text(s.replace(old,new,1))

# Read parameterized target background
path=Path("output/linear_cov_target_00_background.dat")
lines=path.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr))
names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
arr=np.loadtxt(path)
d={n:arr[:,i] for i,n in enumerate(names)}

required=["z","H [1/Mpc]","(.)rho_smg","(.)p_smg"]
missing=[x for x in required if x not in d]
if missing:
    raise RuntimeError(f"missing target background columns {missing}; available={names}")

Nraw=-np.log1p(d["z"])
order=np.argsort(Nraw)
N=Nraw[order]
H=np.asarray(d["H [1/Mpc]"])[order]
rho=np.asarray(d["(.)rho_smg"])[order]
pre=np.asarray(d["(.)p_smg"])[order]
keep=np.r_[True,np.diff(N)>1e-12]
N,H,rho,pre=N[keep],H[keep],rho[keep],pre[keep]

lnHsp=CubicSpline(N,np.log(H))
h=lnHsp(N,1)
Hpa=H*H*h
X=.5*H*H

Nc=-np.log1p(ZC)
xx=(N-Nc)/(2.*WIDTH)
tt=np.tanh(xx)
uu=1.-tt*tt
S=.5*(1.+tt)
S1=uu/(4.*WIDTH)
S2=-tt*uu/(4.*WIDTH*WIDTH)
F=np.exp(AF*S)
F1=F*(AF*S1)
F2=F*((AF*S1)**2+AF*S2)

# Accepted manuscript-compatible linear-G3 reconstruction:
#   G4=F/2, G3=g(phi)X, G2=k1(phi)X+k2(phi)X^2-V(phi)
# and exact No-Slip only on the cosmological trajectory:
#   X G3_X = X g = -G4_phi = -F_phi/2.
g=-F1/(H*H)
gsp=CubicSpline(N,g)
g1=gsp(N,1)

Dtarget=DFLOOR+D0*S**POWER
alphaM=F1/F
alphaB=-2.*alphaM
alphaK_target=Dtarget-1.5*alphaB*alphaB

# hi_class background identities used by the accepted linear-G3 reconstruction.
R=3.*rho + 2.*g1*X*X + 6.*F1*H*H + 3.*(F-1.)*H*H
P=3.*pre + 2.*g1*X*X - 2.*X*F2 \
  - 3.*(F-1.)*H*H - 2.*(F-1.)*Hpa \
  - 2.*F1*(H*H+Hpa)
C=(R+P)/(2.*X)
k2=(F*alphaK_target-C+4.*g1*X-6.*g*H*H)/(4.*X)
k1=C-2.*k2*X
V=.5*(R-P)-k2*X*X

for y,name in [(g,"g"),(k1,"k1"),(k2,"k2"),(V,"V")]:
    if not np.all(np.isfinite(y)):
        raise RuntimeError(f"nonfinite reconstructed {name}")

# Reconstructed alpha_K and D closure on target.
alphaK=(k1+6.*k2*X-4.*g1*X+6.*g*H*H)/F
Dcov=alphaK+1.5*alphaB*alphaB
relD=np.max(np.abs((Dcov-Dtarget)/Dtarget))

# Native sound-speed numerator, when available, is used only as an audit.
if "cs2num" in d:
    Ns=np.asarray(d["cs2num"])[order][keep]
    cs2=Ns/Dcov
else:
    cs2=np.full_like(Dcov,np.nan)

# Build spline table for free covariant evolution.
splines={}
for name,y in [("g",g),("k1",k1),("k2",k2),("V",V),("lnH",np.log(H))]:
    splines[name]=CubicSpline(N,y)

def carr(x):
    return ",".join(f"{v:.17e}" for v in np.asarray(x).ravel())

hp=Path("gravity_smg/sdmc_late266_linear_covariant_table.h")
with hp.open("w") as f:
    f.write("#ifndef SDMC_LATE266_LINEAR_COVARIANT_TABLE_H\n#define SDMC_LATE266_LINEAR_COVARIANT_TABLE_H\n")
    f.write(f"#define SDMC_LIN_N {len(N)}\n")
    f.write(f"#define SDMC_LIN_AF {AF:.17e}\n")
    f.write(f"#define SDMC_LIN_ZC {ZC:.17e}\n")
    f.write(f"#define SDMC_LIN_WIDTH {WIDTH:.17e}\n")
    f.write("static const double sdmc_lin_x[SDMC_LIN_N] = {"+carr(N)+"};\n")
    for name,sp in splines.items():
        coeff=np.asarray(sp.c.T).reshape(-1)
        f.write(f"static const double sdmc_lin_{name}[4*(SDMC_LIN_N-1)] = "+"{"+carr(coeff)+"};\n")
    f.write(r"""
static int sdmc_lin_idx(double x){
  int lo=0,hi=SDMC_LIN_N-1;
  if(x<=sdmc_lin_x[0]) return 0;
  if(x>=sdmc_lin_x[SDMC_LIN_N-1]) return SDMC_LIN_N-2;
  while(hi-lo>1){int m=(lo+hi)/2; if(sdmc_lin_x[m]<=x) lo=m; else hi=m;}
  return lo;
}
static void sdmc_lin_eval(const double *coef,double x,
                          double *y,double *yp,double *ypp,double *yppp){
  int i=sdmc_lin_idx(x);
  double dx=x-sdmc_lin_x[i];
  const double *cc=coef+4*i;
  *y=((cc[0]*dx+cc[1])*dx+cc[2])*dx+cc[3];
  *yp=(3.*cc[0]*dx+2.*cc[1])*dx+cc[2];
  *ypp=6.*cc[0]*dx+2.*cc[1];
  *yppp=6.*cc[0];
}
#endif
""")

# Freeze the independent expansion helper constants because the covariant
# gravity model itself has only a potential-offset parameter.
p=Path("source/background.c")
s=p.read_text()
pat=re.compile(
    r"double Ox\s*=\s*pba->parameters_smg\[0\];\s*"
    r"double lambda_e\s*=\s*pba->parameters_smg\[1\];\s*"
    r"double zt\s*=\s*pba->parameters_smg\[2\];\s*"
    r"double dNt\s*=\s*pba->parameters_smg\[3\];\s*"
    r"double A\s*=\s*pba->parameters_smg\[4\];\s*"
    r"double tauA\s*=\s*pba->parameters_smg\[5\];\s*"
    r"double B\s*=\s*pba->parameters_smg\[6\];\s*"
    r"double tauB\s*=\s*pba->parameters_smg\[7\];"
)
repl=f"""double Ox = {OX:.17g};
    double lambda_e = {LAMBDA_E:.17g};
    double zt = {ZT:.17g};
    double dNt = {DNT:.17g};
    double A = {A_LATE:.17g};
    double tauA = {TAUA:.17g};
    double B = {B_LATE:.17g};
    double tauB = {TAUB:.17g};"""
s2,n=pat.subn(repl,s,count=1)
if n != 1:
    raise RuntimeError(f"split tracker helper replacement count={n}")
p.write_text(s2)

replace_once("include/background.h",
             "sdmc_v3_native_noslip,",
             "sdmc_v3_native_noslip, sdmc_v3_covariant_linear_audit,")

replace_once("gravity_smg/gravity_models_smg.c",
             '#include "gravity_models_smg.h"\n',
             '#include "gravity_models_smg.h"\n#include "sdmc_late266_linear_covariant_table.h"\n')

parser_anchor='  if (strcmp(string1,"eft_alphas_power_law") == 0) {'
parser=r'''  if (strcmp(string1,"sdmc_v3_covariant_linear_audit") == 0) {
    pba->gravity_model_smg = sdmc_v3_covariant_linear_audit;
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

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
replace_once("gravity_smg/gravity_models_smg.c",parser_anchor,parser)

ganchor='  else if(pba->gravity_model_smg == galileon){'
gcode=r'''  else if (pba->gravity_model_smg == sdmc_v3_covariant_linear_audit) {
    double Nc=-log(1.+SDMC_LIN_ZC);
    double xx=(phi-Nc)/(2.*SDMC_LIN_WIDTH);
    double tt=tanh(xx);
    double uu=1.-tt*tt;
    double w=SDMC_LIN_WIDTH;
    double S=.5*(1.+tt);
    double S1=uu/(4.*w);
    double S2=-tt*uu/(4.*w*w);
    double S3=-uu*(1.-3.*tt*tt)/(8.*w*w*w);
    double A=SDMC_LIN_AF;
    double F=exp(A*S);
    double F1=F*(A*S1);
    double F2=F*(A*A*S1*S1+A*S2);
    double F3=F*(A*A*A*S1*S1*S1+3.*A*A*S1*S2+A*S3);

    double g,g1,g2,g3;
    double k1,k1p,k1pp,k1ppp;
    double k2,k2p,k2pp,k2ppp;
    double V,Vp,Vpp,Vppp;
    sdmc_lin_eval(sdmc_lin_g,phi,&g,&g1,&g2,&g3);
    sdmc_lin_eval(sdmc_lin_k1,phi,&k1,&k1p,&k1pp,&k1ppp);
    sdmc_lin_eval(sdmc_lin_k2,phi,&k2,&k2p,&k2pp,&k2ppp);
    sdmc_lin_eval(sdmc_lin_V,phi,&V,&Vp,&Vpp,&Vppp);
    V += 3.*pba->H0*pba->H0*pba->parameters_smg[0];

    pgf->G2 = k1*X+k2*X*X-V;
    pgf->G2_X = k1+2.*k2*X;
    pgf->G2_XX = 2.*k2;
    pgf->G2_phi = k1p*X+k2p*X*X-Vp;
    pgf->G2_Xphi = k1p+2.*k2p*X;
    pgf->G2_XXphi = 2.*k2p;
    pgf->G2_phiphi = k1pp*X+k2pp*X*X-Vpp;
    pgf->G2_Xphiphi = k1pp+2.*k2pp*X;

    pgf->G3_X = g;
    pgf->G3_XX = 0.;
    pgf->G3_XXX = 0.;
    pgf->G3_phi = g1*X;
    pgf->G3_Xphi = g1;
    pgf->G3_XXphi = 0.;
    pgf->G3_phiphi = g2*X;
    pgf->G3_Xphiphi = g2;
    pgf->G3_phiphiphi = g3*X;

    pgf->DG4 = (F-1.)/2.;
    pgf->G4 = F/2.;
    pgf->G4_phi = F1/2.;
    pgf->G4_phiphi = F2/2.;
    pgf->G4_phiphiphi = F3/2.;
  }

  else if(pba->gravity_model_smg == galileon){'''
replace_once("gravity_smg/gravity_models_smg.c",ganchor,gcode)

p=Path("gravity_smg/gravity_models_smg.c")
s=p.read_text()
start=s.index("int gravity_models_initial_conditions_smg(")
pos=s.index("case propto_omega:",start)
ic=r'''case sdmc_v3_covariant_linear_audit:
    {
      double ph=log(a);
      double lnh,lnhp,lnhpp,lnhppp;
      sdmc_lin_eval(sdmc_lin_lnH,ph,&lnh,&lnhp,&lnhpp,&lnhppp);
      pvecback_integration[pba->index_bi_phi_smg]=ph;
      pvecback_integration[pba->index_bi_phi_prime_smg]=a*exp(lnh);
    }
    break;

  '''
s=s[:pos]+ic+s[pos:]
p.write_text(s)

i0=int(np.argmin(np.abs(N)))
summary = [
    f"H0_kms={H0_KMS:.15g}",
    f"AF={AF:.15g}",
    f"zc={ZC:.15g}",
    f"width={WIDTH:.15g}",
    f"D0={D0:.15g}",
    f"Dfloor={DFLOOR:.15g}",
    f"Dclosure_max_rel={relD:.12g}",
    f"F0_target={F[i0]:.15g}",
    f"alphaM0_target={alphaM[i0]:.15g}",
    f"alphaB0_target={alphaB[i0]:.15g}",
    f"D0_target={Dcov[i0]:.15g}",
    f"g0={g[i0]:.15g}",
    f"g0_H0sq={g[i0]*(H0_KMS/299792.458)**2:.15g}",
    f"k1_0={k1[i0]:.15g}",
    f"k2X_0={k2[i0]*X[i0]:.15g}",
    f"V_H2_0={V[i0]/(H[i0]*H[i0]):.15g}",
    f"g_range=[{g.min():.9g},{g.max():.9g}]",
    f"k1_range=[{k1.min():.9g},{k1.max():.9g}]",
    f"k2X_range=[{(k2*X).min():.9g},{(k2*X).max():.9g}]",
    f"V_H2_range=[{(V/(H*H)).min():.9g},{(V/(H*H)).max():.9g}]",
]
if np.any(np.isfinite(cs2)):
    summary += [f"cs2_target_min={np.nanmin(cs2):.15g}", f"cs2_target_max={np.nanmax(cs2):.15g}"]
Path("output/linear_cov_reconstruction_summary.txt").write_text("\n".join(summary)+"\n")
print("LINEAR_COVARIANT_RECONSTRUCTION")
print("\n".join(summary))
