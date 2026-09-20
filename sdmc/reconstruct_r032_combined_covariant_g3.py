#!/usr/bin/env python3
"""
Reconstruct a covariant Horndeski trajectory for the combined SDMC candidate.

Assumptions on the background trajectory:
  phi = N = ln a
  X = H^2 / 2
  G4(phi) = F(phi)/2
  G3(phi,X) = g(phi) X
  G2(phi,X) = k1(phi) X + k2(phi) X^2 - V(phi)

The target background/alphas are taken from the validated parameterized
sdmc_v3_lens_braiding model. We reconstruct g from alpha_B exactly, then solve
for k1,k2,V so that the same background and target alpha_K are reproduced.
"""
from pathlib import Path
import numpy as np, re
from scipy.interpolate import CubicSpline

SRC=Path("output/combined_cov_target_00_background.dat")
if not SRC.exists():
    raise FileNotFoundError(SRC)

lines=SRC.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and "1:z" in l][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(SRC)
d={n:a[:,i] for i,n in enumerate(names)}

N=-np.log1p(d["z"])
order=np.argsort(N)
N=N[order]
def vv(name): return np.asarray(d[name])[order]

Hraw=vv("H [1/Mpc]")
rhoraw=vv("(.)rho_smg")
praw=vv("(.)p_smg")
H2raw=Hraw*Hraw

lnHsp=CubicSpline(N,np.log(Hraw))
rrsp=CubicSpline(N,rhoraw/H2raw)
ppsp=CubicSpline(N,praw/H2raw)

ph=N.copy()
H=np.exp(lnHsp(ph))
h=lnHsp(ph,1)
rho=rrsp(ph)*H*H
pre=ppsp(ph)*H*H
Hpa=H*H*h
X=0.5*H*H

AF=0.0715
zc=5.0
width=1.426
Nc=-np.log1p(zc)
xx=(ph-Nc)/(2.*width)
S=.5*(1.+np.tanh(xx))
F=np.exp(AF*S)
Fsp=CubicSpline(ph,F)
F1=Fsp(ph,1)
F2=Fsp(ph,2)
F3=Fsp(ph,3)

# Target EFT functions from the successful localized-braiding model.
alphaB_target=CubicSpline(N,vv("braiding_smg"))(ph)
alphaM_target=CubicSpline(N,vv("M2_running_smg"))(ph)
Dtarget=CubicSpline(N,vv("kin (D)"))(ph)
alphaK_target=Dtarget-1.5*alphaB_target*alphaB_target

# For G4=F/2, G3=g(phi)X, phi=N:
# alpha_B = (-F_phi + H^2 g)/F.
# Hence g=(F alpha_B + F_phi)/H^2.
g=(F*alphaB_target+F1)/(H*H)
gsp=CubicSpline(ph,g)
g1=gsp(ph,1)

# General background identities for G4=F/2 and G3=g(phi) X.
# The No-Slip formulas used earlier had already substituted
# g=-F_phi/H^2; retain the explicit g terms here.
R=3.*rho + 2.*g1*X*X + 3.*F1*H*H - 6.*X*g*H*H + 3.*(F-1.)*H*H
P=3.*pre + 2.*g1*X*X - 2.*X*F2 - 3.*(F-1.)*H*H - 2.*(F-1.)*Hpa - F1*(2.*H*H+Hpa) + 2.*X*g*Hpa
C=(R+P)/(2.*X)

# Covariant alpha_K equation:
# F alpha_K = k1 + 6 k2 X - 4 g_phi X + 6 g H^2,
# and C=k1+2 k2 X.
k2=(F*alphaK_target-C+4.*g1*X-6.*g*H*H)/(4.*X)
k1=C-2.*k2*X
V=.5*(R-P)-k2*X*X

splines={}
for name,y in [("g",g),("k1",k1),("k2",k2),("V",V),("lnH",np.log(H))]:
    splines[name]=CubicSpline(ph,y)

# Reconstruction audits on the target trajectory.
g1a=splines["g"](ph,1)
alphaB_pred=(-F1+H*H*g)/F
alphaK_pred=(k1+6.*k2*X-4.*g1a*X+6.*g*H*H)/F
D_pred=alphaK_pred+1.5*alphaB_pred*alphaB_pred

mask=ph>=-np.log(101.)
audit={
    "max_abs_alphaB_err_z100":float(np.max(np.abs(alphaB_pred[mask]-alphaB_target[mask]))),
    "max_abs_alphaM_analytic_err_z100":float(np.max(np.abs(F1[mask]/F[mask]-alphaM_target[mask]))),
    "max_abs_D_err_z100":float(np.max(np.abs(D_pred[mask]-Dtarget[mask]))),
    "Dmin_z100":float(D_pred[mask].min()),
    "Dmin_all":float(D_pred.min()),
    "alphaB_peak_abs_z100":float(np.max(np.abs(alphaB_target[mask]))),
}
print("COMBINED_COV_RECON_AUDIT",audit,flush=True)

def carr(x):
    return ",".join(f"{v:.17e}" for v in np.asarray(x).ravel())

hp=Path("gravity_smg/sdmc_combined_covariant_table.h")
with hp.open("w") as f:
    f.write("#ifndef SDMC_COMBINED_COVARIANT_TABLE_H\n#define SDMC_COMBINED_COVARIANT_TABLE_H\n")
    f.write(f"#define SDMC_COMB_COV_N {len(ph)}\n")
    f.write("#define SDMC_COMB_COV_AF 0.0715\n")
    f.write("#define SDMC_COMB_COV_ZC 5.0\n")
    f.write("#define SDMC_COMB_COV_WIDTH 1.426\n")
    f.write("static const double sdmc_comb_cov_x[SDMC_COMB_COV_N] = {"+carr(ph)+"};\n")
    for name,sp in splines.items():
        coeff=np.asarray(sp.c.T).reshape(-1)
        f.write(f"static const double sdmc_comb_cov_{name}[4*(SDMC_COMB_COV_N-1)] = "+"{"+carr(coeff)+"};\n")
    f.write(r'''
static int sdmc_comb_cov_idx(double x){
  int lo=0,hi=SDMC_COMB_COV_N-1;
  if(x<=sdmc_comb_cov_x[0]) return 0;
  if(x>=sdmc_comb_cov_x[SDMC_COMB_COV_N-1]) return SDMC_COMB_COV_N-2;
  while(hi-lo>1){int m=(lo+hi)/2; if(sdmc_comb_cov_x[m]<=x) lo=m; else hi=m;}
  return lo;
}
static void sdmc_comb_cov_eval(const double *coef,double x,
                               double *y,double *yp,double *ypp,double *yppp){
  int i=sdmc_comb_cov_idx(x);
  double dx=x-sdmc_comb_cov_x[i];
  const double *cc=coef+4*i;
  *y=((cc[0]*dx+cc[1])*dx+cc[2])*dx+cc[3];
  *yp=(3.*cc[0]*dx+2.*cc[1])*dx+cc[2];
  *ypp=6.*cc[0]*dx+2.*cc[1];
  *yppp=6.*cc[0];
}
#endif
''')
print("WROTE",hp,flush=True)
