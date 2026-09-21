#!/usr/bin/env python3
from pathlib import Path
import numpy as np, re
from scipy.interpolate import CubicSpline

# Combined best ordinary cosmology + localized braiding target.
AF=0.0715
ZC=5.00
WIDTH=1.426
D0=0.34231919445927034
DFLOOR=0.045
ALENS=0.024
ZLENS=1.5
SIGLENS=0.8
OX=0.6953970160452582
LAMBDA_E=17.925
ZT=17.775
H0_KMS=69.7095092787552

src=Path("output/combined_cov_target_00_background.dat")
if not src.exists():
    raise RuntimeError("target background missing")
lines=src.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and "1:z" in l][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(src)
d={n:a[:,i] for i,n in enumerate(names)}

N=-np.log1p(d["z"])
order=np.argsort(N)
N=N[order]
def vv(name):
    return np.asarray(d[name])[order]

Hraw=vv("H [1/Mpc]")
rhoraw=vv("(.)rho_smg")
praw=vv("(.)p_smg")
H2raw=Hraw*Hraw

# Smooth dimensionless ratios, following the validated r032 covariant replay.
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

Nc=-np.log1p(ZC)
xx=(ph-Nc)/(2.*WIDTH)
S=.5*(1.+np.tanh(xx))
F=np.exp(AF*S)
Fsp=CubicSpline(ph,F)
F1=Fsp(ph,1)
F2=Fsp(ph,2)
F3=Fsp(ph,3)

# Localized EFT braiding deformation, expressed as a covariant G3=g(phi) X.
z=np.exp(-ph)-1.
Delta=ALENS*np.exp(-0.5*((z-ZLENS)/SIGLENS)**2)

# For phi=N and X=H^2/2:
# alpha_B = 2/F[-F_phi/2 + X g].
# Solving alpha_B = -2 alpha_M + Delta gives:
# g = (F*Delta - F_phi)/H^2.
g=(F*Delta-F1)/(H*H)
gsp=CubicSpline(ph,g)
g1=gsp(ph,1)

Dtarget=DFLOOR+D0*S
alphaM=F1/F
alphaB_target=-2.*alphaM+Delta
alphaK_target=Dtarget-1.5*alphaB_target*alphaB_target

# Exact hi_class background identities specialized to
# G4=F/2, G3=g(phi) X, G2=k1 X+k2 X^2-V, phi=N.
R=3.*rho + 2.*g1*X*X + 6.*F1*H*H + 3.*(F-1.)*H*H
P=3.*pre + 2.*g1*X*X - 2.*X*F2 \
   - 3.*(F-1.)*H*H - 2.*(F-1.)*Hpa \
   - 2.*F1*(H*H+Hpa)
C=(R+P)/(2.*X)

# Covariant alpha_K equation:
# F alpha_K = k1 + 6 k2 X - 4 g_phi X + 6 g H^2.
k2=(F*alphaK_target-C+4.*g1*X-6.*g*H*H)/(4.*X)
k1=C-2.*k2*X
V=.5*(R-P)-k2*X*X

# Reconstruction audit on the target trajectory.
k1sp=CubicSpline(ph,k1)
k2sp=CubicSpline(ph,k2)
Vsp=CubicSpline(ph,V)
alphaK_rec=(k1+6.*k2*X-4.*g1*X+6.*g*H*H)/F
alphaB_rec=2.*(-0.5*F1+X*g)/F
Drec=alphaK_rec+1.5*alphaB_rec*alphaB_rec

for arr,name in [(g,"g"),(k1,"k1"),(k2,"k2"),(V,"V"),(Drec,"Drec"),(alphaB_rec,"alphaB_rec")]:
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"nonfinite {name}")

mask100=ph>=-np.log(101.)
print("COMBINED_COV_TARGET_AUDIT",
      "max_abs_alphaB_err_all=",float(np.max(np.abs(alphaB_rec-alphaB_target))),
      "max_abs_D_err_all=",float(np.max(np.abs(Drec-Dtarget))),
      "Dmin_z100=",float(Drec[mask100].min()),
      "Dmin_all=",float(Drec.min()),
      flush=True)

splines={
    "g":gsp,
    "k1":k1sp,
    "k2":k2sp,
    "V":Vsp,
    "lnH":lnHsp,
}

def carr(x):
    return ",".join(f"{v:.17e}" for v in np.asarray(x).ravel())

hp=Path("gravity_smg/sdmc_covariant_braided_table.h")
with hp.open("w") as f:
    f.write("#ifndef SDMC_COVARIANT_BRAIDED_TABLE_H\n#define SDMC_COVARIANT_BRAIDED_TABLE_H\n")
    f.write(f"#define SDMC_CB_N {len(ph)}\n")
    f.write(f"#define SDMC_CB_AF {AF:.17g}\n")
    f.write(f"#define SDMC_CB_ZC {ZC:.17g}\n")
    f.write(f"#define SDMC_CB_WIDTH {WIDTH:.17g}\n")
    f.write("static const double sdmc_cb_x[SDMC_CB_N] = {"+carr(ph)+"};\n")
    for name,sp in splines.items():
        coeff=np.asarray(sp.c.T).reshape(-1)
        f.write(f"static const double sdmc_cb_{name}[4*(SDMC_CB_N-1)] = "+"{"+carr(coeff)+"};\n")
    f.write(r'''
static int sdmc_cb_idx(double x){
  int lo=0,hi=SDMC_CB_N-1;
  if(x<=sdmc_cb_x[0]) return 0;
  if(x>=sdmc_cb_x[SDMC_CB_N-1]) return SDMC_CB_N-2;
  while(hi-lo>1){int m=(lo+hi)/2; if(sdmc_cb_x[m]<=x) lo=m; else hi=m;}
  return lo;
}
static void sdmc_cb_eval(const double *coef,double x,
                         double *y,double *yp,double *ypp,double *yppp){
  int i=sdmc_cb_idx(x);
  double dx=x-sdmc_cb_x[i];
  const double *cc=coef+4*i;
  *y=((cc[0]*dx+cc[1])*dx+cc[2])*dx+cc[3];
  *yp=(3.*cc[0]*dx+2.*cc[1])*dx+cc[2];
  *ypp=6.*cc[0]*dx+2.*cc[1];
  *yppp=6.*cc[0];
}
#endif
''')
print("WROTE",hp,flush=True)

# Freeze the already-established tracker/background constants because the
# covariant gravity model uses parameters_smg for its own shooting parameter.
p=Path("source/background.c"); s=p.read_text()
hs=s.index("static void sdmc_tracker_fluid_target")
he=s.index("static double sdmc_tracker_fluid_w_only",hs)
block=s[hs:he]
repls=[
    (r"double\\s+Ox\\s*=\\s*pba->parameters_smg\\[0\\]\\s*;",f"double Ox = {OX:.17g};"),
    (r"double\\s+lambda_e\\s*=\\s*pba->parameters_smg\\[1\\]\\s*;",f"double lambda_e = {LAMBDA_E:.17g};"),
    (r"double\\s+zt\\s*=\\s*pba->parameters_smg\\[2\\]\\s*;",f"double zt = {ZT:.17g};"),
    (r"double\\s+dNt\\s*=\\s*pba->parameters_smg\\[3\\]\\s*;","double dNt = 0.5;"),
    (r"double\\s+A\\s*=\\s*pba->parameters_smg\\[4\\]\\s*;","double A = 0.01105624999;"),
    (r"double\\s+tauA\\s*=\\s*pba->parameters_smg\\[5\\]\\s*;","double tauA = 0.25;"),
    (r"double\\s+B\\s*=\\s*pba->parameters_smg\\[6\\]\\s*;","double B = 0.01951933685;"),
    (r"double\\s+tauB\\s*=\\s*pba->parameters_smg\\[7\\]\\s*;","double tauB = 1.5;"),
]
for pat,val in repls:
    block,n=re.subn(pat,val,block,count=1)
    if n!=1:
        raise RuntimeError("tracker helper constant pattern missing: "+pat)
s=s[:hs]+block+s[he:]
p.write_text(s)

# Register covariant braided model after the two EFT diagnostic models.
p=Path("include/background.h"); s=p.read_text()
old="    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip, sdmc_v3_independent_kinetic, sdmc_v3_lens_braiding,\n"
new="    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip, sdmc_v3_independent_kinetic, sdmc_v3_lens_braiding, sdmc_covariant_braided_replay,\n"
if old not in s:
    raise RuntimeError("enum anchor missing")
p.write_text(s.replace(old,new,1))

p=Path("gravity_smg/gravity_models_smg.c"); s=p.read_text()
inc='#include "gravity_models_smg.h"\n'
if inc not in s:
    raise RuntimeError("include anchor missing")
s=s.replace(inc,inc+'#include "sdmc_covariant_braided_table.h"\n',1)

anchor='''  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
parser='''  if (strcmp(string1,"sdmc_covariant_braided_replay") == 0) {
     pba->gravity_model_smg = sdmc_covariant_braided_replay;
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
if s.count(anchor)!=1:
    raise RuntimeError("parser anchor count "+str(s.count(anchor)))
s=s.replace(anchor,parser,1)

ganchor='''  else if(pba->gravity_model_smg == galileon){'''
gcode=r'''  else if (pba->gravity_model_smg == sdmc_covariant_braided_replay) {
    double Nc=-log(1.+SDMC_CB_ZC);
    double xx=(phi-Nc)/(2.*SDMC_CB_WIDTH);
    double tt=tanh(xx);
    double uu=1.-tt*tt;
    double w=SDMC_CB_WIDTH;
    double S=.5*(1.+tt);
    double S1=uu/(4.*w);
    double S2=-tt*uu/(4.*w*w);
    double S3=-uu*(1.-3.*tt*tt)/(8.*w*w*w);
    double A=SDMC_CB_AF;
    double F=exp(A*S);
    double F1=F*(A*S1);
    double F2=F*(A*A*S1*S1+A*S2);
    double F3=F*(A*A*A*S1*S1*S1+3.*A*A*S1*S2+A*S3);

    double g,g1,g2,g3;
    double k1,k1p,k1pp,k1ppp;
    double k2,k2p,k2pp,k2ppp;
    double V,Vp,Vpp,Vppp;
    sdmc_cb_eval(sdmc_cb_g,phi,&g,&g1,&g2,&g3);
    sdmc_cb_eval(sdmc_cb_k1,phi,&k1,&k1p,&k1pp,&k1ppp);
    sdmc_cb_eval(sdmc_cb_k2,phi,&k2,&k2p,&k2pp,&k2ppp);
    sdmc_cb_eval(sdmc_cb_V,phi,&V,&Vp,&Vpp,&Vppp);
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
    pgf->G3_phi = g1*X;
    pgf->G3_Xphi = g1;
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
if s.count(ganchor)!=1:
    raise RuntimeError("G anchor count "+str(s.count(ganchor)))
s=s.replace(ganchor,gcode,1)

start=s.index("int gravity_models_initial_conditions_smg(")
pos=s.index("case propto_omega:",start)
ic='''case sdmc_covariant_braided_replay:
    {
      double ph=log(a);
      double lnh,lnhp,lnhpp,lnhppp;
      sdmc_cb_eval(sdmc_cb_lnH,ph,&lnh,&lnhp,&lnhpp,&lnhppp);
      pvecback_integration[pba->index_bi_phi_smg]=ph;
      pvecback_integration[pba->index_bi_phi_prime_smg]=a*exp(lnh);
    }
    break;

  '''
s=s[:pos]+ic+s[pos:]
p.write_text(s)
print("INSTALLED sdmc_covariant_braided_replay",flush=True)
