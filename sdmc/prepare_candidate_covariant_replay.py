#!/usr/bin/env python3
"""
Build a candidate-specific full covariant replay model from a freshly generated
SDMC target background.

The target MUST have been generated at the same ordinary cosmology and frozen
SDMC structural point that will be used in the replay. This prevents reuse of
the old r032 covariant table after H0/omega_b/omega_c/... are changed.
"""
from pathlib import Path
import argparse
import numpy as np
import re
from scipy.interpolate import CubicSpline

ap=argparse.ArgumentParser()
ap.add_argument("--target",required=True)
ap.add_argument("--OmegaX",required=True,type=float)
ap.add_argument("--lambda-e",dest="lambda_e",default=17.925,type=float)
ap.add_argument("--z-t",dest="zt",default=17.775,type=float)
ap.add_argument("--dNt",default=0.5,type=float)
ap.add_argument("--A",default=0.01105624999,type=float)
ap.add_argument("--tauA",default=0.25,type=float)
ap.add_argument("--B",default=0.01951933685,type=float)
ap.add_argument("--tauB",default=1.5,type=float)
ap.add_argument("--AF",default=0.0715,type=float)
ap.add_argument("--zc",default=5.0,type=float)
ap.add_argument("--width",default=1.426,type=float)
ap.add_argument("--Dfloor",default=0.045,type=float)
ap.add_argument("--D0",default=0.34231919445927034,type=float)
args=ap.parse_args()

src=Path(args.target)
lines=src.read_text().splitlines()
hdr=[l for l in lines if l.startswith('#') and '1:z' in l][-1].lstrip('#').strip()
marks=list(re.finditer(r'(\d+)\s*:\s*',hdr)); names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(src)
d={n:a[:,i] for i,n in enumerate(names)}
N=-np.log1p(d['z'])
order=np.argsort(N); N=N[order]
def vv(name): return np.asarray(d[name])[order]

Hraw=vv('H [1/Mpc]')
rhoraw=vv('(.)rho_smg')
praw=vv('(.)p_smg')
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

Nc=-np.log1p(args.zc)
xx=(ph-Nc)/(2.*args.width)
S=.5*(1.+np.tanh(xx))
F=np.exp(args.AF*S)
Fsp=CubicSpline(ph,F)
F1=Fsp(ph,1); F2=Fsp(ph,2); F3=Fsp(ph,3)

# Exact No-Slip relation on target trajectory.
g=-F1/(H*H)
gsp=CubicSpline(ph,g)
g1=gsp(ph,1)

Dtarget=args.Dfloor+args.D0*S
alphaM=F1/F
alphaB=-2.*alphaM
alphaK_target=Dtarget-1.5*alphaB*alphaB

R=3.*rho + 2.*g1*X*X + 6.*F1*H*H + 3.*(F-1.)*H*H
P=3.*pre + 2.*g1*X*X - 2.*X*F2 \
   - 3.*(F-1.)*H*H - 2.*(F-1.)*Hpa \
   - 2.*F1*(H*H+Hpa)
C=(R+P)/(2.*X)
k2=(F*alphaK_target-C+4.*g1*X-6.*g*H*H)/(4.*X)
k1=C-2.*k2*X
V=.5*(R-P)-k2*X*X

spl={}
for name,y in [('g',g),('k1',k1),('k2',k2),('V',V)]:
    sp=CubicSpline(ph,y)
    spl[name]=(y,sp(ph,1),sp(ph,2),sp(ph,3))

alphaK=(k1+6.*k2*X-4.*spl['g'][1]*X+6.*g*H*H)/F
D=alphaK+1.5*alphaB*alphaB
nssp=CubicSpline(N,vv('cs2num'))
Ns=nssp(ph)
cs2=Ns/D
print("CANDIDATE_COV_TARGET_AUDIT",
      "Dmin_all",float(D.min()),
      "Dmax_all",float(D.max()),
      "cs2min_all",float(cs2.min()),
      "cs2max_all",float(cs2.max()),
      flush=True)

splines={}
for name,y in [('g',g),('k1',k1),('k2',k2),('V',V),('lnH',np.log(H))]:
    splines[name]=CubicSpline(ph,y)

def carr(x):
    return ','.join(f'{v:.17e}' for v in np.asarray(x).ravel())

hp=Path('gravity_smg/sdmc_covariant_replay_table.h')
with hp.open('w') as f:
    f.write('#ifndef SDMC_COVARIANT_REPLAY_TABLE_H\n#define SDMC_COVARIANT_REPLAY_TABLE_H\n')
    f.write(f'#define SDMC_COV_N {len(ph)}\n')
    f.write(f'#define SDMC_COV_AF {args.AF:.17g}\n')
    f.write(f'#define SDMC_COV_ZC {args.zc:.17g}\n')
    f.write(f'#define SDMC_COV_WIDTH {args.width:.17g}\n')
    f.write('static const double sdmc_cov_x[SDMC_COV_N] = {'+carr(ph)+'};\n')
    for name,sp in splines.items():
        coeff=np.asarray(sp.c.T).reshape(-1)
        f.write(f'static const double sdmc_cov_{name}[4*(SDMC_COV_N-1)] = '+'{'+carr(coeff)+'};\n')
    f.write(r'''
static int sdmc_cov_idx(double x){
  int lo=0,hi=SDMC_COV_N-1;
  if(x<=sdmc_cov_x[0]) return 0;
  if(x>=sdmc_cov_x[SDMC_COV_N-1]) return SDMC_COV_N-2;
  while(hi-lo>1){int m=(lo+hi)/2; if(sdmc_cov_x[m]<=x) lo=m; else hi=m;}
  return lo;
}
static void sdmc_cov_eval(const double *coef,double x,
                          double *y,double *yp,double *ypp,double *yppp){
  int i=sdmc_cov_idx(x);
  double dx=x-sdmc_cov_x[i];
  const double *cc=coef+4*i;
  *y=((cc[0]*dx+cc[1])*dx+cc[2])*dx+cc[3];
  *yp=(3.*cc[0]*dx+2.*cc[1])*dx+cc[2];
  *ypp=6.*cc[0]*dx+2.*cc[1];
  *yppp=6.*cc[0];
}
#endif
''')

# Freeze split-tracker background constants to this exact candidate.
p=Path('source/background.c'); s=p.read_text()
old='''  double Ox = pba->parameters_smg[0];
  double lambda_e = pba->parameters_smg[1];
  double zt = pba->parameters_smg[2];
  double dNt = pba->parameters_smg[3];
  double A = pba->parameters_smg[4];
  double tauA = pba->parameters_smg[5];
  double B = pba->parameters_smg[6];
  double tauB = pba->parameters_smg[7];'''
new=f'''  double Ox = {args.OmegaX:.17g};
  double lambda_e = {args.lambda_e:.17g};
  double zt = {args.zt:.17g};
  double dNt = {args.dNt:.17g};
  double A = {args.A:.17g};
  double tauA = {args.tauA:.17g};
  double B = {args.B:.17g};
  double tauB = {args.tauB:.17g};'''
if old not in s:
    raise RuntimeError('split tracker constant anchor missing')
p.write_text(s.replace(old,new,1))

p=Path('include/background.h'); s=p.read_text()
old='    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip,\n'
new='    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip, sdmc_v3_native_noslip, sdmc_covariant_replay,\n'
if old not in s:
    raise RuntimeError('covariant replay enum anchor missing')
p.write_text(s.replace(old,new,1))

p=Path('gravity_smg/gravity_models_smg.c'); s=p.read_text()
s=s.replace('#include "gravity_models_smg.h"\n',
            '#include "gravity_models_smg.h"\n#include "sdmc_covariant_replay_table.h"\n',1)

anchor='''  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
parser='''  if (strcmp(string1,"sdmc_covariant_replay") == 0) {
     pba->gravity_model_smg = sdmc_covariant_replay;
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
    raise RuntimeError('parser anchor count '+str(s.count(anchor)))
s=s.replace(anchor,parser,1)

ganchor='''  else if(pba->gravity_model_smg == galileon){'''
gcode=r'''  else if (pba->gravity_model_smg == sdmc_covariant_replay) {
    double Nc=-log(1.+SDMC_COV_ZC);
    double xx=(phi-Nc)/(2.*SDMC_COV_WIDTH);
    double tt=tanh(xx);
    double uu=1.-tt*tt;
    double w=SDMC_COV_WIDTH;
    double S=.5*(1.+tt);
    double S1=uu/(4.*w);
    double S2=-tt*uu/(4.*w*w);
    double S3=-uu*(1.-3.*tt*tt)/(8.*w*w*w);
    double A=SDMC_COV_AF;
    double F=exp(A*S);
    double F1=F*(A*S1);
    double F2=F*(A*A*S1*S1+A*S2);
    double F3=F*(A*A*A*S1*S1*S1+3.*A*A*S1*S2+A*S3);

    double g,g1,g2,g3;
    double k1,k1p,k1pp,k1ppp;
    double k2,k2p,k2pp,k2ppp;
    double V,Vp,Vpp,Vppp;
    sdmc_cov_eval(sdmc_cov_g,phi,&g,&g1,&g2,&g3);
    sdmc_cov_eval(sdmc_cov_k1,phi,&k1,&k1p,&k1pp,&k1ppp);
    sdmc_cov_eval(sdmc_cov_k2,phi,&k2,&k2p,&k2pp,&k2ppp);
    sdmc_cov_eval(sdmc_cov_V,phi,&V,&Vp,&Vpp,&Vppp);
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
    raise RuntimeError('G anchor count '+str(s.count(ganchor)))
s=s.replace(ganchor,gcode,1)

start=s.index('int gravity_models_initial_conditions_smg(')
pos=s.index('case propto_omega:',start)
ic='''case sdmc_covariant_replay:
    {
      double ph=log(a);
      double lnh,lnhp,lnhpp,lnhppp;
      sdmc_cov_eval(sdmc_cov_lnH,ph,&lnh,&lnhp,&lnhpp,&lnhppp);
      pvecback_integration[pba->index_bi_phi_smg]=ph;
      pvecback_integration[pba->index_bi_phi_prime_smg]=a*exp(lnh);
    }
    break;

  '''
s=s[:pos]+ic+s[pos:]
p.write_text(s)
print("CANDIDATE_COVARIANT_REPLAY_INSTALLED",flush=True)
