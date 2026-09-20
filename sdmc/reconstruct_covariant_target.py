#!/usr/bin/env python3
from pathlib import Path
import argparse, json, re
import numpy as np
from scipy.interpolate import CubicSpline

ap=argparse.ArgumentParser()
ap.add_argument("--background",required=True)
ap.add_argument("--header-out",default="gravity_smg/sdmc_covariant_replay_table.h")
ap.add_argument("--audit-out",default="output/covariant_target_reconstruction_audit.json")
ap.add_argument("--AF",type=float,default=0.0715)
ap.add_argument("--zc",type=float,default=5.0)
ap.add_argument("--width",type=float,default=1.426)
ap.add_argument("--Dfloor",type=float,default=0.045)
ap.add_argument("--D0",type=float,default=0.34231919445927034)
args=ap.parse_args()

src=Path(args.background)
lines=src.read_text().splitlines()
hdr=[l for l in lines if l.startswith("#") and "1:z" in l][-1].lstrip("#").strip()
marks=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
for i,m in enumerate(marks):
    e=marks[i+1].start() if i+1<len(marks) else len(hdr)
    names.append(hdr[m.end():e].strip())
a=np.loadtxt(src)
d={n:a[:,i] for i,n in enumerate(names)}
N=-np.log1p(d["z"])
order=np.argsort(N); N=N[order]
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

Nc=-np.log1p(args.zc)
xx=(ph-Nc)/(2.*args.width)
S=.5*(1.+np.tanh(xx))
F=np.exp(args.AF*S)
Fsp=CubicSpline(ph,F)
F1=Fsp(ph,1)
F2=Fsp(ph,2)

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
for name,y in [("g",g),("k1",k1),("k2",k2),("V",V)]:
    sp=CubicSpline(ph,y)
    spl[name]=(y,sp(ph,1),sp(ph,2),sp(ph,3))

alphaK=(k1+6.*k2*X-4.*spl["g"][1]*X+6.*g*H*H)/F
D=alphaK+1.5*alphaB*alphaB

nssp=CubicSpline(N,vv("cs2num"))
Ns=nssp(ph)
cs2=Ns/D
mask=ph>=-np.log(101.)

splines={}
for name,y in [("g",g),("k1",k1),("k2",k2),("V",V),("lnH",np.log(H))]:
    splines[name]=CubicSpline(ph,y)

def carr(x):
    return ",".join(f"{v:.17e}" for v in np.asarray(x).ravel())

hp=Path(args.header_out)
hp.parent.mkdir(parents=True,exist_ok=True)
with hp.open("w") as f:
    f.write("#ifndef SDMC_COVARIANT_REPLAY_TABLE_H\n#define SDMC_COVARIANT_REPLAY_TABLE_H\n")
    f.write(f"#define SDMC_COV_N {len(ph)}\n")
    f.write(f"#define SDMC_COV_AF {args.AF:.17e}\n")
    f.write(f"#define SDMC_COV_ZC {args.zc:.17e}\n")
    f.write(f"#define SDMC_COV_WIDTH {args.width:.17e}\n")
    f.write("static const double sdmc_cov_x[SDMC_COV_N] = {"+carr(ph)+"};\n")
    for name,sp in splines.items():
        coeff=np.asarray(sp.c.T).reshape(-1)
        f.write(f"static const double sdmc_cov_{name}[4*(SDMC_COV_N-1)] = "+"{"+carr(coeff)+"};\n")
    f.write(r"""
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
""")

audit={
    "Dmin_z100":float(D[mask].min()),
    "Dmax_z100":float(D[mask].max()),
    "cs2min_z100":float(cs2[mask].min()),
    "cs2max_z100":float(cs2[mask].max()),
    "Dmin_all":float(D.min()),
    "cs2min_all":float(cs2.min()),
    "cs2max_all":float(cs2.max()),
    "AF":args.AF,"zc":args.zc,"width":args.width,
    "Dfloor":args.Dfloor,"D0":args.D0,
}
Path(args.audit_out).parent.mkdir(parents=True,exist_ok=True)
Path(args.audit_out).write_text(json.dumps(audit,indent=2)+"\n")
print("COV_TARGET_AUDIT",audit,flush=True)
print("WROTE",hp,flush=True)
