#!/usr/bin/env python3
from pathlib import Path
import re, numpy as np, json

def parse_bg(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)},names

def arr(name,x):
    vals=",".join(f"{v:.17e}" for v in x)
    return f"static const double {name}[L300_N] = {{{vals}}};\n"

def main():
    src=Path("inputs/late300_action/linear_cov_free_00_background.dat")
    d,names=parse_bg(src)
    z=np.asarray(d["z"],float); a=1/(1+z)
    H=np.asarray(d["H [1/Mpc]"],float)
    key=lambda s: np.asarray(d[s],float)
    order=np.argsort(a); a=a[order]; H=H[order]
    alphaK=key("kineticity_smg")[order]
    alphaB=key("braiding_smg")[order]
    alphaM=key("M2_running_smg")[order]
    M2=key("M*^2_smg")[order]
    Geff=key("G_eff_smg")[order]
    # ReACT spherical collapse starts at AMIN=3e-5.  The interpolation table must
    # cover that range; holding H constant below 1e-4 badly distorts collapse.
    amin=max(float(a.min()),2.5e-5)
    aa=np.geomspace(amin,1.0,900)
    def I(y): return np.interp(aa,a,y)
    Hs=I(H)
    Efull=Hs/Hs[-1]

    # ReACT's late-time spherical-collapse equations do not include radiation in
    # the stock LCDM/EFT backgrounds.  Strip the known photon+massless-neutrino
    # contribution from the hi_class background before supplying H(a) to collapse.
    # This leaves the accepted-action matter+structural background and avoids
    # importing the radiation-era H~a^-2 behaviour into a matter-era halo model.
    c_km_s=299792.458
    h0=float(Hs[-1]*c_km_s/100.0)
    omega_gamma=2.469e-5
    N_ur=3.046
    omega_r=omega_gamma*(1.0+0.22710731766*N_ur)
    Omega_r0=omega_r/(h0*h0)
    esc2=(Efull*Efull-Omega_r0/aa**4)/(1.0-Omega_r0)
    if np.any(esc2<=0):
        raise RuntimeError(f"radiation-stripped E^2 became non-positive: min={esc2.min()}")
    E=np.sqrt(esc2)
    ak,ab,am,m2,ge=[I(x) for x in (alphaK,alphaB,alphaM,M2,Geff)]
    dak=np.gradient(ak,aa,edge_order=2); dab=np.gradient(ab,aa,edge_order=2); dam=np.gradient(am,aa,edge_order=2)
    d2ak=np.gradient(dak,aa,edge_order=2); d2ab=np.gradient(dab,aa,edge_order=2); d2am=np.gradient(dam,aa,edge_order=2)
    dE=np.gradient(E,aa,edge_order=2)
    HA1=aa*E*dE
    out=Path("output/react_tables"); out.mkdir(parents=True,exist_ok=True)
    h=out/"late300_table.h"
    text=f"""#pragma once
#include <cmath>
static const int L300_N = {len(aa)};
""" + arr("L300_A",aa)+arr("L300_E",E)+arr("L300_HA1",HA1)+arr("L300_AK",ak)+arr("L300_AB",ab)+arr("L300_AM",am)+arr("L300_M2",m2)+arr("L300_GE",ge)+arr("L300_DAK",dak)+arr("L300_DAB",dab)+arr("L300_DAM",dam)+arr("L300_D2AK",d2ak)+arr("L300_D2AB",d2ab)+arr("L300_D2AM",d2am)+"""
inline double l300_interp(double x,const double *y){
  if(x<=L300_A[0]) return y[0];
  if(x>=L300_A[L300_N-1]) return y[L300_N-1];
  int lo=0,hi=L300_N-1;
  while(hi-lo>1){ int m=(lo+hi)/2; if(L300_A[m]<=x) lo=m; else hi=m; }
  double t=(x-L300_A[lo])/(L300_A[hi]-L300_A[lo]);
  return y[lo]+t*(y[hi]-y[lo]);
}
inline double late300_E(double a){return l300_interp(a,L300_E);}
inline double late300_HA1(double a){return l300_interp(a,L300_HA1);}
inline double late300_ak(double a){return l300_interp(a,L300_AK);}
inline double late300_ab(double a){return l300_interp(a,L300_AB);}
inline double late300_am(double a){return l300_interp(a,L300_AM);}
inline double late300_m2(double a){return l300_interp(a,L300_M2);}
inline double late300_ge(double a){return l300_interp(a,L300_GE);}
inline double late300_dak(double a){return l300_interp(a,L300_DAK);}
inline double late300_dab(double a){return l300_interp(a,L300_DAB);}
inline double late300_dam(double a){return l300_interp(a,L300_DAM);}
inline double late300_d2ak(double a){return l300_interp(a,L300_D2AK);}
inline double late300_d2ab(double a){return l300_interp(a,L300_D2AB);}
inline double late300_d2am(double a){return l300_interp(a,L300_D2AM);}
"""
    h.write_text(text)
    ztest=np.array([0,0.25,0.5,1,1.5,2,2.5])
    at=1/(1+ztest)
    diag=[]
    for zz,x in zip(ztest,at):
        diag.append(dict(z=float(zz),a=float(x),
                         E=float(np.interp(x,aa,E)),
                         E_full=float(np.interp(x,aa,Efull)),
                         Omega_r0=float(Omega_r0),
                         alphaK=float(np.interp(x,aa,ak)),alphaB=float(np.interp(x,aa,ab)),
                         alphaM=float(np.interp(x,aa,am)),M2=float(np.interp(x,aa,m2)),
                         Geff=float(np.interp(x,aa,ge))))
    (out/"late300_table_diagnostic.json").write_text(json.dumps(diag,indent=2)+"\n")
    print("LATE300_REACT_TABLE",json.dumps(diag,sort_keys=True))
if __name__=="__main__": main()
