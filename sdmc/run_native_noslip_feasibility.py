#!/usr/bin/env python3
"""Test whether the frozen NKp-v2 split-sector background admits an exact-native
No-Slip B0 completion with a positive prescribed scalar-gradient numerator.

This is a diagnostic only.  It does not change Kp, NKp-v2, the expansion
history, or any observable configuration.

For alpha_B=-2 alpha_M and alpha_H=alpha_T=0, the native hi_class numerator
reduces (up to the tiny numerical No-Slip residual) to

    N_s = -2 (1+alpha_M) h - 3 E_o/F - 2 alpha_M,N,

where h = H_N/H and E_o=(rho_o+p_o)/H^2 is the enthalpy of every non-SMG
sector.  Hence an exact-native B0 shooting equation is

    F_N = alpha_M F,
    alpha_M,N = -1/2 [N_target + 2(1+alpha_M)h + 3 E_o/F].

The earlier recovered B0 ODE instead used

    -2(1+alpha_M)(h+alpha_M) - (3 Omega_m+4 Omega_r)/F,

which differs by both an alpha_M term and, in the split-sector completion, the
tracker-fluid enthalpy.  This script uses h and E_o directly from the frozen
hi_class diagnostic background and asks whether the exact-native ODE can keep
alpha_M(0)=0 while retaining the recovered F(0).
"""
from __future__ import annotations

import bisect
import csv
import math
import re
import sys
from pathlib import Path


def read_background(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    hdr = [x for x in lines if x.startswith("#") and "1:z" in x][-1].lstrip("#").strip()
    marks = list(re.finditer(r"(\d+)\s*:\s*", hdr))
    names = []
    for i, m in enumerate(marks):
        end = marks[i+1].start() if i+1 < len(marks) else len(hdr)
        names.append(hdr[m.end():end].strip())
    rows = []
    for line in lines:
        if line.strip() and not line.startswith("#"):
            vals = [float(v) for v in line.split()]
            rows.append(dict(zip(names, vals)))
    return rows


class FrozenBackground:
    def __init__(self, rows):
        samples = []
        for r in rows:
            z = r["z"]
            if z < 0.0 or z > 1.0e4:
                continue
            N = -math.log1p(z)
            H = r["H [1/Mpc]"]
            h = -1.5*(r["(.)rho_tot"]+r["(.)p_tot"])/(H*H)
            Eo = (r["rho_tot_wo_smg_dbg"]+r["p_tot_wo_smg_dbg"])/(H*H)
            samples.append((N,h,Eo,r["M*^2_smg"],r["M2_running_smg"]))
        samples.sort()
        self.N = [x[0] for x in samples]
        self.h = [x[1] for x in samples]
        self.Eo = [x[2] for x in samples]
        self.F = [x[3] for x in samples]
        self.am = [x[4] for x in samples]
        self.F0 = samples[-1][3]
        self.am0 = samples[-1][4]

    @staticmethod
    def _lerp(xs, ys, x):
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        j = bisect.bisect_right(xs,x)
        x0,x1 = xs[j-1],xs[j]
        y0,y1 = ys[j-1],ys[j]
        t = (x-x0)/(x1-x0)
        return y0+t*(y1-y0)

    def he(self,N):
        return self._lerp(self.N,self.h,N), self._lerp(self.N,self.Eo,N)


def integrate(bg: FrozenBackground, Ns_target: float, zact: float, width: float, dNmax=0.002):
    N0 = max(bg.N[0], -math.log(1.0e4+1.0))
    N1 = 0.0
    nstep = max(1,int(math.ceil((N1-N0)/dNmax)))
    dN = (N1-N0)/nstep
    Nact = -math.log1p(zact)
    F, am = 1.0, 0.0
    min_target = float("inf")

    def rhs(N,F,am):
        h,Eo = bg.he(N)
        Ns_early = -2.0*h-3.0*Eo
        W = 0.5*(1.0-math.tanh((N-Nact)/width))
        Ns_eff = W*Ns_early+(1.0-W)*Ns_target
        Fp = am*F
        amp = -0.5*(Ns_eff+2.0*(1.0+am)*h+3.0*Eo/F)
        return Fp,amp,Ns_eff

    N=N0
    for _ in range(nstep):
        k1F,k1m,n1 = rhs(N,F,am)
        k2F,k2m,n2 = rhs(N+0.5*dN,F+0.5*dN*k1F,am+0.5*dN*k1m)
        k3F,k3m,n3 = rhs(N+0.5*dN,F+0.5*dN*k2F,am+0.5*dN*k2m)
        k4F,k4m,n4 = rhs(N+dN,F+dN*k3F,am+dN*k3m)
        F += dN*(k1F+2*k2F+2*k3F+k4F)/6.0
        am += dN*(k1m+2*k2m+2*k3m+k4m)/6.0
        min_target = min(min_target,n1,n2,n3,n4)
        N += dN
        if (not math.isfinite(F)) or (not math.isfinite(am)) or F <= 0.0 or abs(am) > 1.0e6:
            return float("nan"),float("nan"),min_target
    return F,am,min_target


def find_endpoint_root(bg,Ns_target,width):
    # Search zact logarithmically for alpha_M(0)=0.
    zs = [10.0**(-4.0+i*(math.log10(100.0)+4.0)/80.0) for i in range(81)]
    vals=[]
    for z in zs:
        F,m,_=integrate(bg,Ns_target,z,width)
        vals.append((z,F,m))
    for a,b in zip(vals[:-1],vals[1:]):
        if not (math.isfinite(a[2]) and math.isfinite(b[2])):
            continue
        if a[2] == 0.0:
            return a
        if a[2]*b[2] < 0.0:
            lo,hi=a[0],b[0]
            for _ in range(45):
                mid=math.sqrt(lo*hi)
                _,mm,_=integrate(bg,Ns_target,mid,width)
                _,ml,_=integrate(bg,Ns_target,lo,width)
                if ml*mm <= 0.0:
                    hi=mid
                else:
                    lo=mid
            z=math.sqrt(lo*hi)
            F,m,_=integrate(bg,Ns_target,z,width)
            return z,F,m
    return None


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_native_noslip_feasibility.py BACKGROUND.dat")
    bg = FrozenBackground(read_background(Path(sys.argv[1])))
    Ns_values=[1e-4,1e-3,1e-2]
    widths=[0.05,0.10,0.15,0.20,0.50]
    zacts=[0.05,0.10,0.20,0.50,1.0,2.0,5.0,10.0,20.0,30.0]

    out=Path("output/sdmc_native_noslip_feasibility.csv")
    out.parent.mkdir(exist_ok=True)
    fields=["Ns_target","width","zact","F0","alphaM0","min_target","delta_F0_vs_recovered","delta_alphaM0_vs_zero"]
    rows=[]
    for Ns in Ns_values:
        for w in widths:
            for z in zacts:
                F,m,mint=integrate(bg,Ns,z,w)
                rows.append(dict(Ns_target=Ns,width=w,zact=z,F0=F,alphaM0=m,min_target=mint,
                                 delta_F0_vs_recovered=F-bg.F0 if math.isfinite(F) else float("nan"),
                                 delta_alphaM0_vs_zero=m))
    with out.open("w",newline="",encoding="utf-8") as f:
        cw=csv.DictWriter(f,fieldnames=fields); cw.writeheader(); cw.writerows(rows)

    roots=Path("output/sdmc_native_noslip_endpoint_roots.csv")
    rfields=["Ns_target","width","root_found","zact_root","F0_at_root","recovered_F0","delta_F0","alphaM0_at_root"]
    rrows=[]
    for Ns in Ns_values:
        for w in widths:
            ans=find_endpoint_root(bg,Ns,w)
            if ans is None:
                rrows.append(dict(Ns_target=Ns,width=w,root_found="NO",zact_root="",F0_at_root="",
                                  recovered_F0=bg.F0,delta_F0="",alphaM0_at_root=""))
            else:
                z,F,m=ans
                rrows.append(dict(Ns_target=Ns,width=w,root_found="YES",zact_root=z,F0_at_root=F,
                                  recovered_F0=bg.F0,delta_F0=F-bg.F0,alphaM0_at_root=m))
    with roots.open("w",newline="",encoding="utf-8") as f:
        cw=csv.DictWriter(f,fieldnames=rfields); cw.writeheader(); cw.writerows(rrows)

    print("RECOVERED_ENDPOINT",{"F0":bg.F0,"alphaM0":bg.am0})
    print("EXACT_NATIVE_ENDPOINT_ROOTS")
    for r in rrows:
        print(r)
    print("WROTE",out,roots)


if __name__ == "__main__":
    main()
