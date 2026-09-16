#!/usr/bin/env python3
"""Deep solver diagnostics for the NKp-v2 No-Slip perturbation sector.

This script keeps the frozen NKp-v2 expansion history fixed and probes three
questions:

1. Does hi_class's built-in quasi-static (QS) machinery remove the very large
   late-time fully-dynamic response?
2. What growth would the standard No-Slip QS limit, mu=1/F, predict on the same
   background, independent of the propagating scalar solver?
3. Is the fully-dynamic excess strongly controlled by the width of the smooth
   No-Slip activation/handoff once alpha_M(0)=0 is re-shot for each width?

These are diagnostics, not likelihood fits. Failed or obviously pathological
QS runs are retained explicitly rather than interpreted as physical model
rejections.
"""
from __future__ import annotations
import bisect
import csv
import glob
import math
import os
import subprocess
from pathlib import Path

TEMPLATE=Path("sdmc/config/nkp_v2_full_clustered_noslip.ini")
DIAG_TEMPLATE=Path("sdmc/config/nkp_v2_full_noslip_background_diagnostic.ini")
REF=Path("output/sdmc_nkp_v2_tracker_clustered_cs0003_00_cl_lensed.dat")
BASE_BG=Path("output/sdmc_nkp_v2_full_noslip_00_background.dat")
BASE_PK0=Path("output/sdmc_nkp_v2_full_noslip_00_z1_pk.dat")
BASE_PK10=Path("output/sdmc_nkp_v2_full_noslip_00_z7_pk.dat")
OUT=Path("output/sdmc_nkp_v2_noslip_qs_scan.csv")
METHODS=["fully_dynamic","automatic","quasi_static","quasi_static_debug"]
AUTO_THRESHOLDS=[1e3,1e2,3e1,1e1,3.,1.,3e-1,1e-1,3e-2,1e-2,3e-3,1e-3]
WIDTHS=[0.25,0.35,0.50,0.75,1.00]

FIELDS=[
    "diagnostic","case","method","trigger","width","zstart","returncode","status",
    "sigma8_z0","sigma8_z10","growth_0_over_10",
    "qs_growth_mu1","qs_growth_mu1_over_F","fd_over_qs_mu1_over_F",
    "phiphi_ratio_l100","phiphi_ratio_l500","phiphi_ratio_l1000","phiphi_ratio_l1500",
    "TT_ratio_l200","TT_ratio_l1000","TT_ratio_l2000",
    "min_cs2","min_D","F0","alphaM0","alphaM_max","z_alphaM_max","error"
]


def blank_row(**kw):
    r={k:"" for k in FIELDS}; r.update(kw); return r


def replace_line(text,key,value):
    lines=text.splitlines(); out=[]; n=0
    for line in lines:
        if line.strip().startswith(key+" ="):
            out.append(f"{key} = {value}"); n+=1
        else: out.append(line)
    if n!=1: raise RuntimeError(f"expected one {key}, found {n}")
    return "\n".join(out)+"\n"


def set_or_insert(text,key,value,before_key="method_qs_smg"):
    lines=text.splitlines()
    hit=[i for i,line in enumerate(lines) if line.strip().startswith(key+" =")]
    if len(hit)==1:
        lines[hit[0]]=f"{key} = {value}"
    elif len(hit)==0:
        pos=next((i for i,line in enumerate(lines) if line.strip().startswith(before_key+" =")),None)
        if pos is None:
            lines.append(f"{key} = {value}")
        else:
            lines.insert(pos,f"{key} = {value}")
    else:
        raise RuntimeError(f"expected at most one {key}, found {len(hit)}")
    return "\n".join(lines)+"\n"


def set_parameters(text,zstart,width):
    value=f"0.001, 0.003, 0.1, {zstart:.12g}, 21.10, {width:.12g}"
    return replace_line(text,"parameters_smg",value)


def table(path):
    rows=[]
    with open(path,encoding="utf-8") as f:
        for line in f:
            s=line.strip()
            if s and not s.startswith("#"):
                rows.append([float(x) for x in s.split()])
    return rows


def find1(pattern):
    x=sorted(glob.glob(pattern))
    if len(x)!=1: raise RuntimeError(f"{pattern}: found {len(x)}")
    return x[0]


def sigma8(rows):
    v=[]
    for r in rows:
        k,p=r[0],r[1]; x=8*k
        w=1. if abs(x)<1e-8 else 3*(math.sin(x)-x*math.cos(x))/x**3
        v.append((math.log(k),k**3*p*w*w/(2*math.pi**2)))
    q=sum(.5*(b[1]+a[1])*(b[0]-a[0]) for a,b in zip(v[:-1],v[1:]))
    return math.sqrt(max(q,0.))


def nearest(rows,ell,col): return min(rows,key=lambda r:abs(r[0]-ell))[col]


def clean(prefix):
    for p in glob.glob(prefix+"*"):
        try: os.remove(p)
        except OSError: pass


def background_metrics(path):
    b=table(path)
    # background.dat columns (0 based): z=0, M*^2=22, alpha_M=27, c_s^2=29, D=30
    z0=min(b,key=lambda r:abs(r[0]))
    ammax=max(b,key=lambda r:r[27])
    return {
        "min_cs2":min(r[29] for r in b),
        "min_D":min(r[30] for r in b),
        "F0":z0[22],
        "alphaM0":z0[27],
        "alphaM_max":ammax[27],
        "z_alphaM_max":ammax[0],
    }


def fill_observables(row,prefix,ref):
    cl=table(find1(prefix+"*_cl_lensed.dat"))
    p0=table(find1(prefix+"*_z1_pk.dat"))
    p10=table(find1(prefix+"*_z2_pk.dat"))
    s0,s10=sigma8(p0),sigma8(p10)
    row.update(sigma8_z0=f"{s0:.12g}",sigma8_z10=f"{s10:.12g}",growth_0_over_10=f"{s0/s10:.12g}")
    for e in (100,500,1000,1500): row[f"phiphi_ratio_l{e}"]=f"{nearest(cl,e,5)/nearest(ref,e,5):.12g}"
    for e in (200,1000,2000): row[f"TT_ratio_l{e}"]=f"{nearest(cl,e,1)/nearest(ref,e,1):.12g}"
    return s0/s10


def run_case(base,ref,label,method,trigger=None):
    prefix=f"output/sdmc_nkp_v2_noslip_qs_{label}_"
    clean(prefix)
    ini=Path(f"output/sdmc_nkp_v2_noslip_qs_{label}.ini")
    log=Path(f"output/sdmc_nkp_v2_noslip_qs_{label}.log")
    text=replace_line(base,"method_qs_smg",method)
    if trigger is not None:
        text=set_or_insert(text,"z_fd_qs_smg","0")
        text=set_or_insert(text,"trigger_mass_qs_smg",f"{trigger:.12g}")
        text=set_or_insert(text,"trigger_rad_qs_smg",f"{trigger:.12g}")
    text=replace_line(text,"z_pk","0,10")
    text=replace_line(text,"root",prefix)
    for k in ("input_verbose","background_verbose","thermodynamics_verbose","output_verbose"):
        text=replace_line(text,k,"0")
    text=replace_line(text,"perturbations_verbose","2" if ("quasi_static" in method or trigger is not None) else "0")
    ini.write_text(text,encoding="utf-8")

    row=blank_row(diagnostic="hi_class_qs",case=label,method=method,trigger="" if trigger is None else f"{trigger:.12g}")
    try:
        cp=subprocess.run(["./class",str(ini)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=300)
        row["returncode"]=str(cp.returncode)
        log.write_text(cp.stdout,encoding="utf-8")
        if cp.returncode:
            row["status"]="FAIL"
            tail=[x for x in cp.stdout.strip().splitlines() if x.strip()][-24:]
            row["error"]=("returncode="+str(cp.returncode)+(" | "+" | ".join(tail) if tail else ""))[:4000]
            print(label,"FAIL",row["error"])
            return row
        row["status"]="OK"
        fill_observables(row,prefix,ref)
        print(label,"OK","sigma8",row["sigma8_z0"],"growth",row["growth_0_over_10"],"phi1000/ref",row["phiphi_ratio_l1000"])
        return row
    except subprocess.TimeoutExpired as exc:
        row["status"]="TIMEOUT"; row["error"]=f"timeout after {exc.timeout} s"
        print(label,"TIMEOUT")
        return row
    finally:
        clean(prefix)
        try: ini.unlink()
        except OSError: pass


def lininterp(xs,ys,x):
    if x<=xs[0]: return ys[0]
    if x>=xs[-1]: return ys[-1]
    i=bisect.bisect_right(xs,x)-1
    t=(x-xs[i])/(xs[i+1]-xs[i])
    return ys[i]*(1.-t)+ys[i+1]*t


def qs_growth(background_path,use_mu_F=True,zini=100.,ztarget=10.):
    b=table(background_path)
    pts=[]
    for r in b:
        N=-math.log1p(r[0])
        pts.append((N,math.log(r[3]),r[17],r[22])) # ln H, Omega_m, F
    pts.sort()
    xs=[p[0] for p in pts]; lnH=[p[1] for p in pts]; Om=[p[2] for p in pts]; F=[p[3] for p in pts]

    def coeff(N):
        eps=1.e-5
        h=(lininterp(xs,lnH,N+eps)-lininterp(xs,lnH,N-eps))/(2.*eps)
        om=lininterp(xs,Om,N)
        mu=1./lininterp(xs,F,N) if use_mu_F else 1.
        return h,om,mu

    def deriv(N,D,V):
        h,om,mu=coeff(N)
        return V,-(2.+h)*V+1.5*mu*om*D

    Ni=-math.log1p(zini); N0=0.; Nt=-math.log1p(ztarget)
    nstep=30000
    dN=(N0-Ni)/nstep
    D=1.; V=1.; N=Ni; Dt=None
    for _ in range(nstep):
        if Dt is None and N>=Nt: Dt=D
        k1D,k1V=deriv(N,D,V)
        k2D,k2V=deriv(N+.5*dN,D+.5*dN*k1D,V+.5*dN*k1V)
        k3D,k3V=deriv(N+.5*dN,D+.5*dN*k2D,V+.5*dN*k2V)
        k4D,k4V=deriv(N+dN,D+dN*k3D,V+dN*k3V)
        D += dN*(k1D+2*k2D+2*k3D+k4D)/6.
        V += dN*(k1V+2*k2V+2*k3V+k4V)/6.
        N += dN
    if Dt is None: Dt=D
    return D/Dt


def shoot_alphaM0(diag_base,width):
    prefix="output/sdmc_nkp_v2_width_shoot_"
    ini=Path("output/sdmc_nkp_v2_width_shoot.ini")

    def evaluate(zstart):
        clean(prefix)
        text=set_parameters(diag_base,zstart,width)
        text=replace_line(text,"root",prefix)
        for key in ("input_verbose","background_verbose","output_verbose"):
            text=replace_line(text,key,"0")
        ini.write_text(text,encoding="utf-8")
        cp=subprocess.run(["./class",str(ini)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
        if cp.returncode:
            raise RuntimeError("shoot run failed: "+" | ".join(cp.stdout.strip().splitlines()[-12:]))
        bg=find1(prefix+"*_background.dat")
        b=table(bg)
        am0=min(b,key=lambda r:abs(r[0]))[27]
        clean(prefix)
        return am0

    try:
        grid=[2.,4.,6.,8.,10.,12.,15.,18.,21.,24.,28.,34.,42.,55.,75.,100.,140.]
        vals=[]
        for z in grid:
            try: vals.append((z,evaluate(z)))
            except Exception: continue
        bracket=None
        for a,b in zip(vals[:-1],vals[1:]):
            if a[1]==0. or a[1]*b[1]<0.:
                bracket=(a,b); break
        if bracket is None:
            raise RuntimeError("no alphaM0 sign-change bracket; samples="+str(vals))
        lo,flo=bracket[0]; hi,fhi=bracket[1]
        for _ in range(28):
            mid=.5*(lo+hi); fm=evaluate(mid)
            if abs(fm)<2.e-8: return mid,fm
            if flo*fm<=0.: hi,fhi=mid,fm
            else: lo,flo=mid,fm
        mid=.5*(lo+hi); fm=evaluate(mid)
        return mid,fm
    finally:
        clean(prefix)
        try: ini.unlink()
        except OSError: pass


def run_width_case(base,diag_base,ref,width):
    label=(f"width_{width:.2f}").replace(".","p")
    row=blank_row(diagnostic="activation_width",case=label,method="fully_dynamic",width=f"{width:.12g}")
    try:
        zstart,amroot=shoot_alphaM0(diag_base,width)
        row["zstart"]=f"{zstart:.12g}"
        prefix=f"output/sdmc_nkp_v2_noslip_{label}_"
        clean(prefix)
        ini=Path(f"output/sdmc_nkp_v2_noslip_{label}.ini")
        log=Path(f"output/sdmc_nkp_v2_noslip_qs_{label}.log")
        text=set_parameters(base,zstart,width)
        text=replace_line(text,"method_qs_smg","fully_dynamic")
        text=replace_line(text,"z_pk","0,10")
        text=replace_line(text,"root",prefix)
        for key in ("input_verbose","background_verbose","thermodynamics_verbose","perturbations_verbose","output_verbose"):
            text=replace_line(text,key,"0")
        ini.write_text(text,encoding="utf-8")
        cp=subprocess.run(["./class",str(ini)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=300)
        row["returncode"]=str(cp.returncode)
        log.write_text(cp.stdout,encoding="utf-8")
        if cp.returncode:
            row["status"]="FAIL"
            row["error"]=(f"shot alphaM0={amroot:.4g} | "+" | ".join(cp.stdout.strip().splitlines()[-20:]))[:4000]
            return row
        row["status"]="OK"
        fill_observables(row,prefix,ref)
        m=background_metrics(find1(prefix+"*_background.dat"))
        for k,v in m.items(): row[k]=f"{v:.12g}"
        print(label,"OK","zstart",zstart,"sigma8",row["sigma8_z0"],"growth",row["growth_0_over_10"],"phi1000/ref",row["phiphi_ratio_l1000"])
        return row
    except Exception as exc:
        row["status"]="FAIL"; row["error"]=str(exc)[:4000]
        print(label,"FAIL",row["error"])
        return row
    finally:
        try:
            clean(f"output/sdmc_nkp_v2_noslip_{label}_")
            Path(f"output/sdmc_nkp_v2_noslip_{label}.ini").unlink()
        except OSError: pass


def main():
    base=TEMPLATE.read_text(encoding="utf-8")
    diag_base=DIAG_TEMPLATE.read_text(encoding="utf-8")
    ref=table(REF)
    rows=[]

    for method in METHODS:
        rows.append(run_case(base,ref,method,method))
    for t in AUTO_THRESHOLDS:
        tag=(f"{t:.0e}").replace("+","").replace("-","m")
        rows.append(run_case(base,ref,"auto_t"+tag,"automatic",t))

    # Independent standard No-Slip QS growth equation on the exact same H(N),
    # using mu=1/F. Integrate from z=100 and normalize at z=10 so the result
    # is insensitive to the arbitrary initial amplitude.
    g_mu1=qs_growth(BASE_BG,False)
    g_noslip=qs_growth(BASE_BG,True)
    fd_growth=sigma8(table(BASE_PK0))/sigma8(table(BASE_PK10))
    row=blank_row(
        diagnostic="analytic_qs_growth",case="mu_equals_1_over_F",method="growth_ode",status="OK",
        growth_0_over_10=f"{fd_growth:.12g}",qs_growth_mu1=f"{g_mu1:.12g}",
        qs_growth_mu1_over_F=f"{g_noslip:.12g}",fd_over_qs_mu1_over_F=f"{fd_growth/g_noslip:.12g}",
        error="Growth ODE: D_NN+(2+dlnH/dN)D_N-(3/2)mu Omega_m D=0; normalized at z=10 after start at z=100."
    )
    rows.append(row)
    print("analytic QS growth: mu=1",g_mu1,"mu=1/F",g_noslip,"FD",fd_growth,"FD/QS",fd_growth/g_noslip)

    for w in WIDTHS:
        rows.append(run_width_case(base,diag_base,ref,w))

    with OUT.open("w",newline="",encoding="utf-8") as f:
        wr=csv.DictWriter(f,fieldnames=FIELDS); wr.writeheader(); wr.writerows(rows)
    print("Wrote",OUT)

if __name__=="__main__": main()
