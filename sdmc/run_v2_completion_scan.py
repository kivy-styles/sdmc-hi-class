#!/usr/bin/env python3
"""Search stable low-response EFT completions of the frozen NKp-v2 background.

The expansion history, matter densities and early tracker are held fixed.  We
scan a two-parameter deformation of the smooth B0 shape basis,

  F = F_B0**s,  alpha_M = s alpha_M,B0,  alpha_B = beta_B alpha_M,

with a fixed positive kinetic normalization D=D0.  Native hi_class stability
is mandatory.  Every stable point is then run through fully-dynamic
TT/TE/EE+lensing+P(k).  This is a structural search, not a likelihood fit.
"""
from __future__ import annotations
import csv, glob, math, os, subprocess
from pathlib import Path

BG_TEMPLATE = Path("sdmc/config/nkp_v2_full_noslip_background_diagnostic.ini")
FULL_TEMPLATE = Path("sdmc/config/nkp_v2_full_clustered_noslip.ini")
REF = Path("output/sdmc_nkp_v2_tracker_clustered_cs0003_00_cl_lensed.dat")
OUT = Path("output/sdmc_nkp_v2_completion_scan.csv")

SCALES = [0.0,0.10,0.25,0.50,0.75,1.00]
BETAS = [-4.0,-3.0,-2.5,-2.25,-2.0,-1.75,-1.5,-1.0,-0.5,0.0,1.0,2.0]
D0 = 1.0
NS_SEED = 0.001
ZACT = 21.45883112
DNACT = 0.5

FIELDS = ["scale","beta_B","D0","bg_status","full_status","min_cs2","min_D","F0",
          "alphaM0","alphaM_max","z_alphaM_max","sigma8_z0","sigma8_z10",
          "growth_0_over_10","phiphi_ratio_l100","phiphi_ratio_l500",
          "phiphi_ratio_l1000","phiphi_ratio_l1500","TT_ratio_l200",
          "TT_ratio_l1000","TT_ratio_l2000","error"]


def replace_line(text,key,value):
    lines=text.splitlines(); n=0; out=[]
    for line in lines:
        if line.strip().startswith(key+" ="):
            out.append(f"{key} = {value}"); n+=1
        else: out.append(line)
    if n!=1: raise RuntimeError(f"expected one {key}, found {n}")
    return "\n".join(out)+"\n"


def remove_line(text,key):
    lines=text.splitlines(); out=[]
    for line in lines:
        if line.strip().startswith(key+" ="): continue
        out.append(line)
    return "\n".join(out)+"\n"


def table(path):
    rows=[]
    with open(path,encoding="utf-8") as f:
        for line in f:
            s=line.strip()
            if s and not s.startswith("#"): rows.append([float(x) for x in s.split()])
    return rows


def find1(pattern):
    x=sorted(glob.glob(pattern))
    if len(x)!=1: raise RuntimeError(f"{pattern}: found {len(x)}")
    return x[0]


def clean(prefix):
    for p in glob.glob(prefix+"*"):
        try: os.remove(p)
        except OSError: pass


def sigma8(rows):
    vals=[]
    for r in rows:
        k,p=r[0],r[1]; x=8*k
        w=1. if abs(x)<1e-8 else 3*(math.sin(x)-x*math.cos(x))/x**3
        vals.append((math.log(k),k**3*p*w*w/(2*math.pi**2)))
    q=sum(.5*(b[1]+a[1])*(b[0]-a[0]) for a,b in zip(vals[:-1],vals[1:]))
    return math.sqrt(max(q,0.))


def nearest(rows,ell,col): return min(rows,key=lambda r:abs(r[0]-ell))[col]


def bg_metrics(path):
    b=table(path)
    z0=min(b,key=lambda r:abs(r[0])); ammax=max(b,key=lambda r:r[27])
    return dict(min_cs2=min(r[29] for r in b),min_D=min(r[30] for r in b),
                F0=z0[22],alphaM0=z0[27],alphaM_max=ammax[27],z_alphaM_max=ammax[0])


def params(s,b): return f"{s:.12g}, {b:.12g}, {D0:.12g}, {NS_SEED:.12g}, {ZACT:.12g}, {DNACT:.12g}"


def main():
    bg_base=BG_TEMPLATE.read_text(encoding="utf-8")
    bg_base=remove_line(bg_base,"skip_stability_tests_smg")
    full_base=FULL_TEMPLATE.read_text(encoding="utf-8")
    ref=table(REF)
    rows=[]

    for s in SCALES:
      for beta in BETAS:
        row={k:"" for k in FIELDS}; row.update(scale=f"{s:.12g}",beta_B=f"{beta:.12g}",D0=f"{D0:.12g}")
        tag=(f"s{s:.2f}_b{beta:+.2f}").replace(".","p").replace("+","p").replace("-","m")
        prefix=f"output/sdmc_v2_comp_{tag}_"
        ini=Path(f"output/sdmc_v2_comp_{tag}.ini")
        clean(prefix)
        try:
            text=replace_line(bg_base,"gravity_model","sdmc_v2_completion")
            text=replace_line(text,"parameters_smg",params(s,beta))
            text=replace_line(text,"root",prefix)
            for k in ("input_verbose","background_verbose","output_verbose"):
                text=replace_line(text,k,"0")
            ini.write_text(text,encoding="utf-8")
            cp=subprocess.run(["./class",str(ini)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
            if cp.returncode:
                row["bg_status"]="FAIL"; row["full_status"]="SKIP"
                row["error"]=" | ".join(cp.stdout.strip().splitlines()[-12:])[:2000]
                rows.append(row); print(s,beta,"BG_FAIL"); continue

            row["bg_status"]="OK"
            m=bg_metrics(find1(prefix+"*_background.dat"))
            row.update(**{k:f"{v:.12g}" for k,v in m.items()})
            clean(prefix)

            text=replace_line(full_base,"gravity_model","sdmc_v2_completion")
            text=replace_line(text,"parameters_smg",params(s,beta))
            text=replace_line(text,"method_qs_smg","fully_dynamic")
            text=replace_line(text,"z_pk","0,10")
            text=replace_line(text,"root",prefix)
            for k in ("input_verbose","background_verbose","thermodynamics_verbose","perturbations_verbose","output_verbose"):
                text=replace_line(text,k,"0")
            ini.write_text(text,encoding="utf-8")
            cp=subprocess.run(["./class",str(ini)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=240)
            if cp.returncode:
                row["full_status"]="FAIL"
                row["error"]=" | ".join(cp.stdout.strip().splitlines()[-16:])[:2000]
                rows.append(row); print(s,beta,"FULL_FAIL"); continue

            row["full_status"]="OK"
            cl=table(find1(prefix+"*_cl_lensed.dat"))
            p0=table(find1(prefix+"*_z1_pk.dat")); p10=table(find1(prefix+"*_z2_pk.dat"))
            s0,s10=sigma8(p0),sigma8(p10)
            row.update(sigma8_z0=f"{s0:.12g}",sigma8_z10=f"{s10:.12g}",growth_0_over_10=f"{s0/s10:.12g}")
            for ell in (100,500,1000,1500): row[f"phiphi_ratio_l{ell}"]=f"{nearest(cl,ell,5)/nearest(ref,ell,5):.12g}"
            for ell in (200,1000,2000): row[f"TT_ratio_l{ell}"]=f"{nearest(cl,ell,1)/nearest(ref,ell,1):.12g}"
            rows.append(row)
            print(s,beta,"OK","cs2min",row["min_cs2"],"sig8",row["sigma8_z0"],"phi1000",row["phiphi_ratio_l1000"])
        except subprocess.TimeoutExpired as exc:
            row["full_status"]="TIMEOUT"; row["error"]=f"timeout {exc.timeout}s"; rows.append(row)
        finally:
            clean(prefix)
            try: ini.unlink()
            except OSError: pass

    with OUT.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    print("Wrote",OUT)

if __name__=="__main__": main()
