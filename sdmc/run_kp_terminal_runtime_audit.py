#!/usr/bin/env python3
from pathlib import Path
import re, subprocess, numpy as np, json

OUT=Path("output/kp_terminal_runtime_audit"); OUT.mkdir(parents=True,exist_ok=True)
CONFIG=Path("sdmc/config/kp_tracker_clustered_cs0003.ini")

def make(tag,enabled,cmin,x1,x2,width):
    txt=CONFIG.read_text()
    txt=re.sub(r"(?m)^\s*root\s*=.*$",f"root = {OUT/(tag+'_')}",txt)
    txt += f"""
kp_terminal_enabled = {1 if enabled else 0}
kp_terminal_Cmin = {cmin}
kp_terminal_x1 = {x1}
kp_terminal_x2 = {x2}
kp_terminal_width = {width}
"""
    p=OUT/(tag+".ini"); p.write_text(txt); return p

def run(tag,enabled,cmin=1.0,x1=.002,x2=.998,width=.02):
    p=make(tag,enabled,cmin,x1,x2,width)
    cp=subprocess.run(["./class",str(p)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=600)
    diag=[l for l in cp.stdout.splitlines() if "KPTERM_RUNTIME_" in l]
    print("KPTERM_AUDIT_CLASS",tag,cp.returncode,"diag_lines",len(diag),flush=True)
    for l in diag[:120]: print(l,flush=True)
    if cp.returncode:
        print(cp.stdout[-4000:],flush=True)
        raise SystemExit(tag+" failed")
    return np.loadtxt(OUT/(tag+"_00_cl_lensed.dat"))

null=run("null",False)
strong=run("strong",True,0.35,.002,.998,.02)
ell=null[:,0]
rec={}
for j,nm in [(1,"TT"),(2,"TE"),(3,"EE")]:
    den=null[:,j]
    good=np.abs(den)>1e-30
    frac=np.zeros_like(den)
    frac[good]=strong[good,j]/den[good]-1.
    rec[nm+"_max_abs_frac"]=float(np.max(np.abs(frac[good])))
    for L in [200,500,1000,1500,2000,2500]:
        i=int(np.argmin(np.abs(ell-L)))
        rec[f"{nm}_ratio_{L}"]=float(strong[i,j]/null[i,j]) if abs(null[i,j])>1e-30 else None
print("KPTERM_AUDIT_SPECTRA",json.dumps(rec,sort_keys=True),flush=True)
(OUT/"kp_terminal_runtime_audit.json").write_text(json.dumps(rec,indent=2))
