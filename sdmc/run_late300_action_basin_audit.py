#!/usr/bin/env python3
from pathlib import Path
import sys,re,json,numpy as np

VARIANTS=[
 ("base",0.0,0.0),
 ("phi_p1e4",1e-4,0.0),("phi_m1e4",-1e-4,0.0),
 ("phi_p1e3",1e-3,0.0),("phi_m1e3",-1e-3,0.0),
 ("phi_p1e2",1e-2,0.0),("phi_m1e2",-1e-2,0.0),
 ("vel_p1e4",0.0,1e-4),("vel_m1e4",0.0,-1e-4),
 ("vel_p1e3",0.0,1e-3),("vel_m1e3",0.0,-1e-3),
 ("vel_p1e2",0.0,1e-2),("vel_m1e2",0.0,-1e-2),
 ("both_p1e3",1e-3,1e-3),("both_m1e3",-1e-3,-1e-3),
]

def patch_source():
    p=Path("gravity_smg/gravity_models_smg.c")
    s=p.read_text()
    marker="pba->gravity_model_smg = sdmc_v3_covariant_linear_audit;"
    i=s.index(marker)
    j=s.index('if (strcmp(string1,"eft_alphas_power_law")',i)
    block=s[i:j]
    if "pba->parameters_size_smg = 1;" not in block:
        raise RuntimeError("accepted-action parser size=1 anchor missing")
    block=block.replace("pba->parameters_size_smg = 1;","pba->parameters_size_smg = 3;",1)
    s=s[:i]+block+s[j:]
    old="double ph=log(a);"
    new="double ph=log(a)+pba->parameters_smg[1];"
    if old not in s: raise RuntimeError("phi IC anchor missing")
    s=s.replace(old,new,1)
    old2="pvecback_integration[pba->index_bi_phi_prime_smg]=a*exp(lnh);"
    new2="pvecback_integration[pba->index_bi_phi_prime_smg]=a*exp(lnh)*(1.+pba->parameters_smg[2]);"
    if old2 not in s: raise RuntimeError("phi-prime IC anchor missing")
    s=s.replace(old2,new2,1)
    p.write_text(s)
    print("ACTION_BASIN_PATCHED parameters_size=3")

def prepare():
    base=Path("output/linear_cov_free.ini").read_text()
    out=Path("output"); out.mkdir(exist_ok=True)
    for tag,dphi,dv in VARIANTS:
        txt=base
        txt=txt.replace("parameters_smg = 0.0",f"parameters_smg = 0.0, {dphi:.17g}, {dv:.17g}")
        txt=txt.replace("root = output/linear_cov_free_",f"root = output/basin_{tag}_")
        lines=[]
        for line in txt.splitlines():
            if line.strip().startswith("output = "): continue
            if line.strip().startswith("write_thermodynamics"): line="write_thermodynamics = no"
            if line.strip().startswith("lensing = "): line="lensing = no"
            lines.append(line)
        Path(f"output/basin_{tag}.ini").write_text("\n".join(lines)+"\n")
    Path("output/action_basin_variants.json").write_text(json.dumps(
      [{"tag":t,"delta_phi":p,"delta_velocity_fraction":v} for t,p,v in VARIANTS],
      indent=2,sort_keys=True)+"\n")
    print("ACTION_BASIN_PREPARED",len(VARIANTS))

def read_bg(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    return {n:a[:,i] for i,n in enumerate(names)}

def interp_at(d,name,z):
    N=-np.log1p(d["z"]); o=np.argsort(N)
    return float(np.interp(-np.log1p(z),N[o],d[name][o]))

def analyze():
    outdir=Path("output")
    base=read_bg(outdir/"basin_base_00_background.dat")
    phi_name=next((k for k in base if k.lower().startswith("phi") and "prime" not in k.lower()),None)
    phip_name=next((k for k in base if "phi" in k.lower() and "prime" in k.lower()),None)
    recs=[]
    refs=[1e6,1e4,100.,10.,1.,0.]
    for tag,dphi,dv in VARIANTS:
        rcpath=outdir/f"basin_{tag}.rc"
        rc=int(rcpath.read_text().strip()) if rcpath.exists() else 999
        bgpath=outdir/f"basin_{tag}_00_background.dat"
        r={"tag":tag,"delta_phi":dphi,"delta_velocity_fraction":dv,"returncode":rc}
        if rc!=0 or not bgpath.exists():
            r["status"]="failed"
            recs.append(r); continue
        d=read_bg(bgpath)
        N=-np.log1p(d["z"]); o=np.argsort(N)
        Nb=-np.log1p(base["z"]); ob=np.argsort(Nb)
        H0=np.interp(N[o],Nb[ob],base["H [1/Mpc]"][ob])
        F0=np.interp(N[o],Nb[ob],base["M*^2_smg"][ob])
        hrel=np.abs(d["H [1/Mpc]"][o]/H0-1)
        frel=np.abs(d["M*^2_smg"][o]/F0-1)
        z=d["z"][o]
        m100=z<=100
        slip=np.abs(d["braiding_smg"][o]+2*d["M2_running_smg"][o])
        r.update({
          "status":"success",
          "max_abs_dH_H_z100":float(np.max(hrel[m100])),
          "max_abs_dF_F_z100":float(np.max(frel[m100])),
          "min_D_z100":float(np.min(d["kin (D)"][o][m100])),
          "min_cs2_z100":float(np.min(d["c_s^2"][o][m100])),
          "max_cs2_z100":float(np.max(d["c_s^2"][o][m100])),
          "max_abs_noslip_z100":float(np.max(slip[m100])),
          "snapshots":{}
        })
        for zr in refs:
            q={
              "dH_H":interp_at(d,"H [1/Mpc]",zr)/interp_at(base,"H [1/Mpc]",zr)-1,
              "dF_F":interp_at(d,"M*^2_smg",zr)/interp_at(base,"M*^2_smg",zr)-1,
            }
            if phi_name:
                q["dphi"]=interp_at(d,phi_name,zr)-interp_at(base,phi_name,zr)
            if phip_name:
                bp=interp_at(base,phip_name,zr)
                vp=interp_at(d,phip_name,zr)
                q["dphi_prime_rel"]=vp/bp-1 if bp!=0 else float("nan")
            r["snapshots"][str(zr)]=q
        # Conservative operational classification, reported transparently.
        r["late300_recovery_gate"]=bool(
          r["min_D_z100"]>0 and r["min_cs2_z100"]>0 and r["max_cs2_z100"]<=1.0+1e-8
          and abs(r["snapshots"]["0.0"]["dH_H"])<1e-3
          and abs(r["snapshots"]["0.0"]["dF_F"])<1e-3
        )
        recs.append(r)
    successful=[r for r in recs if r.get("status")=="success"]
    passing=[r for r in successful if r.get("late300_recovery_gate")]
    out={
      "status":"accepted frozen-action initial-condition basin audit",
      "definition":{
        "delta_phi":"additive perturbation to initial phi=ln(a)",
        "delta_velocity_fraction":"fractional perturbation to initial phi_prime=aH",
        "action_retuned":False,
        "potential_offset":0.0,
        "ordinary_cosmology_changed":False
      },
      "n_variants":len(recs),"n_successful":len(successful),"n_recovery_gate":len(passing),
      "all_recovery_gate":len(passing)==len(recs),
      "phi_column":phi_name,"phi_prime_column":phip_name,
      "variants":recs
    }
    Path("output/late300_action_basin_audit.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("LATE300_ACTION_BASIN_AUDIT",json.dumps(out,sort_keys=True))

if __name__=="__main__":
    if len(sys.argv)!=2: raise SystemExit("usage: ... patch|prepare|analyze")
    {"patch":patch_source,"prepare":prepare,"analyze":analyze}[sys.argv[1]]()
