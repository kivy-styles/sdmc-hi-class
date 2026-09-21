#!/usr/bin/env python3
from pathlib import Path
import json, re, subprocess
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/edge004_kp_terminal_screen"); OUT.mkdir(parents=True,exist_ok=True)
BASE_INI=Path("input_artifact/unpacked/candidate.ini")
if not BASE_INI.exists(): raise SystemExit("missing exact-covariant candidate.ini")
TCMB=2.7255; CAL_SIGMA=.0025

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

CASES=[
 ("null",0,1.0,0.0207,0.9793,0.012),
 ("kp_fixed",1,0.539812,0.0207,0.9793,0.012),
 ("kp_fixed_narrow",1,0.539812,0.0207,0.9793,0.006),
 ("kp_fixed_wide",1,0.539812,0.0207,0.9793,0.025),
 ("kp_c0995",1,0.995,0.0207,0.9793,0.012),
 ("kp_c099",1,0.99,0.0207,0.9793,0.012),
 ("kp_c0975",1,0.975,0.0207,0.9793,0.012),
]

def setv(t,k,v):
    p=re.compile(rf"(?m)^\s*{re.escape(k)}\s*=.*$")
    if p.search(t): return p.sub(f"{k} = {v}",t)
    return t+"\n"+f"{k} = {v}\n"

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return a,np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    raw,harr,dls=load_cls(path)
    def pc(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return raw,dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                    chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

base=BASE_INI.read_text()
rows=[]; spectra={}
for tag,en,cmin,x1,x2,width in CASES:
    root=str(OUT/(tag+"_"))
    t=base
    for k,v in [
      ("kp_terminal_enabled",en),("kp_terminal_Cmin",cmin),
      ("kp_terminal_x1",x1),("kp_terminal_x2",x2),
      ("kp_terminal_width",width),("root",root),
      ("write background","yes"),("write thermodynamics","yes"),
    ]:
        t=setv(t,k,v)
    ip=OUT/(tag+".ini"); ip.write_text(t)
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=360)
    log=OUT/(tag+".log"); log.write_text(cp.stdout)
    rec=dict(id=tag,enabled=bool(en),Cmin=float(cmin),x1=float(x1),x2=float(x2),
             width=float(width),returncode=int(cp.returncode),status="FAIL")
    clp=Path(root+"00_cl_lensed.dat")
    rec["runtime_pre_count"]=cp.stdout.count("KPTERM_RUNTIME_PRE")
    rec["runtime_x_count"]=cp.stdout.count("KPTERM_RUNTIME_X")
    rec["runtime_terms_count"]=cp.stdout.count("KPTERM_RUNTIME_TERMS")
    rec["runtime_call_count"]=cp.stdout.count("KPTERM_RUNTIME_CALL")
    # Collect observed runtime coordinate / C range from instrumentation.
    xs=[]; Cs=[]; Fks=[]
    for line in cp.stdout.splitlines():
        if "KPTERM_RUNTIME_TERMS" in line:
            m=re.search(r"x=([+\-0-9.eE]+).*?C=([+\-0-9.eE]+).*?Fk=([+\-0-9.eE]+)",line)
            if m:
                xs.append(float(m.group(1))); Cs.append(float(m.group(2))); Fks.append(float(m.group(3)))
    if xs:
        rec.update(runtime_x_min=min(xs),runtime_x_max=max(xs),
                   runtime_C_min=min(Cs),runtime_C_max=max(Cs),
                   runtime_max_abs_Fk=max(abs(x) for x in Fks))
    if cp.returncode==0 and clp.exists():
        raw,sc=pscore(clp); spectra[tag]=raw; rec.update(sc); rec["status"]="OK"
    else:
        rec["error"]=cp.stdout[-1000:].replace("\n"," | ")
    print("KP_EDGE004_POINT",json.dumps(rec,sort_keys=True),flush=True)
    rows.append(rec)

# Spectral ratios versus exact null, before A_planck.
if "null" in spectra:
    b=spectra["null"]; ell=b[:,0].astype(int)
    for rec in rows:
        tag=rec["id"]
        if tag not in spectra or tag=="null": continue
        a=spectra[tag]
        for L in [200,500,1000,1500,2000,2500]:
            i=int(np.argmin(np.abs(ell-L)))
            rec[f"TT_ratio_{L}"]=float(a[i,1]/b[i,1])
            rec[f"TE_ratio_{L}"]=float(a[i,2]/b[i,2]) if b[i,2]!=0 else float("nan")
            rec[f"EE_ratio_{L}"]=float(a[i,3]/b[i,3])
        m=(ell>=30)&(ell<=2500)
        for name,col in [("TT",1),("TE",2),("EE",3)]:
            den=np.maximum(np.abs(b[m,col]),1e-30)
            rec[f"{name}_max_abs_frac"]=float(np.max(np.abs(a[m,col]-b[m,col])/den))

df=pd.DataFrame(rows)
null=df[df.id=="null"].iloc[0]
if null.status!="OK": raise SystemExit("null case failed")
for i,r in df.iterrows():
    if r.status=="OK":
        df.loc[i,"delta_planck_vs_null"]=float(r.chi2_planck-null.chi2_planck)
df.to_csv(OUT/"edge004_kp_terminal_screen.csv",index=False)

ok=df[df.status=="OK"].copy()
best=ok.nsmallest(len(ok),"chi2_planck").to_dict("records")
fixed=ok[ok.id=="kp_fixed"]
summary={
 "null_chi2":float(null.chi2_planck),
 "best":best[0] if best else None,
 "kp_fixed":fixed.iloc[0].to_dict() if len(fixed) else None,
 "n_ok":int(len(ok)),
}
(OUT/"edge004_kp_terminal_screen_summary.json").write_text(json.dumps(summary,indent=2,default=float))
print("KP_EDGE004_BEST",json.dumps(best,sort_keys=True,default=float),flush=True)
print("KP_EDGE004_SUMMARY",json.dumps(summary,sort_keys=True,default=float),flush=True)
