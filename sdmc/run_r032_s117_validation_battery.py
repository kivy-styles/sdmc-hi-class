#!/usr/bin/env python3
from pathlib import Path
import json,re,math
import numpy as np
import pandas as pd

OUT=Path("output/s117_validation")
OUT.mkdir(parents=True,exist_ok=True)
BAO=Path("bao_data")
NEFF=3.046
ALPHA_NU=0.22710731766

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def interp(df,col,z):
    x=df["z"].to_numpy(); y=df[col].to_numpy(); o=np.argsort(x)
    return float(np.interp(z,x[o],y[o]))

def find_drag_z(th):
    z=th["z"].to_numpy(float)
    # CLASS thermodynamics output uses kappa_b for baryon drag optical depth.
    cands=[c for c in th.columns if c.strip()=="kappa_b" or c.strip().startswith("kappa_b")]
    if not cands: raise RuntimeError(f"kappa_b not found in thermodynamics columns: {list(th.columns)}")
    k=th[cands[0]].to_numpy(float)
    roots=[]
    for i in range(len(z)-1):
        a,b=k[i]-1.0,k[i+1]-1.0
        if a==0: roots.append(z[i])
        elif a*b<0:
            t=(1-k[i])/(k[i+1]-k[i])
            roots.append(z[i]+t*(z[i+1]-z[i]))
    if not roots:
        i=int(np.argmin(np.abs(k-1.0))); return float(z[i])
    # physical drag root is around z~10^3
    return float(min(roots,key=lambda q:abs(q-1059.0)))

def sigma8_from_pk(path):
    a=np.loadtxt(path)
    k=a[:,0]; P=a[:,1]
    # CLASS pk.dat is k[h/Mpc], P[(Mpc/h)^3] in this workflow.
    R=8.0
    x=k*R
    W=np.ones_like(x)
    m=np.abs(x)>1e-5
    W[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    # integrate in ln k for numerical stability: sigma^2=1/(2pi^2) int dlnk k^3 P W^2
    y=k**3*P*W**2/(2*np.pi**2)
    s2=np.trapz(y,np.log(k))
    return float(np.sqrt(s2))

def read_pk(prefix):
    files=sorted(Path("output").glob(prefix+"*pk.dat"))
    if not files: raise RuntimeError(f"no pk.dat for {prefix}")
    # z_pk=0 only => single pk file
    return files[0]

def bao_pack():
    rows=[]
    for ln in (BAO/"desi_2024_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
        if not ln.strip() or ln.startswith("#"): continue
        z,v,q=ln.split(); rows.append((float(z),float(v),q))
    obs=np.array([r[1] for r in rows])
    cov=np.loadtxt(BAO/"desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
    return rows,obs,np.linalg.inv(cov)

bao_rows,bao_obs,bao_inv=bao_pack()

def metrics(prefix,H0):
    bg=table(f"output/{prefix}_00_background.dat")
    th=table(f"output/{prefix}_00_thermodynamics.dat")
    zdrag=find_drag_z(th)
    rd=interp(bg,"comov.snd.hrz.",zdrag)
    pred=[]
    pred_rows=[]
    for z,obs,q in bao_rows:
        DM=interp(bg,"comov. dist.",z)
        H=interp(bg,"H [1/Mpc]",z)
        DH=1/H
        DV=(z*DM*DM*DH)**(1/3)
        val={"DM_over_rs":DM/rd,"DH_over_rs":DH/rd,"DV_over_rs":DV/rd}[q]
        pred.append(val)
        pred_rows.append(dict(z=z,quantity=q,observed=obs,predicted=val,residual=val-obs))
    d=np.asarray(pred)-bao_obs
    bao=float(d@bao_inv@d)
    # z=0 is last row in these CLASS outputs; robustly interpolate.
    age=interp(bg,"proper time [Gyr]",0.0)
    Om=interp(bg,"Omega_m(z)",0.0)
    sig8=sigma8_from_pk(read_pk(prefix))
    S8=sig8*math.sqrt(Om/0.3)
    f0=interp(bg,"gr.fac. f",0.0) if "gr.fac. f" in bg.columns else float("nan")
    fs8=f0*sig8
    return dict(H0=H0,age_Gyr=age,Omega_m=Om,sigma8=sig8,S8=S8,f0=f0,fsigma8_0=fs8,
                z_drag=zdrag,r_drag_Mpc=rd,bao_chi2=bao,bao_rows=pred_rows,
                min_D=float(bg["kin (D)"].min()) if "kin (D)" in bg.columns else None,
                min_cs2=float(bg["c_s^2"].min()) if "c_s^2" in bg.columns else None,
                max_cs2=float(bg["c_s^2"].max()) if "c_s^2" in bg.columns else None)

s117=metrics("s117_validation",69.71482083084993)
local=metrics("local021_validation",68.56858744695782)
std=metrics("s117_gr_baseline",69.71482083084993)

# Differential BBN expansion screen at the earliest common tabulated epoch.
bs=table("output/s117_validation_00_background.dat")
bg=table("output/s117_gr_baseline_00_background.dat")
zmax=min(float(bs["z"].max()),float(bg["z"].max()))
Hs=interp(bs,"H [1/Mpc]",zmax)
Hg=interp(bg,"H [1/Mpc]",zmax)
speed=Hs/Hg
dNeff=(speed*speed-1.0)*(1.0+ALPHA_NU*NEFF)/ALPHA_NU
eta10=273.9*0.022083219194622913
etaD=eta10-6.0*(speed-1.0)
etaHe=eta10+100.0*(speed-1.0)
yD=2.60*(6.0/etaD)**1.6
Yp=0.2485+0.0016*(etaHe-6.0)
# same fit for GR baseline
etaD0=eta10
etaHe0=eta10
yD0=2.60*(6.0/etaD0)**1.6
Yp0=0.2485+0.0016*(etaHe0-6.0)

bbn=dict(z_screen=zmax,H_ratio_to_same_ordinary_GR=speed,equivalent_delta_Neff=dNeff,
         eta10=eta10,Yp_fit=Yp,Yp_GR_fit=Yp0,delta_Yp=Yp-Yp0,
         D_over_H_times_1e5_fit=yD,D_over_H_GR_times_1e5_fit=yD0,
         delta_D_over_H_times_1e5=yD-yD0,
         note="Differential analytic BBN screen using radiation-era H ratio; not a full nuclear-network likelihood.")

# Weak-lensing summary-statistic screens. These are NOT substitutes for
# a full modified-gravity cosmic-shear likelihood.
wl_refs={
  "DES_Y3":{"S8":0.776,"sigma":0.017},
  "KiDSLegacy_DESY3_plus_external":{"S8":0.814,"sigma":0.0115}
}
wl={}
for name,r in wl_refs.items():
    zs=(s117["S8"]-r["S8"])/r["sigma"]
    zl=(local["S8"]-r["S8"])/r["sigma"]
    wl[name]={
      "reference_S8":r["S8"],"reference_sigma":r["sigma"],
      "s117_z":zs,"local021_z":zl,
      "s117_gaussian_chi2":zs*zs,"local021_gaussian_chi2":zl*zl,
      "delta_chi2_s117_minus_local021":zs*zs-zl*zl
    }
wl["note"]="Summary-level S8 diagnostic only; full weak-lensing shear likelihood in SDMC requires model-consistent lensing kernels/nonlinear modeling."

age_screen={
  "oldest_GC_population_2026_Gyr":13.61,
  "oldest_GC_combined_sigma_Gyr":math.sqrt(0.25**2+0.23**2),
  "GC_inferred_universe_age_2026_Gyr":13.81,
  "GC_inferred_universe_combined_sigma_Gyr":math.sqrt(0.25**2+0.23**2),
  "s117_minus_GC_population_Gyr":s117["age_Gyr"]-13.61,
  "s117_minus_GC_inferred_universe_Gyr":s117["age_Gyr"]-13.81,
  "note":"Chronometer comparison screen; stellar-system ages have model/systematic uncertainties and are not a direct cosmological likelihood."
}

ledger={
 "s117":s117,
 "local021":local,
 "same_ordinary_GR":std,
 "comparisons":{
   "age_s117_minus_local021_Gyr":s117["age_Gyr"]-local["age_Gyr"],
   "bao_chi2_s117_minus_local021":s117["bao_chi2"]-local["bao_chi2"],
   "sigma8_s117_minus_local021":s117["sigma8"]-local["sigma8"],
   "S8_s117_minus_local021":s117["S8"]-local["S8"],
   "fsigma8_0_s117_minus_local021":s117["fsigma8_0"]-local["fsigma8_0"],
 },
 "bbn_screen":bbn,
 "weak_lensing_S8_screen":wl,
 "age_chronometer_screen":age_screen,
 "closure_ledger":{
   "delta_planck_s117":-20.194391439610172,
   "delta_desi_s117":-11.01153007793846,
   "delta_pd_s117":-31.20592151754863,
   "delta_pd_local021":-31.152232235955807,
   "fair_gap_pd":-0.053689281592824045,
   "fair_gap_pd_pantheonplus":-0.1400889900016864,
   "fair_gap_pd_union3":-0.04413438997858066,
   "fair_gap_pd_desy5":-0.12071862364814478
 }
}
(OUT/"s117_validation_ledger.json").write_text(json.dumps(ledger,indent=2))
rows=[]
for name,m in [("s117",s117),("local021",local),("same_ordinary_GR",std)]:
    rows.append({k:v for k,v in {"model":name,**m}.items() if k!="bao_rows"})
pd.DataFrame(rows).to_csv(OUT/"s117_validation_summary.csv",index=False)
pd.DataFrame(s117["bao_rows"]).to_csv(OUT/"s117_bao_rows.csv",index=False)
pd.DataFrame(local["bao_rows"]).to_csv(OUT/"local021_bao_rows.csv",index=False)
print("S117_VALIDATION_LEDGER",json.dumps(ledger,sort_keys=True),flush=True)
