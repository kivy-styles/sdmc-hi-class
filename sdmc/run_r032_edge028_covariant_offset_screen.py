#!/usr/bin/env python3
from pathlib import Path
import json, math, re, subprocess
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/edge028_covariant_offset_screen"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.79120830851635
OB=0.022083219194622913
OC=0.12299536722293603
AS=2.1218418557926507e-9
NS=0.9625227132590487
TAU=0.055202901571989066
OFFSETS=[-0.05,-0.02,-0.01,-0.005,-0.002,-0.001,0.0,0.001,0.002,0.005,0.01,0.02,0.05]

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def derived(bg,th):
    b=bg.sort_values("z")
    imax=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,b["z"],b["comov. dist."]))
    rs=float(np.interp(zstar,b["z"],b["comov.snd.hrz."]))
    return dict(z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=float(np.pi*DM/rs))

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    harr,dls=load_cls(path)
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
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

target=table("output/cov_exact_target_00_background.dat")
Nt=-np.log1p(target["z"].to_numpy()); order=np.argsort(Nt); Nt=Nt[order]
Ht=target["H [1/Mpc]"].to_numpy()[order]
Ft=target["M*^2_smg"].to_numpy()[order]

def ini(root,p):
    return f"""H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={AS:.17e}
n_s={NS:.17g}
tau_reio={TAU:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_covariant_exact
parameters_smg={p:.17g}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
write background=yes
write thermodynamics=yes
root={root}
format=class
input_verbose=0
background_verbose=0
thermodynamics_verbose=0
perturbations_verbose=0
spectra_verbose=0
lensing_verbose=0
output_verbose=0
"""

rows=[]
for i,p in enumerate(OFFSETS):
    tag=f"off_{i:02d}_{p:+.3g}".replace("+","p").replace("-","m").replace(".","p")
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,p))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    rec={"id":tag,"offset":float(p),"returncode":int(cp.returncode),"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); thp=Path(root+"00_thermodynamics.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and thp.exists():
        bg=table(bgp); th=table(thp)
        z=bg["z"].to_numpy(); N=-np.log1p(z)
        Href=np.interp(N,Nt,Ht); Fref=np.interp(N,Nt,Ft)
        m100=z<=100.; m1e4=z<=1e4
        hrel=np.abs(bg["H [1/Mpc]"].to_numpy()/Href-1.)
        frel=np.abs(bg["M*^2_smg"].to_numpy()/Fref-1.)
        rec.update(derived(bg,th))
        rec.update(
            max_abs_dH_H_z100=float(hrel[m100].max()),
            rms_dH_H_z100=float(np.sqrt(np.mean(hrel[m100]**2))),
            max_abs_dH_H_z1e4=float(hrel[m1e4].max()),
            max_abs_dF_F_z100=float(frel[m100].max()),
            min_D=float(bg["kin (D)"].min()),
            min_cs2=float(bg["c_s^2"].min()),
            max_cs2=float(bg["c_s^2"].max()),
            max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
        )
        stable=rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.
        rec["stable_subluminal"]=bool(stable)
        if stable and clp.exists():
            rec.update(pscore(clp)); rec["status"]="OK"
        else:
            rec["status"]="UNSTABLE" if not stable else "NO_CLS"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-800:].replace("\n"," | ")
    print("OFFSET_POINT",json.dumps(rec,sort_keys=True),flush=True)
    rows.append(rec)

df=pd.DataFrame(rows); df.to_csv(OUT/"edge028_covariant_offset_screen.csv",index=False)
ok=df[df.status=="OK"].copy()
summary={
 "n_total":len(df),"n_ok":len(ok),
 "best_background": None if len(ok)==0 else ok.nsmallest(1,"rms_dH_H_z100").iloc[0].to_dict(),
 "best_planck": None if len(ok)==0 else ok.nsmallest(1,"chi2_planck").iloc[0].to_dict(),
}
(OUT/"edge028_covariant_offset_screen_summary.json").write_text(json.dumps(summary,indent=2,default=float))
print("OFFSET_SUMMARY",json.dumps(summary,sort_keys=True,default=float),flush=True)
