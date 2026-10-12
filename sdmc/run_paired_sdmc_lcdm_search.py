#!/usr/bin/env python3
"""Matched exploratory CLASS evolution: SDMC Late300 versus flat LCDM local021.

Every candidate is evolved through the same CLASS/CMB/BAO/SN machinery. SDMC's
stage-1 structural background remains phenomenologically prescribed; passing
candidates must have their action reconstructed, frozen, and replayed separately.
This screen is NOT a new full-Plik, official DESI P(k), nested-evidence or MCMC fit.
"""
import argparse, csv, json, math, re, shutil, subprocess, textwrap, traceback
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

TCMB=2.7255; ORAD=4.17998772e-5; SIGMA_CAL=0.0025
OUT=Path("output/paired_free_evolution")
SNDATA=Path("sn_data"); BAODATA=Path("bao_data")
NAMES=("H0","omega_b","omega_cdm","A_s","n_s","tau_reio")
LCDM_BASE=dict(H0=68.56858744695782,omega_b=.022406369378007947,
    omega_cdm=.11824151052483357,A_s=2.05109266435559e-9,
    n_s=.9649593164240942,tau_reio=.044985382026527077)
SDMC_BASE=dict(H0=69.71482083084993,omega_b=.022083219194622913,
    omega_cdm=.12299536722293603,A_s=2.1047019076951746e-9,
    n_s=.9628416801318527,tau_reio=.053167105823755265,
    A_F=.02048146490100771,z_c=3.927876388467848,
    width=.33114133956842123,D0=.34231919445927034,
    kinetic_power=1.,D_floor=.05181542627513409,
    lambda_e=18.40625,z_t=16.19317622240633,
    A_late=.024822032167576252,B_late=.01264704344701022)
SCALE=np.array([1.15,.00018,.00115,4.5e-11,.0035,.006],float)
STRUCT_AXIS=(("A_F",.0013),("z_c",.18),("width",.027),("D_floor",.0030))
TAUA=.25; TAUB=1.5; DNT=.5

def flat_omega(p):
    h=p["H0"]/100.
    return 1.-(p["omega_b"]+p["omega_cdm"]+ORAD)/(h*h)

def design(family,budget=17):
    if budget<1:raise ValueError("budget must be positive")
    ref=dict(SDMC_BASE if family=="sdmc" else LCDM_BASE)
    yield "baseline",ref.copy()
    if budget==1:return
    # Matched ordinary-cosmology Sobol cloud, identical dimensionless offsets
    # and number of evaluations in both families.
    sob=qmc.Sobol(d=len(NAMES),scramble=True,seed=20261012)
    u=sob.random_base2(m=4)
    for i,v in enumerate(u):
        if i+1>=budget:break
        p=ref.copy()
        for key,delta in zip(NAMES,2*(v-.5)*SCALE):
            p[key]=float(p[key]+delta)
        yield f"ordinary_sobol{i:02d}",p
    # Matched budget with separately motivated axis exploration:
    # LCDM may vary its ordinary six parameters; SDMC can test a few
    # structural directions while keeping the six ordinary inputs frozen.
    left=budget-1-min(16,budget-1)
    if left<=0:return
    if family=="sdmc":
        dirs=[(key,sgn*step) for key,step in STRUCT_AXIS for sgn in (1,-1)]
    else:
        dirs=[(key,sgn*float(step)) for key,step in zip(NAMES,SCALE) for sgn in (1,-1)]
    for j in range(left):
        key,delta=dirs[j%len(dirs)]
        p=ref.copy(); p[key]=float(p[key]+delta)
        yield f"axis_{key}_{'plus' if delta>0 else 'minus'}_{j}",p

def build_ini(family,p,root):
    h=p["H0"]/100.
    common=f"""H0 = {p['H0']:.17g}
omega_b = {p['omega_b']:.17g}
omega_cdm = {p['omega_cdm']:.17g}
N_ncdm = 0
N_ur = 3.046
T_cmb = 2.7255
YHe = 0.2453
A_s = {p['A_s']:.17e}
n_s = {p['n_s']:.17g}
tau_reio = {p['tau_reio']:.17g}
"""
    if family=="lcdm":
        spec=f"""Omega_k = 0
Omega_Lambda = {flat_omega(p):.17g}
"""
    else:
        spec=f"""Omega_Lambda = 0
Omega_fld = 3.1443554e-8
fluid_equation_of_state = SDMC_TRACKER
cs2_fld = 0.003
use_ppf = no
Omega_smg = -1
gravity_model = sdmc_v3_independent_kinetic
parameters_smg = {p['A_F']:.17g}, {p['z_c']:.17g}, {p['width']:.17g}, {p['D0']:.17g}, {p['kinetic_power']:.17g}, {p['D_floor']:.17g}
expansion_model = sdmc_full
expansion_smg = {flat_omega(p):.17g}, {p['lambda_e']:.17g}, {p['z_t']:.17g}, {DNT}, {p['A_late']:.17g}, {TAUA}, {p['B_late']:.17g}, {TAUB}
pert_initial_conditions_smg = zero
method_qs_smg = fully_dynamic
a_ini_over_a_today_default = 1.e-8
a_ini_test_qs_smg = 1.e-8
pert_ic_ini_z_ref_smg = 1.e7
a_min_stability_test_smg = 1.e-8
output_background_smg = 3
"""
    commonout=f"""modes = s
output = tCl,pCl,lCl,mPk
lensing = yes
l_max_scalars = 3000
P_k_max_h/Mpc = 1.0
z_pk = 0
write background = yes
write thermodynamics = yes
root = {root}
format = class
input_verbose = 0
background_verbose = 0
thermodynamics_verbose = 0
perturbations_verbose = 0
spectra_verbose = 0
lensing_verbose = 0
output_verbose = 0
"""
    return common+spec+commonout

def read_table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    marks=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    names=[hdr[m.end():(marks[i+1].start() if i+1<len(marks) else len(hdr))].strip()
            for i,m in enumerate(marks)]
    a=np.loadtxt(path)
    return pd.DataFrame(a,columns=names).sort_values("z")

def cmb_files(path):
    a=np.loadtxt(path)
    ell=a[:,0].astype(int)
    n=int(ell.max())+1
    v=(TCMB*1e6)**2
    tt=np.zeros(n);te=np.zeros(n);ee=np.zeros(n);pp=np.zeros(n)
    tt[ell]=a[:,1]*v;ee[ell]=a[:,2]*v;te[ell]=a[:,3]*v
    pp[ell]=a[:,5]*ell*(ell+1)
    return (np.column_stack((ell,tt[ell],te[ell],ee[ell])),
            dict(tt=tt,te=te,ee=ee,pp=pp))

class Likelihoods:
    def __init__(self,extra_sn=False):
        self.high=TTTEEE_lite_native(packages_path="planck_packages")
        self.lowT=TT(packages_path="planck_packages")
        self.lowE=EE(packages_path="planck_packages")
        self.lens=LensingNative(packages_path="planck_packages")
        self.sn={}
        self.prepare_sn("pantheonplus")
        if extra_sn:
            self.prepare_sn("union3");self.prepare_sn("desy5")
        self.prepare_bao()

    def planck(self,path):
        spec,dls=cmb_files(path)
        def f(A):
            lh=float(self.high.chi_squared(spec,A_planck=A))
            lt=float(-2*self.lowT.log_likelihood(dls["tt"],calib=A))
            le=float(-2*self.lowE.log_likelihood(dls["ee"],calib=A))
            lp={self.lens.calibration_param:A} if getattr(self.lens,"calibration_param",None) else {}
            ll=float(-2*self.lens.log_likelihood(dls,**lp))
            pr=((A-1)/SIGMA_CAL)**2
            return (lh+lt+le+ll+pr)
        opt=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
        return float(opt.fun),float(opt.x)

    @staticmethod
    def readcov(path,n):
        raw=np.asarray(np.loadtxt(path),float).ravel()
        if len(raw)==n*n+1 and int(round(raw[0]))==n:
            raw=raw[1:]
        if len(raw)!=n*n:raise RuntimeError("Bad covariance size for "+str(path))
        return raw.reshape(n,n)

    def prepare_sn(self,kind):
        if kind=="pantheonplus":
            d=pd.read_csv(SNDATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",engine="python")
            mag=d["m_b_corr"].to_numpy(float);z=d["zHD"].to_numpy(float);zh=d["zHEL"].to_numpy(float)
            C=self.readcov(SNDATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(mag))
            m=z>.01; mag,z,zh,C=mag[m],z[m],zh[m],C[np.ix_(m,m)]
        elif kind=="union3":
            p=SNDATA/"Union3/lcparam_full.txt";cols=p.read_text().splitlines()[0].removeprefix("#").split()
            d=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
            cm={c.lower():c for c in d.columns}
            mag=d[cm["mb"]].to_numpy(float);z=d[cm["zcmb"]].to_numpy(float);zh=z
            C=self.readcov(SNDATA/"Union3/mag_covmat.txt",len(mag))
        else:
            d=pd.read_csv(SNDATA/"DESY5/DES-SN5YR_HD.csv");cm={c.lower():c for c in d.columns}
            mag=d[cm["mu"]].to_numpy(float);z=d[cm["zhd"]].to_numpy(float)
            zh=d[cm["zhel"]].to_numpy(float);err=d[cm["muerr_final"]].to_numpy(float)
            C=self.readcov(SNDATA/"DESY5/covsys_000.txt",len(mag))+np.diag(err*err)
        cf=cho_factor(C,lower=True,check_finite=False)
        one=np.ones(len(mag))
        u=cho_solve(cf,one,check_finite=False)
        self.sn[kind]=(mag,z,zh,cf,one,float(one@u))

    def score_sn(self,bg):
        out={}
        for kind,(mag,z,zh,cf,one,den) in self.sn.items():
            dm=np.interp(z,bg.z,bg["comov. dist."])
            DA=dm/(1+z)
            mu=5*np.log10((1+zh)*(1+z)*DA)
            resid=mag-mu
            y=cho_solve(cf,resid,check_finite=False)
            out[kind]=float(resid@y-(one@y)**2/den)
        return out

    def prepare_bao(self):
        rows=[]
        for line in (BAODATA/"desi_2024_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                z,val,which=line.split()
                rows.append((float(z),float(val),which))
        C=np.loadtxt(BAODATA/"desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
        self.baorows=rows
        self.baoinv=np.linalg.inv(C)
        self.baoobs=np.asarray([r[1] for r in rows])

    @staticmethod
    def acoustics(bg,th):
        if "g [Mpc^-1]" not in th.columns or "kappa_b" not in th.columns:
            raise RuntimeError("Missing thermodynamics columns for acoustic inference")
        zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
        zz=th.z.to_numpy(float);kb=th.kappa_b.to_numpy(float)
        crossing=[]
        for j in range(len(zz)-1):
            if kb[j]!=kb[j+1] and (kb[j]-1.)*(kb[j+1]-1.)<=0:
                q=(1.-kb[j])/(kb[j+1]-kb[j])
                crossing.append(float(zz[j]+q*(zz[j+1]-zz[j])))
        if not crossing:raise RuntimeError("No baryon drag crossing")
        zd=crossing[0]
        DM=float(np.interp(zstar,bg.z,bg["comov. dist."]))
        rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
        rd=float(np.interp(zd,bg.z,bg["comov.snd.hrz."]))
        return {"z_star":zstar,"z_drag":zd,"theta_star":rs/DM,
                "theta100":100*rs/DM,"ell_A":np.pi*DM/rs,
                "rs_star_Mpc":rs,"rd_Mpc":rd,"DM_star_Mpc":DM}

    def bao(self,bg,rd):
        pred=[]
        for z,_,what in self.baorows:
            DM=float(np.interp(z,bg.z,bg["comov. dist."]))
            H=float(np.interp(z,bg.z,bg["H [1/Mpc]"]))
            DH=1/H
            DV=(z*DM*DM*DH)**(1./3)
            pred.append({"DM_over_rs":DM/rd,"DH_over_rs":DH/rd,"DV_over_rs":DV/rd}[what])
        d=np.asarray(pred)-self.baoobs
        return float(d@self.baoinv@d)

def evaluate(family,label,p,like,keep_selected=True):
    rootfolder=OUT/family
    rootfolder.mkdir(parents=True,exist_ok=True)
    prefix=str(rootfolder/(label+"_"))
    ini=rootfolder/(label+".ini")
    ini.write_text(build_ini(family,p,prefix))
    rec={"family":family,"label":label,**p,"status":"FAIL","Omega_de0":flat_omega(p)}
    try:
        cp=subprocess.run(["./class",str(ini)],capture_output=True,text=True,timeout=360)
        rec["returncode"]=cp.returncode
        if cp.returncode!=0:raise RuntimeError(cp.stdout[-800:]+cp.stderr[-800:])
        bg=read_table(prefix+"00_background.dat")
        th=read_table(prefix+"00_thermodynamics.dat")
        if family=="sdmc":
            D=bg["kin (D)"].to_numpy(float);cs=bg["c_s^2"].to_numpy(float)
            noslip=np.asarray(bg["braiding_smg"]+2*bg["M2_running_smg"])
            rec.update(D_min=float(D.min()),cs2_min=float(cs.min()),
                cs2_max=float(cs.max()),noslip_max=float(np.max(np.abs(noslip))))
            if rec["D_min"]<=0 or rec["cs2_min"]<=0 or rec["cs2_max"]>1.000001:
                raise RuntimeError("Full-background scalar stability/subluminality gate failed")
        planck,Acal=like.planck(prefix+"00_cl_lensed.dat")
        sn=like.score_sn(bg)
        ac=like.acoustics(bg,th)
        bao=like.bao(bg,ac["rd_Mpc"])
        rec.update(status="OK",chi2_planck=planck,A_planck=Acal,
                   chi2_bao=bao,**{f"chi2_{k}":v for k,v in sn.items()},**ac,
                   chi2_shoes_proxy=((p["H0"]-73.49)/1.048)**2,
                   theta100_residual_vs_Planck=ac["theta100"]-1.0411)
        rec["chi2_primary"]=planck+bao+sn["pantheonplus"]
        for key in ("union3","desy5"):
            if key in sn:rec[f"chi2_planck_bao_{key}"]=planck+bao+sn[key]
        rec["old_target_screen_not_unique_covariant_action"]=family=="sdmc"
        if keep_selected and family=="sdmc":
            # Retain all SDMC screened backgrounds only until promotion succeeds.
            covdir=OUT/"sdmc_screen_cache"
            covdir.mkdir(exist_ok=True)
            shutil.copy2(ini,covdir/(label+".ini"))
            shutil.copy2(prefix+"00_background.dat",covdir/(label+"_00_background.dat"))
        print("PAIRED_CANDIDATE",json.dumps(rec,sort_keys=True),flush=True)
    except Exception as e:
        rec["error"]=repr(e)
        print("PAIRED_CANDIDATE",json.dumps(rec,sort_keys=True),flush=True)
    finally:
        for ext in ("00_background.dat","00_thermodynamics.dat","00_cl.dat",
                    "00_cl_lensed.dat","00_pk.dat","00_pk_cb.dat","00_perturbations.dat"):
            fp=Path(prefix+ext)
            if fp.exists():fp.unlink()
    return rec

def run(budget,extra_sn):
    OUT.mkdir(parents=True,exist_ok=True)
    like=Likelihoods(extra_sn=extra_sn)
    rows=[]
    summaryfile=OUT/"summary.json"
    for family in ("lcdm","sdmc"):
        for label,p in design(family,budget):
            r=evaluate(family,label,p,like)
            rows.append(r)
            pd.DataFrame(rows).to_csv(OUT/"candidate_table.csv",index=False)
            success=[x for x in rows if x["status"]=="OK"]
            by={k:sorted((x for x in success if x["family"]==k),key=lambda v:v["chi2_primary"])
                for k in ("lcdm","sdmc")}
            summary={"status":"screen_in_progress","budget_per_family":budget,
              "screen_likelihoods":["Planck2018 Plik-lite lowT lowE lensing","DESI_DR1_BAO",
                                    "PantheonPlus_profiled"]+
                                    (["Union3_profiled","DESY5_profiled"] if extra_sn else []),
              "evaluation_counts":{k:sum(v["family"]==k for v in rows) for k in ("lcdm","sdmc")},
              "success_counts":{k:len(by[k]) for k in ("lcdm","sdmc")},
              "best":{k:by[k][0] if by[k] else None for k in ("lcdm","sdmc")},
              "verdict":"Screen only; no official DESI full-shape, full Plik nuisance, or independent SDMC action selection."
            }
            summaryfile.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    # Preserve promotion target from selected SDMC target simulation to retest by
    # genuinely evolved frozen reconstructed covariant coefficients.
    top=summary["best"]["sdmc"]
    if top is not None:
        label=top["label"]
        cache=OUT/"sdmc_screen_cache"
        Path("output/linear_cov_target.ini").write_text((cache/(label+".ini")).read_text()
            .replace(f"root = {OUT/'sdmc'/(label+'_')}", "root = output/linear_cov_target_"))
        shutil.copy2(cache/(label+"_00_background.dat"),
                     "output/linear_cov_target_00_background.dat")
        (OUT/"promotion_target.json").write_text(json.dumps(top,indent=2,sort_keys=True))
        summary["sdmc_covariant_promotion_target"]=label
    summary["status"]="screen_complete"
    summaryfile.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    return summary

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--budget",type=int,default=17)
    parser.add_argument("--extra-sn",action="store_true")
    args=parser.parse_args()
    out=run(args.budget,args.extra_sn)
    print("PAIRED_EVOLUTION_SUMMARY",json.dumps({
        "status":out["status"],"evaluation_counts":out["evaluation_counts"],
        "success_counts":out["success_counts"],
        "top":{k:(None if v is None else {key:v.get(key) for key in
                ("label","H0","chi2_primary","chi2_planck","chi2_bao","theta100",
                 "chi2_pantheonplus","chi2_shoes_proxy")})
               for k,v in out["best"].items()}
    },sort_keys=True),flush=True)
if __name__=="__main__":main()
