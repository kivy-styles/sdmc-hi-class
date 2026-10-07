#!/usr/bin/env python3
from pathlib import Path
import importlib.util, json, re
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import minimize_scalar

S=importlib.util.spec_from_file_location("rx","sdmc/run_late300_hsc_react_cutscan.py")
rx=importlib.util.module_from_spec(S); S.loader.exec_module(rx)
h=rx.h

ROOT=Path("cosmosis-standard-library/structure/baryon_power_scaling/data")
BEST=Path("output/hsc_baryon_best")
OUT=Path("output/hsc_baryon_envelope"); OUT.mkdir(parents=True,exist_ok=True)

TABLES={
 "bahamas_t76":"logPkRatio_BAHAMAS_T7.6_WMAP9_L400N1024.dat",
 "bahamas_t80":"logPkRatio_BAHAMAS_T8.0_WMAP9_L400N1024.dat",
 "owls_agn":"logPkRatio_owls_AGN.dat",
 "tng100":"logPkRatio_TNG100.dat",
 "eagle":"logPkRatio_eagle.dat",
}

class Ratio:
    def __init__(self,path):
        with open(path) as f:
            cols=f.readline().strip().split(" ")[1:]
        cols=[c for c in cols if c]
        cols.reverse()
        z=[]
        for c in cols:
            q=re.findall(r"\d+",c)[0]
            z.append(int(q)/10.0**(len(q)-1))
        a=np.loadtxt(path,skiprows=1)
        logk=np.asarray(a[:,0],float)
        y=np.zeros((len(logk),len(z)))
        for i in range(len(z)):
            y[:,i]=a[:,-(i+1)]
        self.logk=logk; self.z=np.asarray(z,float)
        self.spl=RectBivariateSpline(self.logk,self.z,y)
    def logratio(self,kh,z):
        lk=np.clip(np.log10(np.maximum(kh,1e-30)),self.logk.min(),self.logk.max())
        zz=np.clip(z,self.z.min(),self.z.max())
        flat=self.spl.ev(lk.ravel(),zz.ravel()).reshape(lk.shape)
        return flat

class HSCBaryon(rx.HSCReact):
    def __init__(self,m,variant,ratio):
        super().__init__(m,variant)
        self.P0=self.P.copy(); self.Q0=self.Q.copy()
        zm=np.broadcast_to(self.z[None,:],self.P.shape)
        km=(self.ells[:,None]+0.5)/np.maximum(self.chi[None,:],1e-12)
        kh=km/m["h"]
        self.logR=ratio.logratio(kh,zm)
        self.set_amp(0.0)
    def set_amp(self,amp):
        fac=10.0**(float(amp)*self.logR)
        self.P=self.P0*fac
        self.Q=self.Q0*fac
        self.X=np.sqrt(np.maximum(self.P*self.Q,0.0))

def load_best(model):
    hits=list((BEST/model).rglob(f"{model}_us_native.json"))
    if len(hits)!=1:
        raise RuntimeError(f"best-fit json for {model}: {hits}")
    d=json.loads(hits[0].read_text())["fit"]
    p=np.r_[d["dz"],d["m"],d["A1"],d["alpha1"],d["psf_u"]]
    return p,float(d["chi2_total"])

def main():
    models={"late300":h.late,"local021":h.lcdm}
    best={k:load_best(k) for k in models}
    result={
      "status":"HSC corrected-ReACT baryonic factorization robustness envelope",
      "variant":"us_native",
      "assumption":"Apply the public CosmoSIS hydro-simulation P_bary/P_DMO ratio multiplicatively to corrected-ReACT P_delta and native-Weyl Q, holding the Weyl-to-matter response Sigma fixed.",
      "scope":"robustness diagnostic, not an official HSC HMCode/BACCO baryonic marginalization",
      "amplitude_definition":"log10 modulation is scaled by A_bary; A_bary=0 is DMO and A_bary=1 is the tabulated hydro scenario.",
      "amplitude_profile_bounds":[0.0,1.5],
      "tables":{}
    }

    # Baseline exact replay check.
    base={}
    for name,m in models.items():
        p,chi_saved=best[name]
        # Any ratio object at amp=0 must exactly reduce to the current us_native fit.
        rr=Ratio(ROOT/TABLES["owls_agn"])
        obj=HSCBaryon(m,"us_native",rr)
        chi=float(obj.objective(p))
        base[name]={"saved":chi_saved,"recomputed":chi,"abs_diff":abs(chi-chi_saved)}
        if abs(chi-chi_saved)>1e-6:
            raise RuntimeError(f"baseline replay mismatch {name}: {chi} vs {chi_saved}")
    result["baseline_replay"]=base

    for tag,fn in TABLES.items():
        rr=Ratio(ROOT/fn)
        rec={}
        for name,m in models.items():
            p,_=best[name]
            obj=HSCBaryon(m,"us_native",rr)
            obj.set_amp(1.0)
            chi1=float(obj.objective(p))
            def fun(a):
                obj.set_amp(float(a))
                return float(obj.objective(p))
            q=minimize_scalar(fun,bounds=(0.0,1.5),method="bounded",
                              options={"xatol":2e-4,"maxiter":80})
            rec[name]={
              "chi2_fixed_nuisance_Abary1":chi1,
              "Abary_profile_fixed_other_nuisance":float(q.x),
              "chi2_profile_Abary_fixed_other_nuisance":float(q.fun),
              "table_k_hMpc_range":[float(10**rr.logk.min()),float(10**rr.logk.max())],
              "table_z_range":[float(rr.z.min()),float(rr.z.max())],
            }
        rec["delta_chi2_Abary1"]=rec["late300"]["chi2_fixed_nuisance_Abary1"]-rec["local021"]["chi2_fixed_nuisance_Abary1"]
        rec["delta_chi2_profile_Abary"]=rec["late300"]["chi2_profile_Abary_fixed_other_nuisance"]-rec["local021"]["chi2_profile_Abary_fixed_other_nuisance"]
        result["tables"][tag]=rec
        print("HSC_BARYON_ENVELOPE_POINT",tag,json.dumps(rec,sort_keys=True),flush=True)

    vals=[v["delta_chi2_profile_Abary"] for v in result["tables"].values()]
    result["profiled_delta_chi2_range"]=[float(min(vals)),float(max(vals))]
    result["verdict"]="If the entire profiled delta-chi2 range remains negative, the current late300 HSC preference survives this public hydro-ratio baryonic factorization envelope at fixed non-baryonic nuisance coordinates. Full joint nuisance+baryon refits are still required before calling the HSC closure final."
    (OUT/"hsc_baryon_envelope.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("HSC_BARYON_ENVELOPE_SUMMARY",json.dumps(result,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
