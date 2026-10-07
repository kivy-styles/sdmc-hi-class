#!/usr/bin/env python3
from pathlib import Path
import importlib.util, json, numpy as np
from scipy.interpolate import RegularGridInterpolator

S=importlib.util.spec_from_file_location("hsc","sdmc/run_late300_hscy3_weyl.py")
h=importlib.util.module_from_spec(S); S.loader.exec_module(h)

REACT=Path("output/react_dense/react_late300/late300_react_target_spectrum.npz")
BEST=Path("output/hsc_best")
OUT=Path("output/hsc_react_cutscan"); OUT.mkdir(parents=True,exist_ok=True)
R=np.load(REACT)
ZR=np.asarray(R["z"],float); KR=np.asarray(R["k"],float)
KMIN=float(KR.min()); KMAX=float(KR.max()); ZMAX=2.5

def bi(arr):
    return RegularGridInterpolator((ZR,np.log(KR)),np.log(np.maximum(arr,1e-12)),
                                   bounds_error=False,fill_value=None)
BOOST={
 ("late300","us"):bi(R["Pnl_us"]/R["Ptarget"]),
 ("late300","ss"):bi(R["Pnl_ss"]/R["Ptarget"]),
 ("local021","gr"):bi(R["Pnl_local021"]/R["Ptarget_local021"]),
}
def eval_boost(I,z,kh):
    z=np.asarray(z,float); kh=np.asarray(kh,float)
    pts=np.column_stack([np.clip(z,0,ZMAX).ravel(),np.log(np.clip(kh,KMIN,KMAX).ravel())])
    q=np.exp(I(pts)).reshape(z.shape)
    w=np.clip((3.0-z)/0.5,0.0,1.0)
    return 1.0+w*(q-1.0)

class HSCReact(h.HSCNativeWeyl):
    def __init__(self,m,variant):
        super().__init__(m)
        self.variant=variant
        if variant=="linear": return
        Pl=self.P.copy(); Ql=self.Q.copy()
        zm=np.broadcast_to(self.z[None,:],Pl.shape)
        km=(self.ells[:,None]+0.5)/np.maximum(self.chi[None,:],1e-12)
        kh=km/m["h"]
        if m["name"]=="late300":
            tag="ss" if variant.startswith("ss_") else "us"
            I=BOOST[("late300",tag)]
        else:
            tag="gr"; I=BOOST[("local021","gr")]
        boost=eval_boost(I,zm,kh)
        P=Pl*boost
        H0=100*m["h"]/h.b.C_KMS
        Agr=(1.5*m["Om"]*H0**2*(1+self.z))[None,:]
        sig2=np.divide(Ql,Agr**2*Pl,out=np.ones_like(Ql),where=(Agr**2*Pl)>0)
        sig=np.sqrt(np.clip(sig2,0.2,2.0))
        if variant.endswith("_native"):
            sigeff=sig
        elif variant.endswith("_screen03"):
            x=kh/0.3
            sigeff=1.0+(sig-1.0)/(1.0+x*x)
        elif variant.endswith("_fullscreen"):
            sigeff=np.ones_like(sig)
        else:
            raise ValueError(variant)
        self.P=P
        self.Q=Agr**2*P*sigeff**2
        self.X=np.sqrt(np.maximum(self.Q*self.P,0.0))

def load_best(model):
    p=next((BEST/model).rglob(f"hsc_weyl_result_{model}.json"))
    d=json.loads(p.read_text())["results"][model]
    return np.r_[d["dz"],d["m"],d["A1"],d["alpha1"],d["psf_u"]]

def prior_chi2(p):
    return float((p[0]/h.DZ_SIG[0])**2+(p[1]/h.DZ_SIG[1])**2+
                 np.sum((np.asarray(p[4:8])/h.M_SIG)**2)+np.sum(np.asarray(p[10:14])**2))

def main():
    variants=["linear","us_native","ss_native","us_screen03","ss_screen03","us_fullscreen","ss_fullscreen"]
    cuts=[600,800,1000,1200,1600,1800]
    models={"late300":h.late,"local021":h.lcdm}
    raw_cov=h.COV/float(h.HARTLAP)
    out={"status":"HSC fixed-full-nuisance scale and corrected-ReACT diagnostic",
         "note":"diagnostic only: nuisance coordinates held at each model's completed full-band linear optimum",
         "variants":{}}
    for v in variants:
        vr={}
        cache={}
        for name,m in models.items():
            p=load_best(name)
            obj=HSCReact(m,v)
            th=obj.theory(p)
            cache[name]=(p,th)
        for ellmax in cuts:
            ix=np.array([i for i,pt in enumerate(h.POINTS) if pt["ell"]<=ellmax],int)
            n=len(ix); x=(1404.-1.)/(1404.-n-2.)
            C=raw_cov[np.ix_(ix,ix)]*x
            L=np.linalg.cholesky(C)
            rec={"ndata":n}
            for name in models:
                p,th=cache[name]
                y=np.linalg.solve(L,h.DATA[ix]-th[ix])
                data_chi=float(y@y); pr=prior_chi2(p)
                rec[name]={"data_chi2":data_chi,"prior_chi2":pr,"total_fixed":data_chi+pr}
            rec["delta_total_fixed"]=rec["late300"]["total_fixed"]-rec["local021"]["total_fixed"]
            vr[str(ellmax)]=rec
            print("HSC_REACT_CUT",v,ellmax,json.dumps(rec,sort_keys=True),flush=True)
        out["variants"][v]=vr
    (OUT/"summary.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("HSC_REACT_CUTSCAN_SUMMARY",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
