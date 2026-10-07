#!/usr/bin/env python3
from pathlib import Path
import sys, importlib.util, json, numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from scipy.special import jv

S=importlib.util.spec_from_file_location("base","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
REACT=Path("output/react_dense/react_late300/late300_react_target_spectrum.npz")
OUT=Path("output/desy3_react_nuisance"); OUT.mkdir(parents=True,exist_ok=True)

DZ_SIG=np.array([0.018,0.015,0.011,0.017])
M_MU=np.array([-0.0063,-0.0198,-0.0241,-0.0369])
M_SIG=np.array([0.0091,0.0078,0.0076,0.0076])
ZPIV=0.62
ELL=b.ELL
ZREACT_MAX=2.5

hL=0.6971482083084993; OmL=(0.022083219194622913+0.12299536722293603)/hL**2
hC=0.6856859; OmC=(0.02240637+0.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)

R=np.load(REACT)
ZR=np.asarray(R["z"],float); KR=np.asarray(R["k"],float)
KMIN=float(KR.min()); KMAX=float(KR.max())

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
    pts=np.column_stack([
        np.clip(z,0,ZREACT_MAX).ravel(),
        np.log(np.clip(kh,KMIN,KMAX).ravel())
    ])
    q=np.exp(I(pts)).reshape(z.shape)
    w=np.clip((3.0-z)/0.5,0.0,1.0)
    return 1.0+w*(q-1.0)

def subset(des,xp,xm):
    use=[i for i,r in enumerate(des["rows"]) if r["ang"] >= (xp if r["typ"]=="xip" else xm)]
    ix=np.array(use,int)
    return dict(rows=[des["rows"][i] for i in use],data=des["data"][ix],
                cov=des["cov"][np.ix_(ix,ix)],zsrc=des["zsrc"],nz=des["nz"])

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

class ReactWeylModel:
    def __init__(self,m,des,variant):
        self.m=m; self.des=des; self.variant=variant
        zmax=min(m["zq"].max(),m["zp"].max(),float(np.max(des["zsrc"])))
        self.z=np.linspace(0.01,zmax,300)
        self.H=np.interp(self.z,m["bg"]["z"],m["bg"]["H"])
        self.chi=np.interp(self.z,m["bg"]["z"],m["bg"]["chi"])
        zm=np.broadcast_to(self.z[None,:],(len(ELL),len(self.z)))
        km=(ELL[:,None]+0.5)/np.maximum(self.chi[None,:],1e-12)
        kh=km/m["h"]
        Q=b.eval_cube(m["Iq"],zm,km)
        Pl=b.eval_cube(m["Ip"],zm,km)
        ok=np.isfinite(Q)&np.isfinite(Pl)&(Pl>0)
        Q=np.where(ok,Q,0.0); Pl=np.where(ok,Pl,0.0)

        if variant=="linear":
            self.P=Pl; self.Q=Q
            self.matter_tag="linear"; self.weyl_mode="linear"
        else:
            if m["name"]=="late300":
                tag="ss" if variant.startswith("ss_") else "us"
                I=BOOST[("late300",tag)]
            else:
                tag="gr"; I=BOOST[("local021","gr")]
            boost=eval_boost(I,zm,kh)
            self.P=Pl*boost
            H0=100*m["h"]/b.C_KMS
            Agr=(1.5*m["Om"]*H0**2*(1+self.z))[None,:]
            sig2=np.divide(Q,Agr**2*Pl,out=np.ones_like(Q),where=(Agr**2*Pl)>0)
            sig=np.sqrt(np.clip(sig2,0.2,2.0))
            if variant.endswith("_native"):
                sigeff=sig; mode="native_linear_Sigma"
            elif variant.endswith("_screen03"):
                x=kh/0.3
                sigeff=1.0+(sig-1.0)/(1.0+x*x)
                mode="No-Slip transition k_s=0.3 h/Mpc"
            elif variant.endswith("_fullscreen"):
                sigeff=np.ones_like(sig); mode="full nonlinear Weyl screen to GR Sigma=1"
            else:
                raise ValueError(variant)
            self.Q=Agr**2*self.P*sigeff**2
            self.Q=np.where(ok,self.Q,0.0)
            self.matter_tag=tag; self.weyl_mode=mode

        self.X=np.sqrt(np.maximum(self.Q*self.P,0.0))
        pref=np.exp(m["Ip"](np.column_stack([self.z,np.full_like(self.z,np.log(0.05))])))
        p0=float(np.exp(m["Ip"]([[0.0,np.log(0.05)]])[0]))
        self.D=np.maximum(np.sqrt(pref/p0),1e-4)
        self.L=np.linalg.cholesky(des["cov"])
        self.bw=[]
        for r in des["rows"]:
            th=np.deg2rad(r["ang"]/60.0)
            self.bw.append(ELL/(2*np.pi)*jv(0 if r["typ"]=="xip" else 4,ELL*th))
        self.bw=np.asarray(self.bw)

    def theory(self,p):
        dz=np.asarray(p[:4]); mm=np.asarray(p[4:8]); A=float(p[8]); alpha=float(p[9])
        nz=[]
        for i in range(4):
            f=interp1d(self.des["zsrc"],self.des["nz"][i],bounds_error=False,fill_value=0.0)
            y=np.maximum(f(self.z-dz[i]),0.0)
            norm=np.trapezoid(y,self.z)
            if not np.isfinite(norm) or norm<=1e-12: return None
            nz.append(y/norm)
        nz=np.asarray(nz)
        g=[]
        for y in nz:
            q0=tail(self.z,y); q1=tail(self.z,y/np.maximum(self.chi,1e-12))
            g.append(q0-self.chi*q1)
        g=np.asarray(g)
        nchi=nz*self.H[None,:]
        F=-A*b.C1RHO*self.m["Om"]/self.D*((1+self.z)/(1+ZPIV))**alpha
        pair={}
        for i in range(4):
            for j in range(i,4):
                gg=np.trapezoid(g[i][None,:]*g[j][None,:]*self.Q,x=self.chi,axis=1)
                gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*F/np.maximum(self.chi,1e-12))[None,:]*self.X,x=self.chi,axis=1)
                ii=np.trapezoid((nchi[i]*nchi[j]*F**2/np.maximum(self.chi,1e-12)**2)[None,:]*self.P,x=self.chi,axis=1)
                pair[(i+1,j+1)]=(gg+gi+ii)*(1+mm[i])*(1+mm[j])
        tv=[]
        for q,r in enumerate(self.des["rows"]):
            cl=pair[(min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))]
            tv.append(np.trapezoid(self.bw[q]*cl,x=ELL))
        return np.asarray(tv)

    def objective(self,p):
        th=self.theory(p)
        if th is None or np.any(~np.isfinite(th)): return 1e100
        y=np.linalg.solve(self.L,self.des["data"]-th)
        chi=float(y@y)
        chi+=float(np.sum((np.asarray(p[:4])/DZ_SIG)**2))
        chi+=float(np.sum((np.asarray(p[4:8])-M_MU)**2/M_SIG**2))
        return chi

    def fit(self):
        x0=np.r_[np.zeros(4),M_MU,0.2,0.0]
        bounds=[(-0.08,0.08)]*4+[(-0.1,0.1)]*4+[(-5,5),(-5,5)]
        starts=[x0,x0.copy(),x0.copy()]
        starts[1][8]=1.0; starts[2][8]=-1.0
        best=None
        for x in starts:
            r=minimize(self.objective,x,method="L-BFGS-B",bounds=bounds,
                       options=dict(maxiter=220,ftol=1e-9,maxls=30))
            if best is None or r.fun<best.fun: best=r
        p=best.x
        return dict(chi2_total=float(best.fun),success=bool(best.success),message=str(best.message),
                    dz=p[:4].tolist(),m=p[4:8].tolist(),A1=float(p[8]),alpha1=float(p[9]),
                    nfev=int(best.nfev),nit=int(best.nit),
                    matter_tag=self.matter_tag,weyl_mode=self.weyl_mode)

def main():
    variant=sys.argv[1] if len(sys.argv)>1 else "linear"
    allowed={"linear","us_native","ss_native","us_screen03","ss_screen03","us_fullscreen","ss_fullscreen"}
    if variant not in allowed: raise SystemExit(f"unknown variant {variant}")
    full=b.load_des(DATA,CUTS,False)
    cases={"official":subset(full,0,0),"c40_150":subset(full,40,150),"c100_250":subset(full,100,250)}
    out={"status":"DES-Y3 ReACT nonlinear native-Weyl nuisance-profile audit",
         "variant":variant,
         "priors":{"dz_sigma":DZ_SIG.tolist(),"m_mean":M_MU.tolist(),"m_sigma":M_SIG.tolist(),
                   "A1":[-5,5],"alpha1":[-5,5],"z_piv":ZPIV},
         "nonlinear_limits":{"ReACT_zmax":ZREACT_MAX,"kmax_h_Mpc":KMAX,
                             "high_z":"boost tapered to unity z=2.5..3",
                             "high_k":"boost frozen at validated kmax"},
         "missing":["baryonic-feedback nuisance","TATT A2/eta2 sector"]}
    out["results"]={}
    for label,des in cases.items():
        rec={"ndata":len(des["data"])}
        for m in (late,lcdm):
            print("REACT_NUISANCE_FIT_START",variant,label,m["name"],flush=True)
            fit=ReactWeylModel(m,des,variant).fit()
            rec[m["name"]]=fit
            print("REACT_NUISANCE_FIT_DONE",variant,label,m["name"],json.dumps(fit,sort_keys=True),flush=True)
        rec["delta_chi2_late300_minus_local021"]=rec["late300"]["chi2_total"]-rec["local021"]["chi2_total"]
        out["results"][label]=rec
        print("REACT_NUISANCE_RESULT",variant,label,json.dumps(rec,sort_keys=True),flush=True)
    path=OUT/f"desy3_react_nuisance_{variant}.json"
    path.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("DESY3_REACT_NUISANCE_SUMMARY",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
