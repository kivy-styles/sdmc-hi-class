#!/usr/bin/env python3
import sys, json, importlib.util
from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from fastpt import FASTPT

S=importlib.util.spec_from_file_location("rx","sdmc/run_late300_hsc_react_cutscan.py")
rx=importlib.util.module_from_spec(S); S.loader.exec_module(rx)
h=rx.h
ZPIV=h.ZPIV

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

class HSCReactTATT(rx.HSCReact):
    """
    HSC corrected-ReACT theory with the same TATT construction used by the
    CosmoSIS standard-library TATT module.

    FAST-PT terms are generated from the model's linear matter spectrum.
    The NLA piece uses the corrected-ReACT nonlinear matter spectrum already
    loaded by HSCReact.  Matter-intrinsic power is converted to the native-Weyl
    shear-intrinsic power using X/P = sqrt(Q/P), matching the existing NLA
    convention exactly when A2=bias_ta=0.
    """
    def __init__(self,m,variant):
        super().__init__(m,variant)
        self._build_tatt_basis()

    def _build_tatt_basis(self):
        m=self.m
        h0=float(m["h"])
        # load_model stores k in 1/Mpc and P in Mpc^3. FAST-PT convention here
        # is h/Mpc and (Mpc/h)^3.
        ksrc=np.asarray(m["kp"],float)/h0
        Psrc=np.asarray(m["P"][0],float)*h0**3
        good=np.isfinite(ksrc)&np.isfinite(Psrc)&(ksrc>0)&(Psrc>0)
        ksrc=ksrc[good]; Psrc=Psrc[good]
        kmin=max(float(ksrc.min())*1.001,1e-4)
        kmax=min(float(ksrc.max())*.999,50.0)
        if not (kmax>10*kmin):
            raise RuntimeError(f"insufficient FAST-PT k range {kmin} {kmax}")
        k=np.geomspace(kmin,kmax,256)
        P0=np.exp(np.interp(np.log(k),np.log(ksrc),np.log(Psrc)))
        fpt=FASTPT(k,to_do=["IA_all"],low_extrap=-5,high_extrap=3,n_pad=len(k))
        tt=fpt.IA_tt(P0)
        ta=fpt.IA_ta(P0)
        mix=fpt.IA_mix(P0)
        base={
          "tt_EE":tt[0],"tt_BB":tt[1],
          "ta_dE1":ta[0],"ta_dE2":ta[1],"ta_EE":ta[2],"ta_BB":ta[3],
          "mix_A":mix[0],"mix_B":mix[1],"mix_D_EE":mix[2],"mix_D_BB":mix[3],
        }
        # Evaluate each z=0 FAST-PT term at the Limber k used by every HSC ell,
        # then evolve with D^4 as in the standard-library interface.
        zm=np.broadcast_to(self.z[None,:],self.P.shape)
        km=(self.ells[:,None]+0.5)/np.maximum(self.chi[None,:],1e-12) # 1/Mpc
        kh=km/h0
        lk=np.log(k)
        lq=np.log(np.clip(kh,kmin,kmax))
        D4=self.D[None,:]**4
        self.tatt={}
        for name,val in base.items():
            val=np.asarray(val,float)
            q=np.interp(lq.ravel(),lk,val).reshape(lq.shape)
            # (Mpc/h)^3 -> Mpc^3
            self.tatt[name]=q*D4/h0**3
        self.fastpt_k_bounds=[kmin,kmax]

    def theory_tatt(self,p):
        # p = dz4, m4, A1, alpha1, A2, alpha2, bias_ta, psf4
        dz=np.asarray(p[:4]); mm=np.asarray(p[4:8])
        A1=float(p[8]); alpha1=float(p[9]); A2=float(p[10])
        alpha2=float(p[11]); bias_ta=float(p[12]); u=np.asarray(p[13:17])
        nz=[]
        for i,name in enumerate(["wl_0","wl_1","wl_2","wl_3"]):
            zz,nn=h.TRACERS[name]
            f=interp1d(zz,nn,bounds_error=False,fill_value=0.0)
            y=np.maximum(f(self.z-dz[i]),0.0)
            norm=np.trapezoid(y,self.z)
            if not np.isfinite(norm) or norm<=1e-12: return None
            nz.append(y/norm)
        nz=np.asarray(nz)
        g=np.asarray([tail(self.z,y)-self.chi*tail(self.z,y/np.maximum(self.chi,1e-12)) for y in nz])
        nchi=nz*self.H[None,:]

        red=(1+self.z)/(1+ZPIV)
        C1=-A1*h.b.C1RHO*self.m["Om"]/self.D*red**alpha1
        Cdel=-(bias_ta*A1)*h.b.C1RHO*self.m["Om"]/self.D*red**alpha1
        C2=5.0*A2*h.b.C1RHO*self.m["Om"]/(self.D**2)*red**alpha2
        C1=C1[None,:]; Cdel=Cdel[None,:]; C2=C2[None,:]

        T=self.tatt
        mi=(C1*self.P
            + Cdel*(T["ta_dE1"]+T["ta_dE2"])
            + C2*(T["mix_A"]+T["mix_B"]))
        ii=(C1*C1*self.P
            + Cdel*Cdel*T["ta_EE"]
            + 2.0*C1*Cdel*(T["ta_dE1"]+T["ta_dE2"])
            + C2*C2*T["tt_EE"]
            + 2.0*C2*(C1*T["mix_A"]+C1*T["mix_B"]+Cdel*T["mix_D_EE"]))

        # Convert density-intrinsic to native-Weyl-intrinsic.  This reduces
        # algebraically to the existing NLA GI term when A2=bias_ta=0.
        wi=np.divide(self.X*mi,self.P,out=np.zeros_like(mi),where=self.P>0)

        cls={}
        for i in range(4):
            for j in range(i,4):
                gg=np.trapezoid(g[i][None,:]*g[j][None,:]*self.Q,x=self.chi,axis=1)
                gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])/np.maximum(self.chi,1e-12))[None,:]*wi,
                                x=self.chi,axis=1)
                iit=np.trapezoid((nchi[i]*nchi[j]/np.maximum(self.chi,1e-12)**2)[None,:]*ii,
                                 x=self.chi,axis=1)
                cls[(i,j)]=(gg+gi+iit)*(1+mm[i])*(1+mm[j])
        psf=self.psf(u)
        tv=[]
        for pt in h.POINTS:
            pair=(min(pt["b1"],pt["b2"]),max(pt["b1"],pt["b2"]))
            w=pt["window"]; idx=self.window_maps[id(w)]
            cl=cls[pair][idx]
            wt=np.asarray(w.weight[:,pt["window_ind"]],float)
            val=float(wt@cl/wt.sum())
            wi0=pt["window_ind"]
            if 0<=wi0<len(psf): val+=float(psf[wi0])
            tv.append(val)
        return np.asarray(tv)

    def objective_tatt(self,p):
        t=self.theory_tatt(p)
        if t is None or np.any(~np.isfinite(t)): return 1e100
        y=np.linalg.solve(self.L,h.DATA-t)
        chi=float(y@y)
        chi += float((p[0]/h.DZ_SIG[0])**2+(p[1]/h.DZ_SIG[1])**2)
        chi += float(np.sum((np.asarray(p[4:8])/h.M_SIG)**2))
        chi += float(np.sum(np.asarray(p[13:17])**2))
        return chi

    def validate_nla_limit(self):
        p_nla=np.r_[0.,0.,.07,.15,np.zeros(4),.8,-4.7,np.zeros(4)]
        p_tatt=np.r_[p_nla[:10],0.0,0.0,0.0,p_nla[10:14]]
        a=self.theory(p_nla); b=self.theory_tatt(p_tatt)
        d=np.asarray(a)-np.asarray(b)
        scale=max(float(np.max(np.abs(a))),1e-30)
        return {
          "max_abs":float(np.max(np.abs(d))),
          "max_rel_to_peak":float(np.max(np.abs(d))/scale),
          "rms":float(np.sqrt(np.mean(d*d))),
          "nla_objective":float(self.objective(p_nla)),
          "tatt_nla_limit_objective":float(self.objective_tatt(p_tatt)),
          "objective_abs_diff":float(abs(self.objective(p_nla)-self.objective_tatt(p_tatt))),
          "fastpt_k_hMpc_bounds":self.fastpt_k_bounds
        }

    def fit_tatt(self):
        # HSC standard-library starting coordinates for TATT, plus alternative
        # NLA-like and sign-flipped starts.
        base=np.r_[0.,0.,.07,.15,np.zeros(4), .8,-4.7,-.9,-4.0,.05, np.zeros(4)]
        starts=[base.copy() for _ in range(4)]
        starts[1][8:13]=[0.,0.,0.,0.,0.]
        starts[2][8:13]=[-1.,0.,1.,0.,.5]
        starts[3][8:13]=[1.,0.,-1.,0.,1.]
        bounds=[(-1,1)]*4+[(-.1,.1)]*4+[(-6,6),(-6,6),(-6,6),(-6,6),(0,2)]+[(-5,5)]*4
        best=None
        for i,x in enumerate(starts):
            r=minimize(self.objective_tatt,x,method="L-BFGS-B",bounds=bounds,
                       options=dict(maxiter=360,ftol=2e-9,maxls=35))
            print("HSC_TATT_START_DONE",i,float(r.fun),bool(r.success),str(r.message),flush=True)
            if best is None or r.fun<best.fun: best=r
        p=best.x
        return dict(chi2_total=float(best.fun),success=bool(best.success),message=str(best.message),
                    dz=p[:4].tolist(),m=p[4:8].tolist(),A1=float(p[8]),alpha1=float(p[9]),
                    A2=float(p[10]),alpha2=float(p[11]),bias_ta=float(p[12]),
                    psf_u=p[13:17].tolist(),nfev=int(best.nfev),nit=int(best.nit))

def main():
    mode=sys.argv[1]
    model_name=sys.argv[2]
    variant=sys.argv[3] if len(sys.argv)>3 else "us_native"
    m={"late300":h.late,"local021":h.lcdm}[model_name]
    obj=HSCReactTATT(m,variant)
    outdir=Path("output/hsc_tatt"); outdir.mkdir(parents=True,exist_ok=True)
    val=obj.validate_nla_limit()
    print("HSC_TATT_NLA_LIMIT",model_name,variant,json.dumps(val,sort_keys=True),flush=True)
    if val["max_rel_to_peak"]>5e-7 or val["objective_abs_diff"]>1e-5:
        raise RuntimeError("TATT NLA-limit validation failed")
    if mode=="validate":
        out={"status":"HSC corrected-ReACT TATT NLA-limit validation","model":model_name,
             "variant":variant,"validation":val}
    elif mode=="fit":
        fit=obj.fit_tatt()
        out={"status":"HSC-Y3 corrected-ReACT full TATT nuisance refit",
             "model":model_name,"variant":variant,"ndata":len(h.DATA),"fit":fit,
             "validation":val,
             "nuisance":{"TATT":["A1","alpha1","A2","alpha2","bias_ta"],
                         "photoz_dz":4,"shear_m":4,"psf_eigenmodes":4},
             "remaining_missing":["baryonic-feedback nuisance"],
             "qualification":"FAST-PT TATT terms use each model's linear matter spectrum; native-Weyl GI is obtained with sqrt(Q/P), exactly reproducing the existing NLA convention in the A2=bias_ta=0 limit."}
    else:
        raise SystemExit("mode must be validate or fit")
    p=outdir/f"{mode}_{model_name}_{variant}.json"
    p.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("HSC_TATT_RESULT",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
