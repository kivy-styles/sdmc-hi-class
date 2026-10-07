#!/usr/bin/env python3
from pathlib import Path
import sys
import importlib.util, json, numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
import sacc

S=importlib.util.spec_from_file_location("base","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

ROOT=Path("cosmosis-standard-library")
SACCFILE=ROOT/"likelihood/hsc_cosmic_shear/hsc_y3_fourier_shear.sacc"
PSFFILE=ROOT/"likelihood/hsc_cosmic_shear/ppcorr_psf_all_ells_lmax_1800_catalog2.npz"
PSFTRANS=ROOT/"likelihood/hsc_cosmic_shear/psf_transform_matrix_lmax_1800_catalog2.npz"
OUT=Path("output/hsc_weyl"); OUT.mkdir(parents=True,exist_ok=True)

hL=.6971482083084993; OmL=(.022083219194622913+.12299536722293603)/hL**2
hC=.6856859; OmC=(.02240637+.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)
ZPIV=.62

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def load_hsc():
    s=sacc.Sacc.load_fits(str(SACCFILE))
    for dt in list(s.get_data_types()):
        if dt!="cl_ee": s.remove_selection(data_type=dt)
    for tr in list(s.get_tracer_combinations("cl_ee")):
        s.remove_selection("cl_ee",tr,ell__lt=300.)
        s.remove_selection("cl_ee",tr,ell__gt=1800.)
    data=np.asarray(s.get_mean(),float)
    cov=np.asarray(s.covariance.dense,float)
    r=1404.; p=len(data); hartlap_cov_factor=(r-1.)/(r-p-2.)
    cov=cov*hartlap_cov_factor
    tracers={}
    for name in ["wl_0","wl_1","wl_2","wl_3"]:
        t=s.get_tracer(name)
        z=np.asarray(t.z,float); nz=np.asarray(t.nz,float)
        nz=nz/np.trapezoid(nz,z)
        tracers[name]=(z,nz)
    psft=np.load(PSFFILE); psfx=np.load(PSFTRANS)
    psf_template=np.asarray(psft["arr_1"],float)
    psf_transform=np.asarray(psfx["arr_0"],float)
    psf_means=np.asarray(psfx["arr_1"],float)
    points=[]
    for d in s.data:
        b1,b2=d.tracers
        points.append(dict(b1=int(b1.split("_")[-1]),b2=int(b2.split("_")[-1]),
                           ell=float(d["ell"]),window=d["window"],window_ind=int(d["window_ind"])))
    return s,data,cov,tracers,psf_template,psf_transform,psf_means,points,hartlap_cov_factor

SACC,DATA,COV,TRACERS,PSFT,PSFX,PSFMEAN,POINTS,HARTLAP=load_hsc()

# HSC priors: dz1,dz2 Gaussian; dz3,dz4 top-hat only.
DZ_SIG=np.array([.024,.022])
M_SIG=.01

class HSCNativeWeyl:
    def __init__(self,m):
        self.m=m
        self.L=np.linalg.cholesky(COV)
        zmax=min(m["zq"].max(),m["zp"].max(),4.0)
        self.z=np.linspace(.01,zmax,320)
        self.H=np.interp(self.z,m["bg"]["z"],m["bg"]["H"])
        self.chi=np.interp(self.z,m["bg"]["z"],m["bg"]["chi"])
        # All HSC pairs in the retained sample share compatible window ell grids.
        vals=[]
        for p in POINTS:
            vals.extend(np.asarray(p["window"].values,float).tolist())
        self.ells=np.unique(np.asarray(vals,float))
        zm=np.broadcast_to(self.z[None,:],(len(self.ells),len(self.z)))
        km=(self.ells[:,None]+.5)/np.maximum(self.chi[None,:],1e-12)
        self.Q=b.eval_cube(m["Iq"],zm,km)
        self.P=b.eval_cube(m["Ip"],zm,km)
        ok=np.isfinite(self.Q)&np.isfinite(self.P)
        self.Q=np.where(ok,self.Q,0.0); self.P=np.where(ok,self.P,0.0)
        self.X=np.sqrt(np.maximum(self.Q*self.P,0.0))
        pref=np.exp(m["Ip"](np.column_stack([self.z,np.full_like(self.z,np.log(.05))])))
        p0=float(np.exp(m["Ip"]([[0.,np.log(.05)]])[0]))
        self.D=np.maximum(np.sqrt(pref/p0),1e-4)
        # index maps for each distinct SACC window grid
        self.window_maps={}
        for pt in POINTS:
            w=pt["window"]; key=id(w)
            if key not in self.window_maps:
                vals=np.asarray(w.values,float)
                self.window_maps[key]=np.searchsorted(self.ells,vals)

    def psf(self,u):
        p=np.linalg.inv(PSFX)@np.asarray(u,float) + PSFMEAN
        out=np.zeros(len(PSFT[0][0]))
        for i in range(4):
            for j in range(4):
                out += p[i]*p[j]*PSFT[i][j]
        return out

    def theory(self,p):
        dz=np.asarray(p[:4]); mm=np.asarray(p[4:8]); A=float(p[8]); alpha=float(p[9]); u=np.asarray(p[10:14])
        nz=[]
        for i,name in enumerate(["wl_0","wl_1","wl_2","wl_3"]):
            zz,nn=TRACERS[name]
            f=interp1d(zz,nn,bounds_error=False,fill_value=0.0)
            y=np.maximum(f(self.z-dz[i]),0.0)
            norm=np.trapezoid(y,self.z)
            if not np.isfinite(norm) or norm<=1e-12: return None
            nz.append(y/norm)
        nz=np.asarray(nz)
        g=[]
        for y in nz:
            g.append(tail(self.z,y)-self.chi*tail(self.z,y/np.maximum(self.chi,1e-12)))
        g=np.asarray(g)
        nchi=nz*self.H[None,:]
        F=-A*b.C1RHO*self.m["Om"]/self.D*((1+self.z)/(1+ZPIV))**alpha
        cls={}
        for i in range(4):
            for j in range(i,4):
                gg=np.trapezoid(g[i][None,:]*g[j][None,:]*self.Q,x=self.chi,axis=1)
                gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*F/np.maximum(self.chi,1e-12))[None,:]*self.X,x=self.chi,axis=1)
                ii=np.trapezoid((nchi[i]*nchi[j]*F**2/np.maximum(self.chi,1e-12)**2)[None,:]*self.P,x=self.chi,axis=1)
                cls[(i,j)]=(gg+gi+ii)*(1+mm[i])*(1+mm[j])
        psf=self.psf(u)
        tv=[]
        for pt in POINTS:
            pair=(min(pt["b1"],pt["b2"]),max(pt["b1"],pt["b2"]))
            w=pt["window"]; idx=self.window_maps[id(w)]
            cl=cls[pair][idx]
            wt=np.asarray(w.weight[:,pt["window_ind"]],float)
            val=float(wt@cl/wt.sum())
            wi=pt["window_ind"]
            if 0 <= wi < len(psf): val += float(psf[wi])
            tv.append(val)
        return np.asarray(tv)

    def objective(self,p):
        t=self.theory(p)
        if t is None or np.any(~np.isfinite(t)): return 1e100
        y=np.linalg.solve(self.L,DATA-t)
        chi=float(y@y)
        chi += float((p[0]/DZ_SIG[0])**2 + (p[1]/DZ_SIG[1])**2)
        chi += float(np.sum((np.asarray(p[4:8])/M_SIG)**2))
        chi += float(np.sum(np.asarray(p[10:14])**2))
        return chi

    def fit(self):
        x0=np.r_[0.,0.,.07,.15, np.zeros(4), .8,-4.7, np.zeros(4)]
        bounds=[(-1,1)]*4+[(-.1,.1)]*4+[(-6,6),(-6,6)]+[(-5,5)]*4
        starts=[x0,x0.copy(),x0.copy()]
        starts[1][8]=0.; starts[1][9]=0.
        starts[2][8]=-1.; starts[2][9]=0.
        best=None
        for x in starts:
            r=minimize(self.objective,x,method="L-BFGS-B",bounds=bounds,
                       options=dict(maxiter=260,ftol=2e-9,maxls=30))
            if best is None or r.fun<best.fun: best=r
        p=best.x
        return dict(chi2_total=float(best.fun),success=bool(best.success),message=str(best.message),
                    dz=p[:4].tolist(),m=p[4:8].tolist(),A1=float(p[8]),alpha1=float(p[9]),
                    psf_u=p[10:14].tolist(),nfev=int(best.nfev),nit=int(best.nit))

def main():
    requested=sys.argv[1] if len(sys.argv)>1 else "both"
    model_map={"late300":late,"local021":lcdm}
    if requested=="both":
        models=[late,lcdm]
    elif requested in model_map:
        models=[model_map[requested]]
    else:
        raise SystemExit(f"unknown model {requested}")

    out={"status":"HSC-Y3 native-Weyl bandpower NLA-z nuisance pilot",
         "ndata":len(DATA),"hartlap_cov_factor":float(HARTLAP),
         "ell_nominal_min":float(min(p["ell"] for p in POINTS)),
         "ell_nominal_max":float(max(p["ell"] for p in POINTS)),
         "priors":{"dz1_sigma":.024,"dz2_sigma":.022,"dz3_bounds":[-1,1],"dz4_bounds":[-1,1],
                   "m_sigma":.01,"A1_bounds":[-6,6],"alpha1_bounds":[-6,6],"psf_u_sigma":1.0},
         "missing":["full TATT A2/alpha2/bias_ta sector","validated nonlinear modified-gravity/baryonic prescription"]}
    out["results"]={}
    for m in models:
        print("HSC_WEYL_FIT_START",m["name"],flush=True)
        out["results"][m["name"]]=HSCNativeWeyl(m).fit()
        print("HSC_WEYL_FIT_DONE",m["name"],json.dumps(out["results"][m["name"]],sort_keys=True),flush=True)

    if len(models)==2:
        out["delta_chi2_late300_minus_local021"]=out["results"]["late300"]["chi2_total"]-out["results"]["local021"]["chi2_total"]
        outfile=OUT/"hsc_weyl_result.json"
    else:
        outfile=OUT/f"hsc_weyl_result_{models[0]['name']}.json"

    print("HSC_WEYL_RESULT",json.dumps(out,sort_keys=True),flush=True)
    outfile.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")

if __name__=="__main__": main()
