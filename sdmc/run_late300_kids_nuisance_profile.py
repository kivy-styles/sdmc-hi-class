#!/usr/bin/env python3
"""
KiDS-Legacy native-Weyl nuisance-profile audit for late300 vs local021.

Lensing GG always uses native linear hi_class Weyl power Q_W=k^4 P_Weyl.
Nuisances reproduce the KiDS-Legacy correlated prior construction:
  - 6 correlated photo-z shifts;
  - mass-dependent IA: A, beta, 6 correlated log10<M*> terms;
  - fixed f_r and log10 M_piv from the public KiDS-Legacy config.
The COSEBI transform uses the public precomputed W_n(l) tables for 2'-300',
n=1..6, validated against the official C++ transformer before this workflow.

Two IA variants are reported:
  linear_ia: linear matter power in BK-corrected-like IA terms.
  halofit_ia: model-specific CLASS/HALOFIT matter power in IA terms only.
The GG lensing term remains native linear Weyl in both variants; therefore this
is still not a final nonlinear modified-gravity KiDS likelihood.
"""
from pathlib import Path
import importlib.util,json,time
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from astropy.io import fits

S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
KDIR=Path("cosmosis-standard-library/likelihood/KiDS-Legacy")
WDIR=Path("cosmosis-standard-library/shear/cosebis/WnLog")
OUT=Path("output/kids_nuisance"); OUT.mkdir(parents=True,exist_ok=True)

ELL=np.geomspace(2.0,2.0e4,420)
NZGRID=300
FR=np.array([0.158,0.198,0.206,0.258,0.207,0.026])
LOGMPIV=13.5

# The public standard-library uncorrelated coordinates and Gaussian means.
DZ_U_MEAN=np.array([2.5374,-2.44484,-1.37982,-0.494275,1.40234,6.26323])
IA_U_MEAN=np.array([
  19.77069603580704893,30.74866003109325874,
  44.30836735636642,-32.06136674281973,143.92123862428525,
  -137.819672863503,-202.85133168708262,419.8363708097586
])
DZCOV=np.loadtxt(KDIR/"Nz_covariance.txt")
IACOV=np.loadtxt(KDIR/"massdep_cov.txt")
LDZ=np.linalg.cholesky(DZCOV)
LIA=np.linalg.cholesky(IACOV)

def c1rho():
    C1_M_sun=5e-14; M_sun=1.9891e30; Mpc=3.0857e22
    C1_SI=C1_M_sun/M_sun*Mpc**3
    G=6.67384e-11; H=100*1000.0/Mpc
    rho=3*H**2/(8*np.pi*G)
    return C1_SI*rho
C1RHO=c1rho()

with fits.open(DATA) as f:
    dd=f["En"].data
    D=np.asarray(dd["VALUE"],float)
    ROWS=[(int(r["BIN1"])-1,int(r["BIN2"])-1,int(r["ANGBIN"])-1) for r in dd]
    COV=np.asarray(f["COVMAT"].data,float)
    nz=f["NZ_SOURCE"].data
    ZSRC=np.asarray(nz["Z_MID"],float)
    NZFID=np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])
LCOV=np.linalg.cholesky(COV)

# trapezoidal ell integration weights
DX=np.diff(ELL)
TRAP=np.empty_like(ELL)
TRAP[0]=DX[0]/2; TRAP[-1]=DX[-1]/2
TRAP[1:-1]=(DX[:-1]+DX[1:])/2
WN=np.zeros((6,len(ELL)))
for n in range(1,7):
    tab=np.loadtxt(WDIR/f"WnLog{n}-2.00-300.00.table")
    WN[n-1]=np.interp(np.log(ELL),tab[:,0],tab[:,1],left=0.0,right=0.0)
KCOSEBI=WN*(ELL*TRAP/(2*np.pi))[None,:]

hL=0.6971482083084993; OmL=(0.022083219194622913+0.12299536722293603)/hL**2
hC=0.6856859; OmC=(0.02240637+0.11824151)/hC**2
MODELS=[
 b.load_model("late300","output/weyl_spectra",hL,OmL),
 b.load_model("local021","output/weyl_spectra",hC,OmC)
]

def tail(z,y):
    # integral from z to zmax
    dz=np.diff(z)
    t=0.5*(y[:-1]+y[1:])*dz
    out=np.zeros_like(y)
    out[:-1]=np.cumsum(t[::-1])[::-1]
    return out

class Likelihood:
    def __init__(self,m,variant):
        self.m=m; self.variant=variant
        zmax=min(float(m["zq"].max()),float(m["zp"].max()),float(ZSRC.max()))
        zlo=max(0.01,float(ZSRC[ZSRC>0].min()))
        self.z=np.linspace(zlo,zmax,NZGRID)
        self.H=np.interp(self.z,m["bg"]["z"],m["bg"]["H"])
        self.chi=np.interp(self.z,m["bg"]["z"],m["bg"]["chi"])
        # trapezoidal chi integration weights
        dc=np.diff(self.chi)
        self.wchi=np.empty_like(self.chi)
        self.wchi[0]=dc[0]/2; self.wchi[-1]=dc[-1]/2
        self.wchi[1:-1]=(dc[:-1]+dc[1:])/2
        zm=np.broadcast_to(self.z[None,:],(len(ELL),len(self.z)))
        km=(ELL[:,None]+0.5)/np.maximum(self.chi[None,:],1e-12)
        self.Q=b.eval_cube(m["Iq"],zm,km)
        self.Plin=b.eval_cube(m["Ip"],zm,km)
        self.Pnl=b.eval_cube(m["Ipnl"],zm,km)
        ok=np.isfinite(self.Q)&np.isfinite(self.Plin)
        self.Q=np.where(ok,self.Q,0.0)
        self.Plin=np.where(np.isfinite(self.Plin),self.Plin,0.0)
        self.Pnl=np.where(np.isfinite(self.Pnl),self.Pnl,self.Plin)
        self.valid_fraction=float(np.mean(ok))

    def components(self,t):
        # t are standardized offsets around the public uncorrelated-prior means.
        udz=DZ_U_MEAN+np.asarray(t[:6])
        uia=IA_U_MEAN+np.asarray(t[6:14])
        dz=LDZ@udz
        phys=LIA@uia
        A=float(phys[0]); beta=float(phys[1]); logM=np.asarray(phys[2:8])
        coeff=FR*10.0**((logM-LOGMPIV)*beta)

        nzs=[]
        for i in range(6):
            ff=interp1d(ZSRC,NZFID[i],bounds_error=False,fill_value=0.0)
            # same additive convention as the DES audit: positive bias raises mean z
            y=np.maximum(ff(self.z-dz[i]),0.0)
            norm=np.trapezoid(y,self.z)
            if not np.isfinite(norm) or norm<=0: return None
            nzs.append(y/norm)
        nzs=np.asarray(nzs)
        g=np.empty_like(nzs)
        for i,y in enumerate(nzs):
            g[i]=tail(self.z,y)-self.chi*tail(self.z,y/np.maximum(self.chi,1e-12))
        nchi=nzs*self.H[None,:]

        PIA=self.Plin if self.variant=="linear_ia" else self.Pnl
        # BK-corrected growth using the safely linear reference k=0.05 1/Mpc.
        I=self.m["Ip"] if self.variant=="linear_ia" else self.m["Ipnl"]
        pref=np.exp(I(np.column_stack([self.z,np.full_like(self.z,np.log(0.05))])))
        p0=float(np.exp(I([[0.0,np.log(0.05)]])[0]))
        growth=np.maximum(np.sqrt(pref/p0),1e-4)
        F=-A*C1RHO*self.m["Om"]/growth

        # Native-Weyl GG and perfectly correlated Weyl-matter GI approximation.
        X=np.sqrt(np.maximum(self.Q*PIA,0.0))
        w=self.wchi
        GG=np.einsum("iz,jz,lz,z->lij",g,g,self.Q,w,optimize=True)
        nc=nchi*coeff[:,None]
        GI=np.einsum("iz,jz,lz,z,z->lij",g,nc,X,F/np.maximum(self.chi,1e-12),w,optimize=True)
        GI=GI+np.swapaxes(GI,1,2)
        II=np.einsum("iz,jz,lz,z,z->lij",nc,nc,PIA,F**2/np.maximum(self.chi,1e-12)**2,w,optimize=True)
        CL=GG+GI+II
        # E[n,i,j] = sum_l K[n,l] C_l[i,j]
        E=np.einsum("nl,lij->nij",KCOSEBI,CL,optimize=True)
        th=np.array([E[n,min(i,j),max(i,j)] for i,j,n in ROWS])
        return th,dict(dz=dz.tolist(),A=A,beta=beta,logM=logM.tolist(),coeff=coeff.tolist())

    def objective(self,t):
        if np.any(np.abs(t)>5): return 1e100
        q=self.components(t)
        if q is None: return 1e100
        th,_=q
        y=np.linalg.solve(LCOV,D-th)
        # -2 log posterior up to constants: chi2_data + Gaussian standardized prior penalty
        return float(y@y + np.dot(t,t))

    def fit(self):
        starts=[np.zeros(14)]
        s=np.zeros(14); s[:6]=0.35*np.array([1,-1,1,-1,1,-1]); starts.append(s)
        s2=np.zeros(14); s2[6:]=0.25*np.array([1,-1,1,-1,1,-1,1,-1]); starts.append(s2)
        best=None
        records=[]
        for q,x0 in enumerate(starts):
            t0=time.time()
            r=minimize(self.objective,x0,method="L-BFGS-B",bounds=[(-5,5)]*14,
                       options=dict(maxiter=240,ftol=2e-9,maxls=30,maxfun=6500))
            records.append(dict(start=q,fun=float(r.fun),success=bool(r.success),nit=int(r.nit),nfev=int(r.nfev),message=str(r.message),seconds=time.time()-t0))
            if best is None or r.fun<best.fun: best=r
        # one derivative-free polish around best, limited budget
        rp=minimize(self.objective,best.x,method="Powell",bounds=[(-5,5)]*14,
                    options=dict(maxiter=45,xtol=2e-4,ftol=2e-5,maxfev=4500))
        records.append(dict(start="powell_polish",fun=float(rp.fun),success=bool(rp.success),nit=int(rp.nit),nfev=int(rp.nfev),message=str(rp.message)))
        if rp.fun<best.fun: best=rp
        th,phys=self.components(best.x)
        y=np.linalg.solve(LCOV,D-th)
        chi_data=float(y@y); penalty=float(np.dot(best.x,best.x))
        return dict(neg2logpost=float(best.fun),chi2_data=chi_data,prior_penalty=penalty,
                    t=best.x.tolist(),physical=phys,success=bool(best.success),
                    optimizer_records=records,valid_fraction=self.valid_fraction)

def main():
    out={"status":"KiDS-Legacy 126-point native-Weyl correlated-nuisance profile",
         "scientific_scope":{
           "GG":"native linear hi_class Weyl",
           "photoz":"6-bin correlated KiDS-Legacy prior",
           "IA":"public mass-dependent IA parameterization; two matter-power variants",
           "cosebi":"public Wn tables 2'-300', n=1..6; fast transform prevalidated against official C++",
           "nonlinear_GG":"not included; no validated nonlinear SDMC Weyl prescription",
           "baryons":"not varied because nonlinear GG is not used"
         },"results":{}}
    for variant in ["linear_ia","halofit_ia"]:
        vr={}
        for m in MODELS:
            print("KIDS_NUISANCE_START",variant,m["name"],flush=True)
            fit=Likelihood(m,variant).fit()
            vr[m["name"]]=fit
            print("KIDS_NUISANCE_DONE",variant,m["name"],json.dumps(fit,sort_keys=True),flush=True)
        vr["delta_chi2_data_late300_minus_local021"]=vr["late300"]["chi2_data"]-vr["local021"]["chi2_data"]
        vr["delta_neg2logpost_late300_minus_local021"]=vr["late300"]["neg2logpost"]-vr["local021"]["neg2logpost"]
        out["results"][variant]=vr
        print("KIDS_NUISANCE_RESULT",variant,json.dumps(vr,sort_keys=True),flush=True)
    (OUT/"late300_kids_nuisance_profile.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("LATE300_KIDS_NUISANCE_SUMMARY",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
