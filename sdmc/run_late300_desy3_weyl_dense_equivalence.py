#!/usr/bin/env python3
"""
Dense-Hankel equivalence test for the late300 DES-Y3 shear calculation.

Purpose: determine whether the discrepancy between the recent direct-native-Weyl
pilot and the independent Sigma^2 P_nl tomographic branch is caused by the
Hankel transform / source-kernel implementation or by the physics prescription.

This script:
  * uses the already-validated native Q_W(k,z) from hi_class,
  * constructs Q_W,NL = Q_W,lin * P_nl/P_lin,
  * projects GG, GI and II in the exact algebraic form equivalent to
    scale-independent No-Slip Sigma^2 P_nl,
  * uses the same dense 50k-point Hankel transform and the same NLA eta/A_IA +
    four shear-calibration nuisance profile as the independent shear branch,
  * compares native-Q_W and reconstructed Sigma^2 P predictions directly.

It is a numerical equivalence/debugging test, not a final nonlinear-MG likelihood.
"""
from pathlib import Path
import configparser,json,re
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.special import j0,jv

ROOT=Path("output/weyl_dense_equivalence"); ROOT.mkdir(parents=True,exist_ok=True)
SPEC=Path("output/weyl_spectra")
DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
C1RHO=0.0134; ZPIV=0.62
M_MEAN=np.array([-0.0063,-0.0198,-0.0241,-0.0369])
M_SIG=np.array([0.0091,0.0078,0.0076,0.0076])
ELL=np.geomspace(2.,3e4,1400)
ELLX=np.geomspace(2.,3e4,50000); LNELLX=np.log(ELLX)

def zhdr(p):
    for line in Path(p).read_text().splitlines()[:10]:
        m=re.search(r"(?:at\s+)?(?:redshift\s+)?z\s*=\s*([0-9.eE+\-]+)",line,re.I)
        if m:return float(m.group(1))
    raise RuntimeError(f"no z header {p}")

def parse_bg(p):
    lines=Path(p).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(p); return pd.DataFrame(a,columns=names).sort_values("z")

def cube(prefix,kind,h):
    if kind=="Q":
        files=sorted(SPEC.glob(prefix+"*pk_weyl.dat"))
    elif kind=="Pl":
        files=sorted(p for p in SPEC.glob(prefix+"*pk.dat") if "pk_weyl" not in p.name and "pk_nl" not in p.name)
    elif kind=="Pn":
        files=sorted(SPEC.glob(prefix+"*pk_nl.dat"))
    rows=[]
    for p in files:
        z=zhdr(p); a=np.loadtxt(p)
        k=a[:,0]*h
        y=a[:,1]*h if kind=="Q" else a[:,1]/h**3
        rows.append((z,k,y))
    rows.sort(key=lambda t:t[0])
    if not rows: raise RuntimeError((prefix,kind))
    z=np.array([r[0] for r in rows]); k=rows[0][1]
    P=np.array([np.exp(np.interp(np.log(k),np.log(r[1]),np.log(r[2]))) for r in rows])
    return z,k,P

def interp(z,k,P):
    return RegularGridInterpolator((z,np.log(k)),np.log(np.maximum(P,1e-300)),
                                   bounds_error=False,fill_value=np.nan)

def evalit(I,zmat,kmat):
    pts=np.column_stack([zmat.ravel(),np.log(kmat.ravel())])
    return np.exp(I(pts)).reshape(kmat.shape)

def load_model(name,h,Om):
    p00=SPEC/f"{name}_00_background.dat"
    bg=parse_bg(p00)
    zq,kq,Q=cube(f"{name}_00_","Q",h)
    zp,kp,Pl=cube(f"{name}_00_","Pl",h)
    zn,kn,Pn=cube(f"{name}_nl_00_","Pn",h)
    return dict(name=name,h=h,Om=Om,bg=bg,zq=zq,zp=zp,zn=zn,
                Iq=interp(zq,kq,Q),Ip=interp(zp,kp,Pl),In=interp(zn,kn,Pn))

def load_des():
    with fits.open(DATA) as f:
        xp=f["xip"].data; xm=f["xim"].data; nz=f["nz_source"].data
        zsrc=np.asarray(nz["Z_MID"],float)
        nzsrc=np.asarray([nz[f"BIN{i}"] for i in range(1,5)],float)
        cov=np.asarray(f["COVMAT"].data,float); hdr=f["COVMAT"].header
        starts={}
        for ii in range(20):
            nk=f"NAME_{ii}"; sk=f"STRT_{ii}"
            if nk in hdr and sk in hdr: starts[str(hdr[nk]).strip().lower()]=int(hdr[sk])
    cp=configparser.ConfigParser(); cp.read(CUTS); sec=cp["2pt_like"]
    rows=[]; idx=[]
    for kind,dat in [("xip",xp),("xim",xm)]:
        st=starts[kind]
        for q,r in enumerate(dat):
            i,j=int(r["BIN1"]),int(r["BIN2"])
            key=f"angle_range_{kind}_{min(i,j)}_{max(i,j)}"
            lo,hi=[float(v) for v in sec[key].split()]
            if lo<=float(r["ANG"])<=hi:
                rows.append(dict(kind=kind,b1=i,b2=j,ang=float(r["ANG"]),val=float(r["VALUE"])))
                idx.append(st+q)
    idx=np.asarray(idx,int)
    return dict(rows=rows,data=np.array([r["val"] for r in rows]),cov=cov[np.ix_(idx,idx)],
                zsrc=zsrc,nzsrc=nzsrc/np.trapezoid(nzsrc,zsrc,axis=1)[:,None])

DES=load_des()
THETA=np.array(sorted(set(r["ang"] for r in DES["rows"])))
THRAD=THETA*np.pi/(180*60)
B0=np.array([j0(ELLX*t) for t in THRAD])
B4=np.array([jv(4,ELLX*t) for t in THRAD])

def xiform(cl,order):
    cli=np.interp(LNELLX,np.log(ELL),cl)
    B=B0 if order==0 else B4
    return np.trapezoid(B*(ELLX[None,:]**2*cli[None,:]/(2*np.pi)),LNELLX,axis=1)

def setup(model):
    bg=model["bg"]; zmax=min(3.0,float(model["zq"].max()),float(model["zp"].max()),float(model["zn"].max()),float(DES["zsrc"].max()))
    z=np.linspace(0.01,zmax,150)
    chi=np.interp(z,bg.z,bg["comov. dist."])
    H=np.interp(z,bg.z,bg["H [1/Mpc]"])
    Om0=float(np.interp(0,bg.z,bg["Omega_m(z)"]))
    M2=np.interp(z,bg.z,bg["M*^2_smg"]) if "M*^2_smg" in bg else np.ones_like(z)
    Sigma=1/M2
    nz=np.array([np.interp(z,DES["zsrc"],n,left=0,right=0) for n in DES["nzsrc"]])
    g=np.zeros_like(nz)
    for i in range(4):
        for iz,zz in enumerate(z):
            m=DES["zsrc"]>=zz
            if not np.any(m):continue
            chis=np.interp(DES["zsrc"][m],bg.z,bg["comov. dist."])
            g[i,iz]=np.trapezoid(DES["nzsrc"][i,m]*np.maximum(chis-chi[iz],0)/np.maximum(chis,1e-12),DES["zsrc"][m])
    nchi=nz*H[None,:]
    zmat=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    k=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
    Ql=evalit(model["Iq"],zmat,k)
    Pl=evalit(model["Ip"],zmat,k)
    Pn=evalit(model["In"],zmat,k)
    valid=np.isfinite(Ql)&np.isfinite(Pl)&np.isfinite(Pn)
    Ql=np.where(valid,Ql,0); Pl=np.where(valid,Pl,0); Pn=np.where(valid,Pn,Pl)
    boost=np.divide(Pn,Pl,out=np.ones_like(Pn),where=Pl>0)
    boost=np.clip(boost,0.2,50)
    Qn=Ql*boost
    Cross=np.sqrt(np.maximum(Qn*Pn,0))
    # Explicit Sigma^2 Pnl identity, for debugging only.
    H0=100*model["h"]/299792.458
    A=(1.5*Om0*H0*H0*(1+z)*Sigma)
    Qrecon=A[None,:]**2*Pn
    ratio=np.divide(Qn,Qrecon,out=np.full_like(Qn,np.nan),where=Qrecon>0)
    finite=ratio[np.isfinite(ratio)&(k>0.005)&(k<1.0)]
    kref=0.05
    pref=np.exp(model["Ip"](np.column_stack([z,np.full_like(z,np.log(kref))])))
    p0=float(np.exp(model["Ip"]([[0.,np.log(kref)]])[0])); D=np.sqrt(pref/p0)
    return dict(z=z,chi=chi,H=H,Om0=Om0,Sigma=Sigma,nz=nz,nchi=nchi,g=g,Pn=Pn,Qn=Qn,Cross=Cross,D=D,
                qratio_median=float(np.nanmedian(finite)),qratio_maxdev=float(np.nanmax(np.abs(finite-1))),
                valid=float(np.mean(valid)))

def pieces(S,eta):
    z=S["z"]; chi=S["chi"]; g=S["g"]; nchi=S["nchi"]; P=S["Pn"]; Q=S["Qn"]; Cross=S["Cross"]
    F=-C1RHO*S["Om0"]/np.maximum(S["D"],1e-8)*((1+z)/(1+ZPIV))**eta
    out={}
    for i in range(4):
      for j in range(i,4):
        gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
        gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*F/np.maximum(chi,1e-12))[None,:]*Cross,x=chi,axis=1)
        ii=np.trapezoid((nchi[i]*nchi[j]*F**2/np.maximum(chi,1e-12)**2)[None,:]*P,x=chi,axis=1)
        out[(i+1,j+1)]=(gg,gi,ii)
    return out

def vector(pc,A,mcal):
    cache={}; vals=[]
    for r in DES["rows"]:
        pair=tuple(sorted((r["b1"],r["b2"]))); order=0 if r["kind"]=="xip" else 4
        key=(pair,order)
        if key not in cache:
            gg,gi,ii=pc[pair]; cache[key]=xiform(gg+A*gi+A*A*ii,order)
        arr=cache[key]; q=np.argmin(np.abs(THETA-r["ang"]))
        fac=(1+mcal[pair[0]-1])*(1+mcal[pair[1]-1])
        vals.append(arr[q]*fac)
    return np.asarray(vals)

def fit(model,S):
    L=np.linalg.cholesky(DES["cov"])
    best=None
    for eta in np.linspace(-5,5,21):
        pc=pieces(S,float(eta))
        def fun(x):
            th=vector(pc,float(x[0]),np.asarray(x[1:5]))
            y=np.linalg.solve(L,DES["data"]-th)
            return float(y@y+np.sum(((x[1:5]-M_MEAN)/M_SIG)**2))
        res=minimize(fun,np.r_[0.7,M_MEAN],method="L-BFGS-B",bounds=[(-5,5)]+[(-.1,.1)]*4,
                     options=dict(maxiter=150,ftol=1e-10))
        rec=dict(eta=float(eta),chi2=float(res.fun),AIA=float(res.x[0]),m=[float(v) for v in res.x[1:5]],success=bool(res.success))
        if best is None or rec["chi2"]<best["chi2"]: best=rec
    return best

def main():
    hl=69.71482083084993/100; hc=68.56859/100
    late=load_model("late300",hl,(0.022083219194622913+0.12299536722293603)/hl**2)
    loc=load_model("local021",hc,(0.02240637+0.11824151)/hc**2)
    SL=setup(late); SC=setup(loc)
    fl=fit(late,SL); fc=fit(loc,SC)
    result=dict(
      ndata=len(DES["data"]),
      late300=fl,local021=fc,
      delta_chi2_late300_minus_local021=fl["chi2"]-fc["chi2"],
      late300_Qnl_vs_Sigma2Pnl=dict(median=SL["qratio_median"],max_abs_minus1=SL["qratio_maxdev"]),
      local021_Qnl_vs_Sigma2Pnl=dict(median=SC["qratio_median"],max_abs_minus1=SC["qratio_maxdev"]),
      projection="dense Hankel 50k; native QW times model-specific Halofit boost; NLA eta/AIA + DES m_i priors",
      interpretation="If this agrees with the independent Sigma^2 P_nl branch, the earlier small native-Weyl penalty was a coarse-Hankel numerical artefact."
    )
    (ROOT/"dense_equivalence.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("DENSE_WEYL_EQUIVALENCE",json.dumps(result,sort_keys=True),flush=True)

if __name__=="__main__": main()
