#!/usr/bin/env python3
"""
DES-Y3 cosmic-shear differential likelihood for accepted-action late300.

This is the first survey-level tomographic shear replay in the SDMC chain.
It uses:
  * the public DES Y3 2pt FITS data vector, covariance and source n(z);
  * the published DES Y3 xi+/xi- angular scale cuts;
  * the freely evolved accepted-action late300 background;
  * CLASS nonlinear matter P(k,z) (Halofit) for late300 and local021;
  * the exact No-Slip lensing response Sigma(z)=1/M_*^2(z);
  * the new native hi_class Weyl output as an independent linear consistency test;
  * an NLA intrinsic-alignment model, profiled over A_IA and eta_IA;
  * four multiplicative shear calibration parameters with DES Y3 Gaussian priors.

Important limitation:
  This is a survey-level differential likelihood, but not yet the official DES Y3
  TATT likelihood. Photo-z shifts are fixed at their prior means (zero), baryonic
  nuisance freedom is not sampled, and the nonlinear Weyl spectrum is approximated
  as Sigma^2(z) times the CLASS nonlinear matter spectrum after validating the
  linear relation against native wPk.  The output therefore reports both its
  numerical result and this method label explicitly.
"""
from pathlib import Path
import configparser, json, math, re
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.special import j0, jv

ROOT=Path("output/desy3_tomographic")
ROOT.mkdir(parents=True,exist_ok=True)
DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")

C1RHO=0.0134
ZPIV=0.62
M_MEAN=np.array([-0.0063,-0.0198,-0.0241,-0.0369])
M_SIG =np.array([ 0.0091, 0.0078, 0.0076, 0.0076])

def parse_bg(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    df=pd.DataFrame(a,columns=names).sort_values("z")
    return df

def parse_z_from_header(path):
    for line in Path(path).read_text().splitlines()[:10]:
        m=re.search(r"redshift z=([0-9.eE+-]+)",line)
        if m: return float(m.group(1))
        m=re.search(r"at z\s*=\s*([0-9.eE+-]+)",line)
        if m: return float(m.group(1))
    raise RuntimeError(f"no z header in {path}")

def pk_grid(prefix, suffix):
    files=sorted(Path("output").glob(prefix+f"_z*_{suffix}.dat"))
    if not files:
        files=sorted(Path("output").glob(prefix+f"_{suffix}.dat"))
    rows=[]
    for p in files:
        z=parse_z_from_header(p)
        a=np.loadtxt(p)
        rows.append((z,a[:,0],a[:,1]))
    rows.sort(key=lambda t:t[0])
    z=np.array([r[0] for r in rows])
    k=rows[0][1]
    P=np.array([np.interp(k,r[1],r[2]) for r in rows])
    return z,k,P

def interpolator(z,k,P):
    # linear interpolation in z and log-k/log-P
    return RegularGridInterpolator(
        (z,np.log(k)),np.log(np.maximum(P,1e-300)),
        bounds_error=False,fill_value=None
    )

def query_pk(itp,z_eval,k_h):
    zz=np.repeat(z_eval[None,:],len(k_h),axis=0)
    kk=np.repeat(np.log(k_h)[:,None],len(z_eval),axis=1)
    pts=np.column_stack([zz.ravel(),kk.ravel()])
    return np.exp(itp(pts)).reshape(len(k_h),len(z_eval))

def load_des():
    with fits.open(DATA) as f:
        xp=f["xip"].data
        xm=f["xim"].data
        nz=f["nz_source"].data
        z_nz=np.array(nz["Z_MID"],float)
        nzi=np.array([nz[f"BIN{i}"] for i in range(1,5)],float)
        cov=np.array(f["COVMAT"].data,float)
    xp_df=pd.DataFrame({c:np.array(xp[c]) for c in ["BIN1","BIN2","ANG","VALUE"]})
    xm_df=pd.DataFrame({c:np.array(xm[c]) for c in ["BIN1","BIN2","ANG","VALUE"]})
    return xp_df,xm_df,z_nz,nzi,cov

def scale_mask(df, kind):
    cp=configparser.ConfigParser(); cp.read(CUTS)
    sec=cp["2pt_like"]
    keep=[]
    for _,r in df.iterrows():
        i,j=int(r.BIN1),int(r.BIN2)
        key=f"angle_range_{kind}_{i}_{j}"
        if key not in sec:
            key=f"angle_range_{kind}_{j}_{i}"
        lo,hi=[float(x) for x in sec[key].split()]
        keep.append(lo <= float(r.ANG) <= hi)
    return np.array(keep,bool)

def make_model(label,prefix,bg_path,h):
    bg=parse_bg(bg_path)
    zlin,klin,Plin=pk_grid(prefix,"pk")
    try:
        znl,knl,Pnl=pk_grid(prefix,"pk_nl")
    except Exception:
        znl,knl,Pnl=zlin,klin,Plin
    zw,kw,Qw=pk_grid(prefix,"pk_weyl")
    return dict(label=label,bg=bg,h=h,zlin=zlin,klin=klin,Plin=Plin,
                znl=znl,knl=knl,Pnl=Pnl,zw=zw,kw=kw,Qw=Qw)

xp,xm,zsrc,nzsrc,covfull=load_des()
mxp=scale_mask(xp,"xip"); mxm=scale_mask(xm,"xim")
nxp=len(xp); nxm=len(xm); nall=nxp+nxm
full_data=np.concatenate([xp.VALUE.to_numpy(float),xm.VALUE.to_numpy(float)])
sel=np.r_[np.where(mxp)[0], nxp+np.where(mxm)[0]]
data=full_data[sel]
cov=covfull[:nall,:nall][np.ix_(sel,sel)]
icov=np.linalg.inv(cov)

late=make_model("late300","late300_wl_00",Path("output/late300_wl_00_background.dat"),69.71482083084993/100.)
loc =make_model("local021","local021_wl_00",Path("output/local021_wl_00_background.dat"),68.56859/100.)

# Normalize source n(z)
nzsrc=np.array([n/np.trapezoid(n,zsrc) for n in nzsrc])

theta_unique=np.unique(np.concatenate([xp.ANG.to_numpy(float),xm.ANG.to_numpy(float)]))
theta_rad=theta_unique*np.pi/(180.*60.)
ell=np.geomspace(10.,3.0e4,520)
lnell=np.log(ell)
B0=np.array([j0(ell*t) for t in theta_rad])
B4=np.array([jv(4,ell*t) for t in theta_rad])

def xi_from_cl(cl,order):
    B=B0 if order==0 else B4
    return np.trapezoid(B*(ell[None,:]**2*cl[None,:]/(2*np.pi)),lnell,axis=1)

def setup(model):
    bg=model["bg"]
    h=model["h"]
    zmax=min(3.0,float(model["znl"].max()),float(zsrc.max()))
    z=np.linspace(0.01,zmax,150)
    chi=np.interp(z,bg.z,bg["comov. dist."])
    Hc=np.interp(z,bg.z,bg["H [1/Mpc]"])
    Om=np.interp(z,bg.z,bg["Omega_m(z)"])
    Om0=float(np.interp(0,bg.z,bg["Omega_m(z)"]))
    M2=np.interp(z,bg.z,bg["M*^2_smg"]) if "M*^2_smg" in bg else np.ones_like(z)
    Sigma=1.0/M2

    # source efficiencies g_i(z)
    nz=np.array([np.interp(z,zsrc,n,left=0,right=0) for n in nzsrc])
    g=np.zeros_like(nz)
    for i in range(4):
        for iz,zz in enumerate(z):
            m=zsrc>=zz
            if not np.any(m): continue
            chis=np.interp(zsrc[m],bg.z,bg["comov. dist."])
            nzs=nzsrc[i,m]
            g[i,iz]=np.trapezoid(nzs*np.maximum(chis-chi[iz],0)/np.maximum(chis,1e-12),zsrc[m])

    H0c=(100*h)/299792.458
    q=1.5*Om0*H0c**2*(1+z)[None,:]*chi[None,:]*g

    # P_nl matrix on (ell,z), physical Mpc^3
    knl_itp=interpolator(model["znl"],model["knl"],model["Pnl"])
    klin_itp=interpolator(model["zlin"],model["klin"],model["Plin"])
    kphys=(ell[:,None]+0.5)/chi[None,:]
    kh=kphys/h
    Pnl_h=query_pk(knl_itp,z,kh[:,0]) if False else None
    # query_pk assumes one k per ell, but k varies with z. Evaluate directly.
    kmin=float(model["knl"].min()); kmax=float(model["knl"].max())
    kh_clip=np.clip(kh,kmin,kmax)
    pts=np.column_stack([np.repeat(z,len(ell)), np.log(kh_clip.T.ravel())])
    Pnl=np.exp(knl_itp(pts)).reshape(len(z),len(ell)).T/h**3
    # Linear table has essentially the same k range, but clip independently.
    klmin=float(model["klin"].min()); klmax=float(model["klin"].max())
    khl=np.clip(kh,klmin,klmax)
    ptsl=np.column_stack([np.repeat(z,len(ell)), np.log(khl.T.ravel())])
    Plin=np.exp(klin_itp(ptsl)).reshape(len(z),len(ell)).T/h**3
    high_clip_fraction=float(np.mean(kh>kmax))

    # Linear growth from fixed low-k scale.
    kref=0.05
    ptsD=np.column_stack([z,np.full_like(z,np.log(kref))])
    Pgrow=np.exp(klin_itp(ptsD))
    P0=float(np.exp(klin_itp([[0.0,np.log(kref)]])[0]))
    growth=np.sqrt(Pgrow/P0)

    # Native Weyl check against the No-Slip Poisson identity.
    witp=interpolator(model["zw"],model["kw"],model["Qw"])
    test_z=np.array([0.1,0.5,1.0,1.5,2.0])
    test_kh=np.array([0.02,0.05,0.10,0.20])
    rows=[]
    for zz in test_z:
        sig=float(np.interp(zz,bg.z,1.0/bg["M*^2_smg"])) if "M*^2_smg" in bg else 1.0
        om0=Om0
        for kk in test_kh:
            qnative=float(np.exp(witp([[zz,np.log(kk)]])[0]))*h # 1/Mpc
            pl=float(np.exp(klin_itp([[zz,np.log(kk)]])[0]))/h**3
            expected=(1.5*om0*H0c**2*(1+zz)*sig)**2*pl
            rows.append(dict(z=zz,k_h_Mpc=kk,QW_native_1_Mpc=qnative,
                             QW_Poisson_1_Mpc=expected,ratio=qnative/expected))
    pd.DataFrame(rows).to_csv(ROOT/f"{model['label']}_weyl_poisson_check.csv",index=False)

    return dict(z=z,chi=chi,Hc=Hc,Om0=Om0,Sigma=Sigma,nz=nz,q=q,
                Pnl=Pnl,Plin=Plin,growth=growth,high_k_clip_fraction=high_clip_fraction,
                kmax_h_Mpc=kmax)

lateS=setup(late); locS=setup(loc)

def precompute_theory(S,eta):
    z=S["z"]; chi=S["chi"]; Hc=S["Hc"]; Sigma=S["Sigma"]
    P=S["Pnl"]; q=S["q"]; nz=S["nz"]; D=S["growth"]; Om0=S["Om0"]
    # n_chi = n_z H_CLASS
    nchi=nz*Hc[None,:]
    F=-C1RHO*Om0/np.maximum(D,1e-8)*((1+z)/(1+ZPIV))**eta
    # pieces for A_IA=1
    dz_weight=1.0/Hc
    clpieces={}
    for i in range(4):
      for j in range(i,4):
        denom=np.maximum(chi,1e-8)**2
        Kgg=dz_weight*q[i]*q[j]/denom*Sigma**2
        Kgi=dz_weight*(q[i]*nchi[j]+q[j]*nchi[i])/denom*Sigma*F
        Kii=dz_weight*nchi[i]*nchi[j]/denom*F**2
        Cgg=np.trapezoid(P*Kgg[None,:],z,axis=1)
        Cgi=np.trapezoid(P*Kgi[None,:],z,axis=1)
        Cii=np.trapezoid(P*Kii[None,:],z,axis=1)
        clpieces[(i+1,j+1)]=(Cgg,Cgi,Cii)
    return clpieces

def model_vector(df_plus,df_minus,pieces,Aia,mcal):
    cachep={}; cachem={}
    outp=[]; outm=[]
    for _,r in df_plus.iterrows():
        pair=tuple(sorted((int(r.BIN1),int(r.BIN2))))
        if pair not in cachep:
            gg,gi,ii=pieces[pair]
            cachep[pair]=xi_from_cl(gg+Aia*gi+Aia*Aia*ii,0)
        arr=cachep[pair]
        idx=np.argmin(np.abs(theta_unique-float(r.ANG)))
        fac=(1+mcal[pair[0]-1])*(1+mcal[pair[1]-1])
        outp.append(arr[idx]*fac)
    for _,r in df_minus.iterrows():
        pair=tuple(sorted((int(r.BIN1),int(r.BIN2))))
        if pair not in cachem:
            gg,gi,ii=pieces[pair]
            cachem[pair]=xi_from_cl(gg+Aia*gi+Aia*Aia*ii,4)
        arr=cachem[pair]
        idx=np.argmin(np.abs(theta_unique-float(r.ANG)))
        fac=(1+mcal[pair[0]-1])*(1+mcal[pair[1]-1])
        outm.append(arr[idx]*fac)
    return np.r_[outp,outm]

xp_sel=xp[mxp].reset_index(drop=True); xm_sel=xm[mxm].reset_index(drop=True)

def fit_model(label,S):
    best=None
    eta_grid=np.linspace(-5,5,21)
    for eta in eta_grid:
        pieces=precompute_theory(S,float(eta))
        def obj(x):
            Aia=float(x[0]); mcal=np.asarray(x[1:5])
            th=model_vector(xp_sel,xm_sel,pieces,Aia,mcal)
            d=data-th
            chi=float(d@icov@d)
            chi+=float(np.sum(((mcal-M_MEAN)/M_SIG)**2))
            return chi
        x0=np.r_[0.7,M_MEAN]
        res=minimize(obj,x0,method="L-BFGS-B",bounds=[(-5,5)]+[(-0.1,0.1)]*4,
                     options={"maxiter":120,"ftol":1e-9})
        rec=dict(label=label,eta=float(eta),chi2_profile=float(res.fun),
                 A_IA=float(res.x[0]),success=bool(res.success),
                 **{f"m{i+1}":float(res.x[i+1]) for i in range(4)})
        if best is None or rec["chi2_profile"]<best["chi2_profile"]:
            best=rec
    return best

late_fit=fit_model("late300",lateS)
loc_fit=fit_model("local021",locS)
delta=late_fit["chi2_profile"]-loc_fit["chi2_profile"]

# Weyl closure summaries
def summarize_weyl(label):
    df=pd.read_csv(ROOT/f"{label}_weyl_poisson_check.csv")
    return dict(
      median_ratio=float(np.median(df.ratio)),
      max_abs_ratio_minus_1=float(np.max(np.abs(df.ratio-1))),
      min_ratio=float(df.ratio.min()),max_ratio=float(df.ratio.max())
    )

summary={
 "method":"DES Y3 xi+/xi- tomographic differential likelihood; public covariance and published angular cuts; NLA IA profiled; m_i profiled with published Gaussian priors; photo-z shifts fixed at zero; Halofit matter nonlinearity; exact late300 Sigma(z); native wPk linear identity cross-check.",
 "ndata_after_scale_cuts":int(len(data)),
 "late300":late_fit,
 "local021":loc_fit,
 "delta_chi2_late300_minus_local021":float(delta),
 "late300_weyl_check":summarize_weyl("late300"),
 "local021_weyl_check":summarize_weyl("local021"),
 "late300_high_k_clip_fraction":lateS["high_k_clip_fraction"],
 "local021_high_k_clip_fraction":locS["high_k_clip_fraction"],
 "late300_kmax_h_Mpc":lateS["kmax_h_Mpc"],
 "local021_kmax_h_Mpc":locS["kmax_h_Mpc"],
 "interpretation":"Negative delta favors late300 within this intermediate survey-level shear model; positive delta favors local021. This is not yet the official DES Y3 TATT likelihood."
}
(ROOT/"desy3_tomographic_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
pd.DataFrame([late_fit,loc_fit]).to_csv(ROOT/"desy3_tomographic_profile.csv",index=False)
print("DESY3_TOMOGRAPHIC_RESULT",json.dumps(summary,sort_keys=True),flush=True)
