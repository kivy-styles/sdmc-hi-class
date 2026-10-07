#!/usr/bin/env python3
"""
KiDS-Legacy native-Weyl xi+/xi- large-angle audit for late300 vs local021.

This uses the public no-cut KiDS-Legacy 2PCF data vector/covariance and evaluates
native hi_class Weyl predictions directly against the FITS row ordering.

Two theory variants are carried:
  gg    : native linear Weyl shear only
  total : gg + GI + II using the published KiDS-Legacy nuisance central point
          reconstructed in the preceding covariance audit.

Because SDMC still lacks a calibrated nonlinear-Weyl prescription, the primary
scientific outputs are the release safe-cut result and progressively more
conservative angular cuts.  The all-scale result is diagnostic only.
"""
from pathlib import Path
import importlib.util, json, numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from scipy.special import jv
from astropy.io import fits

S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

DATA=Path("shear_release/data/KiDS_Legacy_cosmic_shear_data_release/KiDS_Legacy_xipm.fits")
SAFE=Path("shear_release/scale_cuts/safe_idx_kids.npy")
OUT=Path("output/kids_xipm_weyl"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2.0,2.0e4,1000)

# Published KiDS-Legacy nuisance central point reconstructed from the
# correlated-prior covariance transform in the preceding audit.
DZ=np.array([0.025825256822487035,-0.013466855669952453,0.001512059936331642,
             -0.008220981600339534,0.011414273687594095,0.05428260718820411])
IA_BIN=np.array([0.1452637378320543,0.39532862967317906,0.5602915856507814,
                 0.8349285079980582,0.779915332417869,0.11201047680838332])

hL=0.6971482083084993
OmL=(0.022083219194622913+0.12299536722293603)/hL**2
hC=0.6856859
OmC=(0.02240637+0.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)

with fits.open(DATA) as f:
    cov=np.asarray(f[1].data,float)
    xp=f[2].data.copy()
    xm=f[3].data.copy()
    nz=f[4].data.copy()

zsrc=np.asarray(nz["Z_MID"],float)
bins=np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])
data=np.concatenate([np.asarray(xp["VALUE"],float),np.asarray(xm["VALUE"],float)])

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def cls_for_model(model):
    zmax=min(model["zq"].max(),model["zp"].max(),float(zsrc.max()))
    z=np.linspace(max(0.01,float(zsrc[zsrc>0].min())),zmax,460)
    H=np.interp(z,model["bg"]["z"],model["bg"]["H"])
    chi=np.interp(z,model["bg"]["z"],model["bg"]["chi"])
    nzs=[]
    for i,arr in enumerate(bins):
        f=interp1d(zsrc,arr,bounds_error=False,fill_value=0.0)
        y=np.maximum(f(z-DZ[i]),0.0)
        n=np.trapezoid(y,z)
        if n<=0: raise RuntimeError(f"bad n(z) bin {i+1}")
        nzs.append(y/n)
    nzs=np.asarray(nzs)
    g=[]
    for y in nzs:
        q0=tail(z,y); q1=tail(z,y/np.maximum(chi,1e-12))
        g.append(q0-chi*q1)
    g=np.asarray(g)
    nchi=nzs*H[None,:]

    zm=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    km=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
    Q=b.eval_cube(model["Iq"],zm,km)
    P=b.eval_cube(model["Ip"],zm,km)
    ok=np.isfinite(Q)&np.isfinite(P)
    Q=np.where(ok,Q,0.0); P=np.where(ok,P,0.0)
    X=np.sqrt(np.maximum(Q*P,0.0))

    kref=0.05
    pref=np.exp(model["Ip"](np.column_stack([z,np.full_like(z,np.log(kref))])))
    p0=float(np.exp(model["Ip"]([[0.0,np.log(kref)]])[0]))
    D=np.maximum(np.sqrt(pref/p0),1e-4)
    Fbase=-b.C1RHO*model["Om"]/D
    Fi=IA_BIN[:,None]*Fbase[None,:]

    pairs={"gg":{},"total":{}}
    for i in range(6):
        for j in range(i,6):
            gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
            gi=np.trapezoid(
                ((g[i]*nchi[j]*Fi[j] + g[j]*nchi[i]*Fi[i])/
                 np.maximum(chi,1e-12))[None,:]*X,
                x=chi,axis=1
            )
            ii=np.trapezoid(
                (nchi[i]*nchi[j]*Fi[i]*Fi[j]/np.maximum(chi,1e-12)**2)[None,:]*P,
                x=chi,axis=1
            )
            pair=(i+1,j+1)
            pairs["gg"][pair]=gg
            pairs["total"][pair]=gg+gi+ii
    return pairs,dict(zmin=float(z.min()),zmax=float(z.max()),
                      valid_fraction=float(ok.mean()))

def xi_vector(pair_cls,variant):
    outp=[]; outm=[]
    for r in xp:
        pair=(min(int(r["BIN1"]),int(r["BIN2"])),max(int(r["BIN1"]),int(r["BIN2"])))
        cl=pair_cls[variant][pair]
        th=np.deg2rad(float(r["ANG"])/60.0)
        outp.append(np.trapezoid(ELL/(2*np.pi)*jv(0,ELL*th)*cl,x=ELL))
    for r in xm:
        pair=(min(int(r["BIN1"]),int(r["BIN2"])),max(int(r["BIN1"]),int(r["BIN2"])))
        cl=pair_cls[variant][pair]
        th=np.deg2rad(float(r["ANG"])/60.0)
        outm.append(np.trapezoid(ELL/(2*np.pi)*jv(4,ELL*th)*cl,x=ELL))
    return np.r_[outp,outm]

def safe_indices():
    a=np.load(SAFE)
    if a.dtype==bool:
        return np.flatnonzero(a)
    return np.asarray(a,int).ravel()

def angular_indices(pmin,mmin):
    mp=np.asarray(xp["ANG"],float)>=pmin
    mm=np.asarray(xm["ANG"],float)>=mmin
    return np.flatnonzero(np.r_[mp,mm])

def chi(theory,idx):
    idx=np.asarray(idx,int)
    C=cov[np.ix_(idx,idx)]
    r=data[idx]-theory[idx]
    try:
        L=np.linalg.cholesky(C)
        y=np.linalg.solve(L,r)
        q=float(y@y)
    except np.linalg.LinAlgError:
        q=float(r@np.linalg.pinv(C)@r)
    return q

models={}
meta={}
for m in (late,lcdm):
    pc,md=cls_for_model(m)
    models[m["name"]]={v:xi_vector(pc,v) for v in ("gg","total")}
    meta[m["name"]]=md

cuts={
    "all":np.arange(len(data)),
    "release_safe":safe_indices(),
    "c20_100":angular_indices(20.0,100.0),
    "c40_150":angular_indices(40.0,150.0),
    "c60_180":angular_indices(60.0,180.0),
    "c80_200":angular_indices(80.0,200.0),
    "c100_250":angular_indices(100.0,250.0),
}
out={
  "status":"KiDS-Legacy public 2PCF native-Weyl large-angle audit",
  "nuisance_central":{"dz":DZ.tolist(),"ia_bin_amplitude":IA_BIN.tolist()},
  "meta":meta,
  "data":{"n_xip":len(xp),"n_xim":len(xm),"n_total":len(data),
          "theta_p_unique":sorted(set(map(float,xp["ANG"]))),
          "theta_m_unique":sorted(set(map(float,xm["ANG"])))},
  "cuts":{}
}
for label,idx in cuts.items():
    rec={"ndata":int(len(idx))}
    for variant in ("gg","total"):
        qL=chi(models["late300"][variant],idx)
        qC=chi(models["local021"][variant],idx)
        rec[variant]={
            "late300_chi2":qL,"local021_chi2":qC,
            "delta_chi2_late300_minus_local021":qL-qC,
            "late300_chi2_per_point":qL/len(idx),
            "local021_chi2_per_point":qC/len(idx)
        }
    out["cuts"][label]=rec
    print("KIDS_XIPM_WEYL",label,json.dumps(rec,sort_keys=True),flush=True)

(OUT/"late300_kids_xipm_weyl.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
