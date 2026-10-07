#!/usr/bin/env python3
from pathlib import Path
import importlib.util, json
import numpy as np
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import RegularGridInterpolator

S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
REACT=Path("output/react_dense/react_late300/late300_react_target_spectrum.npz")
OUT=Path("output/kids_react_components"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2,1e5,1800)
C=299792.458
C1RHO=0.0134
ZREACT_MAX=2.5

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def source():
    with fits.open(DATA) as f:
        nz=f["NZ_SOURCE"].data
        return np.asarray(nz["Z_MID"],float), np.array([
            np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)
        ])

def make_boost_interpolator(zr,kr,boost):
    zr=np.asarray(zr,float); kr=np.asarray(kr,float); boost=np.asarray(boost,float)
    return RegularGridInterpolator(
        (zr,np.log(kr)), np.log(np.maximum(boost,1e-12)),
        bounds_error=False, fill_value=None
    )

def eval_boost(I,z,kh,zmax=ZREACT_MAX):
    z=np.asarray(z,float); kh=np.asarray(kh,float)
    zc=np.clip(z,0,zmax)
    # ReACT k range is finite; freeze the correction at the validated boundary.
    # This is explicit and recorded in metadata instead of silently reverting to CLASS HALOFIT.
    pts=np.column_stack([zc.ravel(),np.log(np.clip(kh.ravel(),KMIN,KMAX))])
    out=np.exp(I(pts)).reshape(z.shape)
    # Above the ReACT redshift ceiling, taper the nonlinear correction back to unity by z=3.
    w=np.clip((3.0-z)/0.5,0.0,1.0)
    out=1.0+w*(out-1.0)
    return out

def build_react():
    r=np.load(REACT)
    zr=np.asarray(r["z"],float); kr=np.asarray(r["k"],float)
    out={}
    for tag in ["us","ss"]:
        out[("late300",tag)]=make_boost_interpolator(zr,kr,r[f"Pnl_{tag}"]/r["Ptarget"])
    out[("local021","gr")]=make_boost_interpolator(
        zr,kr,r["Pnl_local021"]/r["Ptarget_local021"]
    )
    return r,zr,kr,out

RDATA,ZR,KR,BOOST=build_react()
KMIN=float(KR.min()); KMAX=float(KR.max())

def components(m,zsrc,nzbins,variant):
    zmax=min(m["zq"].max(),float(np.max(zsrc)))
    z=np.linspace(.01,zmax,420)
    H=np.interp(z,m["bg"]["z"],m["bg"]["H"])
    chi=np.interp(z,m["bg"]["z"],m["bg"]["chi"])

    nz=[]
    for a in nzbins:
        y=np.maximum(np.interp(z,zsrc,a,left=0,right=0),0)
        y/=np.trapezoid(y,z)
        nz.append(y)
    nz=np.asarray(nz); nchi=nz*H[None,:]
    g=np.asarray([tail(z,y)-chi*tail(z,y/np.maximum(chi,1e-12)) for y in nz])

    zm=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    km=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)  # 1/Mpc
    kh=km/m["h"]                                        # h/Mpc

    Q=b.eval_cube(m["Iq"],zm,km)
    Pl=b.eval_cube(m["Ip"],zm,km)
    ok=np.isfinite(Q)&np.isfinite(Pl)&(Pl>0)
    Q=np.where(ok,Q,0.0); Pl=np.where(ok,Pl,0.0)

    if variant=="linear":
        P=Pl; Quse=Q
        matter_tag="linear"; weyl_mode="linear"
    else:
        if m["name"]=="late300":
            matter_tag="ss" if variant.startswith("ss_") else "us"
            I=BOOST[("late300",matter_tag)]
        else:
            matter_tag="gr"
            I=BOOST[("local021","gr")]
        boost=eval_boost(I,zm,kh)
        P=Pl*boost

        H0=100*m["h"]/C
        Agr=(1.5*m["Om"]*H0**2*(1+z))[None,:]
        sig2=np.divide(Q,Agr**2*Pl,out=np.ones_like(Q),where=(Agr**2*Pl)>0)
        sig=np.sqrt(np.clip(sig2,0.2,2.0))

        if variant.endswith("_native"):
            sigeff=sig; weyl_mode="native_linear_Sigma"
        elif variant.endswith("_screen03"):
            x=kh/0.3
            sigeff=1.0+(sig-1.0)/(1.0+x*x)
            weyl_mode="No-Slip transition k_s=0.3 h/Mpc"
        elif variant.endswith("_fullscreen"):
            sigeff=np.ones_like(sig)
            weyl_mode="full nonlinear Weyl screen to GR Sigma=1"
        else:
            raise ValueError(variant)
        Quse=Agr**2*P*sigeff**2
        Quse=np.where(ok,Quse,0.0)

    pref=np.exp(m["Ip"](np.column_stack([z,np.full_like(z,np.log(.05))])))
    p0=float(np.exp(m["Ip"]([[0.0,np.log(.05)]])[0]))
    D=np.maximum(np.sqrt(pref/p0),1e-4)
    F=-C1RHO*m["Om"]/D
    X=np.sqrt(np.maximum(Quse*P,0.0))

    out={q:{"ell":ELL} for q in ["gg","gi","ii"]}
    for i in range(6):
        for j in range(i,6):
            key=f"bin_{j+1}_{i+1}"
            out["gg"][key]=np.trapezoid(g[i][None,:]*g[j][None,:]*Quse,x=chi,axis=1)
            out["gi"][key]=np.trapezoid(
                ((g[i]*nchi[j]+g[j]*nchi[i])*F/np.maximum(chi,1e-12))[None,:]*X,
                x=chi,axis=1
            )
            out["ii"][key]=np.trapezoid(
                (nchi[i]*nchi[j]*F**2/np.maximum(chi,1e-12)**2)[None,:]*P,
                x=chi,axis=1
            )

    meta={
        "model":m["name"],"variant":variant,"matter_tag":matter_tag,
        "weyl_mode":weyl_mode,
        "react_kmin_h_Mpc":KMIN,"react_kmax_h_Mpc":KMAX,
        "react_zmax":ZREACT_MAX,
        "high_k_policy":"freeze ReACT boost at kmax; no legacy CLASS HALOFIT tail",
        "high_z_policy":"taper nonlinear boost to unity between z=2.5 and z=3"
    }
    return out,meta

def main():
    zsrc,nzbins=source()
    hL=.6971482083084993; OmL=(.022083219194622913+.12299536722293603)/hL**2
    hC=.6856859; OmC=(.02240637+.11824151)/hC**2
    ms=[
        b.load_model("late300","output/weyl_spectra",hL,OmL),
        b.load_model("local021","output/weyl_spectra",hC,OmC)
    ]
    variants=[
        "linear",
        "us_native","ss_native",
        "us_screen03","ss_screen03",
        "us_fullscreen","ss_fullscreen"
    ]
    meta={}
    for v in variants:
        meta[v]={}
        for m in ms:
            x,mm=components(m,zsrc,nzbins,v)
            meta[v][m["name"]]=mm
            for comp,d in x.items():
                np.savez(OUT/f"{m['name']}_{v}_{comp}.npz",**d)
            print("KIDS_REACT_COMPONENT",m["name"],v,json.dumps(mm,sort_keys=True),flush=True)
    (OUT/"metadata.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")

if __name__=="__main__":
    main()
