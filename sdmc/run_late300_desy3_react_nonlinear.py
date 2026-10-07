#!/usr/bin/env python3
from pathlib import Path
import importlib.util, json, csv
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import jv

S=importlib.util.spec_from_file_location("base","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

REACT=Path("output/react_dense/react_late300/late300_react_target_spectrum.npz")
OUT=Path("output/desy3_react"); OUT.mkdir(parents=True,exist_ok=True)
ZREACT_MAX=2.5

def boost_interp(zr,kr,arr):
    return RegularGridInterpolator((zr,np.log(kr)),np.log(np.maximum(arr,1e-12)),
                                   bounds_error=False,fill_value=None)

R=np.load(REACT)
ZR=np.asarray(R["z"],float); KR=np.asarray(R["k"],float)
KMIN=float(KR.min()); KMAX=float(KR.max())
BOOST={
    ("late300","us"):boost_interp(ZR,KR,R["Pnl_us"]/R["Ptarget"]),
    ("late300","ss"):boost_interp(ZR,KR,R["Pnl_ss"]/R["Ptarget"]),
    ("local021","gr"):boost_interp(ZR,KR,R["Pnl_local021"]/R["Ptarget_local021"]),
}

def eval_boost(I,z,kh):
    z=np.asarray(z,float); kh=np.asarray(kh,float)
    zc=np.clip(z,0,ZREACT_MAX)
    pts=np.column_stack([zc.ravel(),np.log(np.clip(kh.ravel(),KMIN,KMAX))])
    q=np.exp(I(pts)).reshape(z.shape)
    w=np.clip((3.0-z)/0.5,0.0,1.0)
    return 1.0+w*(q-1.0)

def model_theory(model,des,variant):
    bg=model["bg"]; h=model["h"]; Om=model["Om"]
    zmax=min(model["zq"].max(),model["zp"].max(),float(np.max(des["zsrc"])))
    z,H,chi,nz,g,nchi=b.prepare_sources(des,bg,zmax)
    zmat=np.broadcast_to(z[None,:],(len(b.ELL),len(z)))
    kmat=(b.ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
    kh=kmat/h

    Q=b.eval_cube(model["Iq"],zmat,kmat)
    Pl=b.eval_cube(model["Ip"],zmat,kmat)
    valid=np.isfinite(Q)&np.isfinite(Pl)&(Pl>0)
    Q=np.where(valid,Q,0.0); Pl=np.where(valid,Pl,0.0)

    if variant=="linear":
        P=Pl; Quse=Q
        matter_tag="linear"; weyl_mode="linear"
    else:
        if model["name"]=="late300":
            matter_tag="ss" if variant.startswith("ss_") else "us"
            I=BOOST[("late300",matter_tag)]
        else:
            matter_tag="gr"; I=BOOST[("local021","gr")]
        boost=eval_boost(I,zmat,kh)
        P=Pl*boost

        H0=100*h/b.C_KMS
        Agr=(1.5*Om*H0**2*(1+z))[None,:]
        sig2=np.divide(Q,Agr**2*Pl,out=np.ones_like(Q),where=(Agr**2*Pl)>0)
        sig=np.sqrt(np.clip(sig2,0.2,2.0))
        if variant.endswith("_native"):
            sigeff=sig; weyl_mode="native_linear_Sigma"
        elif variant.endswith("_screen03"):
            x=kh/0.3
            sigeff=1.0+(sig-1.0)/(1.0+x*x)
            weyl_mode="No-Slip transition k_s=0.3 h/Mpc"
        elif variant.endswith("_fullscreen"):
            sigeff=np.ones_like(sig); weyl_mode="full nonlinear Weyl screen to GR Sigma=1"
        else:
            raise ValueError(variant)
        Quse=Agr**2*P*sigeff**2
        Quse=np.where(valid,Quse,0.0)

    kref=0.05
    pref=np.exp(model["Ip"](np.column_stack([z,np.full_like(z,np.log(kref))])))
    p0=float(np.exp(model["Ip"]([[0.0,np.log(kref)]])[0]))
    D=np.maximum(np.sqrt(pref/p0),1e-4)
    Funit=-b.C1RHO*Om/D

    pair_cls={}
    cross=np.sqrt(np.maximum(Quse*P,0.0))
    for i in range(4):
        for j in range(i,4):
            gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Quse,x=chi,axis=1)
            gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*Funit/np.maximum(chi,1e-12))[None,:]*cross,x=chi,axis=1)
            ii=np.trapezoid((nchi[i]*nchi[j]*Funit**2/np.maximum(chi,1e-12)**2)[None,:]*P,x=chi,axis=1)
            pair_cls[(i+1,j+1)]=(gg,gi,ii)

    t0=[];t1=[];t2=[]
    for r in des["rows"]:
        pair=(min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))
        gg,gi,ii=pair_cls[pair]
        th=np.deg2rad(r["ang"]/60.0)
        J=jv(0 if r["typ"]=="xip" else 4,b.ELL*th)
        wt=b.ELL/(2*np.pi)*J
        t0.append(np.trapezoid(wt*gg,x=b.ELL))
        t1.append(np.trapezoid(wt*gi,x=b.ELL))
        t2.append(np.trapezoid(wt*ii,x=b.ELL))
    return np.array(t0),np.array(t1),np.array(t2),{
        "matter_tag":matter_tag,"weyl_mode":weyl_mode,
        "react_kmin_h_Mpc":KMIN,"react_kmax_h_Mpc":KMAX,"react_zmax":ZREACT_MAX
    }

def main():
    data=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
    cuts=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
    hL=.6971482083084993; OmL=(.022083219194622913+.12299536722293603)/hL**2
    hC=.6856859; OmC=(.02240637+.11824151)/hC**2
    late=b.load_model("late300","output/weyl_spectra",hL,OmL)
    lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)
    variants=["linear","us_native","ss_native","us_screen03","ss_screen03","us_fullscreen","ss_fullscreen"]
    out={
      "status":"DES-Y3 xi+/xi- ReACT nonlinear native-Weyl audit; gamma2=gamma3=0 ReACT bracket; one-parameter NLA A_IA profile",
      "nonlinear_limits":{"ReACT_zmax":2.5,"high_z":"boost tapered to unity from z=2.5 to z=3","high_k":"boost frozen at validated kmax"},
      "results":{}
    }
    for label,extra in [("official_cuts",False),("large_scale_cuts",True)]:
        des=b.load_des(data,cuts,extra_conservative=extra)
        rec={"ndata":len(des["data"]),"variants":{}}
        for v in variants:
            rr={}
            for m in [late,lcdm]:
                T=model_theory(m,des,v)
                p=b.chi_profile(des,T[:3]); p["projection"]=T[3]
                rr[m["name"]]=p
            rr["delta_chi2_profile"]=rr["late300"]["chi2_profile"]-rr["local021"]["chi2_profile"]
            rr["delta_chi2_A0"]=rr["late300"]["chi2_A0"]-rr["local021"]["chi2_A0"]
            rec["variants"][v]=rr
        out["results"][label]=rec
    (OUT/"desy3_react_profile.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    with (OUT/"desy3_react_profile.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["cuts","variant","model","ndata","chi2_A0","chi2_profile","AIA_best","delta_chi2_profile"])
        for cut,rec in out["results"].items():
            for v,r in rec["variants"].items():
                for m in ["late300","local021"]:
                    q=r[m]; w.writerow([cut,v,m,rec["ndata"],q["chi2_A0"],q["chi2_profile"],q["AIA_best"],r["delta_chi2_profile"]])
    print("DESY3_REACT_PROFILE",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__":
    main()
