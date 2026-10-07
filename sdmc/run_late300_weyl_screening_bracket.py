#!/usr/bin/env python3
"""
Nonlinear-response bracket for late300 DES-Y3 native-Weyl tomography.

The native hi_class Weyl output is linear.  For local021 we use the standard
HALOFIT matter boost.  For late300 we scan a screened nonlinear-response
exponent lambda_NL:
    Q_W,NL = Q_W,lin * B_HF^lambda_NL
    P_m,NL = P_m,lin * B_HF^lambda_NL
where B_HF=P_HF/P_lin from the same background and linear spectrum.

lambda_NL=0 is the strictly linear late300 response.
lambda_NL=1 is the naive unscreened HALOFIT envelope used in the first pilot.
Intermediate values quantify how much nonlinear enhancement can be tolerated
before late300 loses parity with local021.

This is a sensitivity bracket, not a derived SDMC nonlinear prescription.
"""
from pathlib import Path
import importlib.util, json, numpy as np, csv
from scipy.special import jv

spec=importlib.util.spec_from_file_location("weylbase","sdmc/run_late300_desy3_weyl_tomography.py")
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)

ROOT=Path("output/weyl_screening_bracket"); ROOT.mkdir(parents=True,exist_ok=True)
DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
C1RHO=base.C1RHO; ELL=base.ELL

h_l=69.71482083084993/100.0
Om_l=(0.022083219194622913+0.12299536722293603)/h_l**2
h_c=68.56859/100.0
Om_c=(0.02240637+0.11824151)/h_c**2
late=base.load_model("late300","output/weyl_spectra",h_l,Om_l)
lcdm=base.load_model("local021","output/weyl_spectra",h_c,Om_c)

official=base.load_des(DATA,CUTS,extra_conservative=False)

def subset(des,xip_min,xim_min):
    keep=[i for i,r in enumerate(des["rows"])
          if r["ang"] >= (xip_min if r["typ"]=="xip" else xim_min)]
    ix=np.asarray(keep,dtype=int)
    return dict(rows=[des["rows"][i] for i in keep],
                data=des["data"][ix],cov=des["cov"][np.ix_(ix,ix)],
                zsrc=des["zsrc"],nz=des["nz"],ndata_full=des.get("ndata_full",len(des["data"])))

def theory_lambda(model,des,lam):
    bg=model["bg"]; Om=model["Om"]
    zmax=min(model["zq"].max(),model["zp"].max(),float(np.max(des["zsrc"])))
    z,H,chi,nz,g,nchi=base.prepare_sources(des,bg,zmax)
    zmat=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    kmat=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
    Q=base.eval_cube(model["Iq"],zmat,kmat)
    Pl=base.eval_cube(model["Ip"],zmat,kmat)
    Pn=base.eval_cube(model["Ipnl"],zmat,kmat)
    valid=np.isfinite(Q)&np.isfinite(Pl)
    Q=np.where(valid,Q,0.0); Pl=np.where(valid,Pl,0.0)
    Pn=np.where(np.isfinite(Pn),Pn,Pl)
    boost=np.divide(Pn,Pl,out=np.ones_like(Pn),where=Pl>0)
    boost=np.clip(boost,1.0,50.0)
    fac=boost**float(lam)
    Q=Q*fac
    P=Pl*fac

    kref=0.05
    pref=np.exp(model["Ip"](np.column_stack([z,np.full_like(z,np.log(kref))])))
    p0=float(np.exp(model["Ip"]([[0.0,np.log(kref)]])[0]))
    D=np.maximum(np.sqrt(pref/p0),1e-4)
    Funit=-C1RHO*Om/D

    pair_cls={}
    cross=np.sqrt(np.maximum(Q*P,0.0))
    for i in range(4):
        for j in range(i,4):
            gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
            gi_int=((g[i]*nchi[j]+g[j]*nchi[i])*Funit/np.maximum(chi,1e-12))[None,:]*cross
            gi=np.trapezoid(gi_int,x=chi,axis=1)
            ii_int=(nchi[i]*nchi[j]*Funit**2/np.maximum(chi,1e-12)**2)[None,:]*P
            ii=np.trapezoid(ii_int,x=chi,axis=1)
            pair_cls[(i+1,j+1)]=(gg,gi,ii)

    t0=[];t1=[];t2=[]
    for r in des["rows"]:
        pair=(min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))
        gg,gi,ii=pair_cls[pair]
        th=np.deg2rad(r["ang"]/60.0)
        J=jv(0 if r["typ"]=="xip" else 4,ELL*th)
        wt=ELL/(2*np.pi)*J
        t0.append(np.trapezoid(wt*gg,x=ELL))
        t1.append(np.trapezoid(wt*gi,x=ELL))
        t2.append(np.trapezoid(wt*ii,x=ELL))
    return np.array(t0),np.array(t1),np.array(t2),dict(
        valid_grid_fraction=float(np.mean(valid)),
        boost_median=float(np.median(boost[valid])),
        boost_p95=float(np.percentile(boost[valid],95)))

cuts=[("official",0.0,0.0),("c20_100",20.0,100.0),("c40_150",40.0,150.0),
      ("c60_180",60.0,180.0),("c100_250",100.0,250.0)]
lams=np.round(np.linspace(0,1,21),3)

out={"status":"nonlinear-response sensitivity bracket; not a derived SDMC nonlinear prescription",
     "definition":"QW_NL=QW_lin*(P_HF/P_lin)^lambda_NL; local021 fixed at lambda_NL=1",
     "cuts":{}}
rows=[]

for label,xp,xm in cuts:
    des=subset(official,xp,xm)
    # fixed local021 standard-HALOFIT reference
    TL=theory_lambda(lcdm,des,1.0)
    pl=base.chi_profile(des,TL[:3])
    rec={"ndata":len(des["data"]),"local021_halofit":pl,"lambda_scan":[]}
    dvals=[]
    for lam in lams:
        T=theory_lambda(late,des,float(lam))
        pr=base.chi_profile(des,T[:3])
        delta=pr["chi2_profile"]-pl["chi2_profile"]
        row=dict(lambda_NL=float(lam),chi2_late300=pr["chi2_profile"],
                 AIA_late300=pr["AIA_best"],delta_chi2=delta,
                 valid_grid_fraction=T[3]["valid_grid_fraction"],
                 boost_median=T[3]["boost_median"],boost_p95=T[3]["boost_p95"])
        rec["lambda_scan"].append(row); dvals.append((float(lam),delta))
        rows.append([label,xp,xm,len(des["data"]),float(lam),pr["chi2_profile"],
                     pl["chi2_profile"],delta,pr["AIA_best"],pl["AIA_best"]])
    # interpolate first zero crossing if present
    cross=None
    for (l0,d0),(l1,d1) in zip(dvals[:-1],dvals[1:]):
        if d0==0 or d0*d1<0:
            cross=l0+(0-d0)*(l1-l0)/(d1-d0)
            break
    rec["lambda_parity_linear_interp"]=cross

    # Equal-treatment sensitivity: apply the same nonlinear-response exponent
    # to late300 and local021.  This isolates how the *relative* ranking moves
    # as one goes from linear to a common HALOFIT-shaped nonlinear correction.
    common=[]
    common_crossings=[]
    prev=None
    for late_row,lam in zip(rec["lambda_scan"],lams):
        TLc=theory_lambda(lcdm,des,float(lam))
        plc=base.chi_profile(des,TLc[:3])
        delta_common=late_row["chi2_late300"]-plc["chi2_profile"]
        rr=dict(lambda_NL=float(lam),
                chi2_late300=late_row["chi2_late300"],
                chi2_local021=plc["chi2_profile"],
                delta_chi2=delta_common,
                AIA_late300=late_row["AIA_late300"],
                AIA_local021=plc["AIA_best"])
        common.append(rr)
        if prev is not None and prev["delta_chi2"]*delta_common < 0:
            l0=prev["lambda_NL"]; d0=prev["delta_chi2"]
            lc=l0+(0-d0)*(float(lam)-l0)/(delta_common-d0)
            common_crossings.append(lc)
        prev=rr
    rec["common_lambda_scan"]=common
    rec["common_lambda_crossings_linear_interp"]=common_crossings
    out["cuts"][label]=rec
    print("WEYL_SCREENING_BRACKET",label,json.dumps(rec,sort_keys=True),flush=True)

with (ROOT/"late300_weyl_screening_bracket.csv").open("w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["cuts","xip_min","xim_min","ndata","lambda_NL","chi2_late300",
                "chi2_local021_halofit","delta_chi2","AIA_late300","AIA_local021"])
    w.writerows(rows)
(ROOT/"late300_weyl_screening_bracket.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
