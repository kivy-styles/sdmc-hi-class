#!/usr/bin/env python3
"""
DES-Y3 native-Weyl nuisance-profile audit for accepted-action late300.

This extends the previous tomography pilot by profiling the nuisance parameters
that dominate the DES-Y3 cosmic-shear calibration:
  * 4 source photo-z mean shifts with the public DES-Y3 Gaussian priors;
  * 4 multiplicative shear calibration parameters with the public priors;
  * NLA intrinsic-alignment amplitude A_IA and redshift slope eta_IA.

It compares accepted-action late300 and optimized local021 with exactly the
same nuisance model and the same nonlinear-response exponent lambda_NL.

The nonlinear response is still a sensitivity envelope,
    Q_W -> Q_W * (P_HF/P_lin)^lambda_NL,
not a derived SDMC nonlinear prescription.  The run therefore remains a
nuisance-complete *pilot* rather than the final official DES-Y3 likelihood.
"""
from pathlib import Path
import importlib.util, json, csv
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.special import jv

spec=importlib.util.spec_from_file_location("weylbase","sdmc/run_late300_desy3_weyl_tomography.py")
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)

OUT=Path("output/weyl_nuisance_profile"); OUT.mkdir(parents=True,exist_ok=True)
DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
ELL=base.ELL
C1RHO=base.C1RHO
ZPIV=0.62

DZ_SIG=np.array([0.018,0.015,0.011,0.017])
M_MU=np.array([-0.0063,-0.0198,-0.0241,-0.0369])
M_SIG=np.array([0.0091,0.0078,0.0076,0.0076])

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

def tail_integral(z,y):
    zr=z[::-1]; yr=y[::-1]
    c=np.concatenate([[0.0],np.cumsum((yr[1:]+yr[:-1])*0.5*np.diff(zr))])
    return (-c)[::-1]

def prepare_fixed(model,des,lam):
    zmax=min(model["zq"].max(),model["zp"].max(),float(np.max(des["zsrc"])))
    zlo=max(0.01,float(np.min(des["zsrc"][des["zsrc"]>0])))
    zhi=min(zmax,float(np.max(des["zsrc"])))
    z=np.linspace(zlo,zhi,320)
    H=np.interp(z,model["bg"]["z"],model["bg"]["H"])
    chi=np.interp(z,model["bg"]["z"],model["bg"]["chi"])
    zmat=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    kmat=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)

    Q=base.eval_cube(model["Iq"],zmat,kmat)
    Pl=base.eval_cube(model["Ip"],zmat,kmat)
    Pn=base.eval_cube(model["Ipnl"],zmat,kmat)
    valid=np.isfinite(Q)&np.isfinite(Pl)
    Q=np.where(valid,Q,0.0)
    Pl=np.where(valid,Pl,0.0)
    Pn=np.where(np.isfinite(Pn),Pn,Pl)
    boost=np.divide(Pn,Pl,out=np.ones_like(Pn),where=Pl>0)
    boost=np.clip(boost,1.0,50.0)
    fac=boost**float(lam)
    Q=Q*fac
    P=Pl*fac
    cross=np.sqrt(np.maximum(Q*P,0.0))

    kref=0.05
    pref=np.exp(model["Ip"](np.column_stack([z,np.full_like(z,np.log(kref))])))
    p0=float(np.exp(model["Ip"]([[0.0,np.log(kref)]])[0]))
    D=np.maximum(np.sqrt(pref/p0),1e-4)
    Fbase=-C1RHO*model["Om"]/D

    # Trapezoid weights for chi and ell.
    wx=np.empty_like(chi)
    wx[1:-1]=(chi[2:]-chi[:-2])*0.5
    wx[0]=(chi[1]-chi[0])*0.5; wx[-1]=(chi[-1]-chi[-2])*0.5
    wl=np.empty_like(ELL)
    wl[1:-1]=(ELL[2:]-ELL[:-2])*0.5
    wl[0]=(ELL[1]-ELL[0])*0.5; wl[-1]=(ELL[-1]-ELL[-2])*0.5

    # Bessel transform weights row-by-row.
    Brow=[]
    for r in des["rows"]:
        th=np.deg2rad(r["ang"]/60.0)
        J=jv(0 if r["typ"]=="xip" else 4,ELL*th)
        Brow.append(wl*ELL/(2*np.pi)*J)
    Brow=np.asarray(Brow)

    L=np.linalg.cholesky(des["cov"])
    return dict(z=z,H=H,chi=chi,Q=Q,P=P,cross=cross,D=D,Fbase=Fbase,
                wx=wx,Brow=Brow,L=L,des=des,model=model,lam=float(lam),
                valid_grid_fraction=float(np.mean(valid)),
                boost_p95=float(np.percentile(boost[valid],95)))

def shifted_nz(des,z,dz):
    out=[]
    for i,arr in enumerate(des["nz"]):
        f=interp1d(des["zsrc"],arr,bounds_error=False,fill_value=0.0,kind="linear")
        # CosmoSIS additive photo-z convention: n_biased(z)=n(z-delta_z)
        y=np.maximum(f(z-dz[i]),0.0)
        n=np.trapezoid(y,z)
        if not np.isfinite(n) or n<=0: return None
        out.append(y/n)
    return np.asarray(out)

def theory_from_nuis(fix,x):
    # x = dz1..4, m1..4, AIA, etaIA
    dz=np.asarray(x[:4]); m=np.asarray(x[4:8]); A=float(x[8]); eta=float(x[9])
    z=fix["z"]; H=fix["H"]; chi=fix["chi"]; des=fix["des"]
    nz=shifted_nz(des,z,dz)
    if nz is None: return None
    g=[]
    for y in nz:
        t0=tail_integral(z,y)
        t1=tail_integral(z,y/np.maximum(chi,1e-12))
        g.append(t0-chi*t1)
    g=np.asarray(g)
    nchi=nz*H[None,:]

    Fred=fix["Fbase"]*((1+z)/(1+ZPIV))**eta
    Q=fix["Q"]; P=fix["P"]; cross=fix["cross"]; wx=fix["wx"]
    cls={}
    for i in range(4):
        for j in range(i,4):
            gg=np.sum(g[i][None,:]*g[j][None,:]*Q*wx[None,:],axis=1)
            gi0=((g[i]*nchi[j]+g[j]*nchi[i])*Fred/np.maximum(chi,1e-12))[None,:]*cross
            gi=np.sum(gi0*wx[None,:],axis=1)
            ii0=(nchi[i]*nchi[j]*Fred**2/np.maximum(chi,1e-12)**2)[None,:]*P
            ii=np.sum(ii0*wx[None,:],axis=1)
            cls[(i+1,j+1)]=(gg,A*gi,A*A*ii)

    th=np.empty(len(des["rows"]))
    for q,r in enumerate(des["rows"]):
        pair=(min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))
        gg,gi,ii=cls[pair]
        val=float(np.dot(fix["Brow"][q],gg+gi+ii))
        val*= (1+m[pair[0]-1])*(1+m[pair[1]-1])
        th[q]=val
    return th

def objective(fix,x,return_parts=False):
    th=theory_from_nuis(fix,x)
    if th is None or np.any(~np.isfinite(th)):
        return 1e30 if not return_parts else None
    y=np.linalg.solve(fix["L"],fix["des"]["data"]-th)
    data=float(y@y)
    dz=np.asarray(x[:4]); m=np.asarray(x[4:8])
    prior=float(np.sum((dz/DZ_SIG)**2)+np.sum(((m-M_MU)/M_SIG)**2))
    total=data+prior
    if return_parts:
        return dict(total=total,data_chi2=data,prior_chi2=prior,theory=th)
    return total

BOUNDS=[(-0.09,0.09),(-0.075,0.075),(-0.055,0.055),(-0.085,0.085),
        (-0.10,0.10),(-0.10,0.10),(-0.10,0.10),(-0.10,0.10),
        (-5.0,5.0),(-5.0,5.0)]

def fit_model(model,des,lam):
    fix=prepare_fixed(model,des,lam)
    starts=[]
    x=np.zeros(10); x[4:8]=M_MU; x[8]=0.2; x[9]=0.0
    starts.append(x)
    best=None
    for x0 in starts:
        r=minimize(lambda v: objective(fix,v),x0,method="L-BFGS-B",bounds=BOUNDS,
                   options=dict(maxiter=55,ftol=1e-8,gtol=3e-5,maxls=20))
        parts=objective(fix,r.x,return_parts=True)
        rec=dict(success=bool(r.success),message=str(r.message),nit=int(r.nit),nfev=int(r.nfev),
                 total_chi2=float(parts["total"]),data_chi2=float(parts["data_chi2"]),
                 prior_chi2=float(parts["prior_chi2"]),params=[float(v) for v in r.x],
                 valid_grid_fraction=fix["valid_grid_fraction"],boost_p95=fix["boost_p95"])
        if best is None or rec["total_chi2"]<best["total_chi2"]: best=rec
    return best

cuts=[("official",0.0,0.0),("c100_250",100.0,250.0)]
lams=[0.55]
out={"status":"DES-Y3 nuisance-profile native-Weyl pilot",
     "nuisance_order":["dz1","dz2","dz3","dz4","m1","m2","m3","m4","AIA","etaIA"],
     "priors":{"dz_sigma":DZ_SIG.tolist(),"m_mean":M_MU.tolist(),"m_sigma":M_SIG.tolist(),
               "AIA":"top-hat [-5,5]","etaIA":"top-hat [-5,5]"},
     "nonlinear_note":"common lambda_NL sensitivity envelope, not a derived SDMC nonlinear model",
     "cuts":{}}
rows=[]
for label,xp,xm in cuts:
    des=subset(official,xp,xm)
    rr={"ndata":len(des["data"]),"xip_min":xp,"xim_min":xm,"lambda":{}}
    for lam in lams:
        a=fit_model(late,des,lam)
        b=fit_model(lcdm,des,lam)
        delta=a["total_chi2"]-b["total_chi2"]
        dd=a["data_chi2"]-b["data_chi2"]
        rr["lambda"][str(lam)]=dict(late300=a,local021=b,
                                      delta_profile_total_chi2=float(delta),
                                      delta_data_chi2=float(dd))
        rows.append([label,len(des["data"]),lam,a["total_chi2"],b["total_chi2"],delta,
                     a["data_chi2"],b["data_chi2"],dd,a["prior_chi2"],b["prior_chi2"],
                     a["params"][8],b["params"][8],a["params"][9],b["params"][9]])
        print("WEYL_NUISANCE_FASTCHECK",label,lam,json.dumps(rr["lambda"][str(lam)],sort_keys=True),flush=True)
    out["cuts"][label]=rr

with (OUT/"late300_weyl_nuisance_fastcheck.csv").open("w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["cuts","ndata","lambda_NL","late300_total","local021_total","delta_total",
                "late300_data","local021_data","delta_data","late300_prior","local021_prior",
                "late300_AIA","local021_AIA","late300_etaIA","local021_etaIA"])
    w.writerows(rows)
(OUT/"late300_weyl_nuisance_fastcheck.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("WEYL_NUISANCE_FASTCHECK_DONE",json.dumps(out,sort_keys=True))
