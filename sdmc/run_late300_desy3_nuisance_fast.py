#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,csv,numpy as np
from scipy.optimize import minimize
from scipy.special import jv

S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)
b.ELL=np.geomspace(2.0,1.0e4,300)

DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
OUT=Path("output/weyl_nuisance_fast"); OUT.mkdir(parents=True,exist_ok=True)

hL=0.6971482083084993
OmL=(0.022083219194622913+0.12299536722293603)/hL**2
hC=0.6856859
OmC=(0.02240637+0.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)
official=b.load_des(DATA,CUTS,False)

DZ_SIG=np.array([0.018,0.015,0.011,0.017])
M_MEAN=np.array([-0.0063,-0.0198,-0.0241,-0.0369])
M_SIG=np.array([0.0091,0.0078,0.0076,0.0076])
PIVOT=0.62

def subset(des,xp,xm):
    keep=[i for i,r in enumerate(des["rows"]) if r["ang"] >= (xp if r["typ"]=="xip" else xm)]
    ix=np.asarray(keep,int)
    return dict(rows=[des["rows"][i] for i in keep],data=des["data"][ix],
                cov=des["cov"][np.ix_(ix,ix)],zsrc=des["zsrc"],nz=des["nz"])

def prep_model(m,des):
    zmax=min(m["zq"].max(),m["zp"].max(),float(np.max(des["zsrc"])))
    zlo=max(0.01,float(np.min(des["zsrc"][des["zsrc"]>0])))
    zhi=min(zmax,float(np.max(des["zsrc"])))
    z=np.linspace(zlo,zhi,180)
    H=np.interp(z,m["bg"]["z"],m["bg"]["H"])
    chi=np.interp(z,m["bg"]["z"],m["bg"]["chi"])
    zm=np.broadcast_to(z[None,:],(len(b.ELL),len(z)))
    km=(b.ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
    Q=b.eval_cube(m["Iq"],zm,km); P=b.eval_cube(m["Ip"],zm,km)
    ok=np.isfinite(Q)&np.isfinite(P)
    Q=np.where(ok,Q,0.0); P=np.where(ok,P,0.0)
    pref=np.exp(m["Ip"](np.column_stack([z,np.full_like(z,np.log(0.05))])))
    p0=float(np.exp(m["Ip"]([[0.0,np.log(0.05)]])[0]))
    D=np.maximum(np.sqrt(pref/p0),1e-4)
    F0=-b.C1RHO*m["Om"]/D
    L=np.linalg.cholesky(des["cov"])
    return dict(z=z,H=H,chi=chi,Q=Q,P=P,X=np.sqrt(np.maximum(Q*P,0.0)),
                F0=F0,L=L,des=des,m=m,valid=float(np.mean(ok)))

def shifted_sources(ctx,dz):
    z=ctx["z"]; chi=ctx["chi"]; H=ctx["H"]; des=ctx["des"]
    nz=[]
    for i,arr in enumerate(des["nz"]):
        y=np.interp(z+dz[i],des["zsrc"],arr,left=0.0,right=0.0)
        y=np.maximum(y,0.0)
        n=np.trapezoid(y,z)
        if n<=0: return None
        nz.append(y/n)
    nz=np.asarray(nz)
    g=[]
    for y in nz:
        t0=b.tail_integral(z,y)
        t1=b.tail_integral(z,y/np.maximum(chi,1e-12))
        g.append(t0-chi*t1)
    return nz,np.asarray(g),nz*H[None,:]

def theory(ctx,p):
    A,eta=p[0],p[1]; dz=np.asarray(p[2:6]); mm=np.asarray(p[6:10])
    src=shifted_sources(ctx,dz)
    if src is None: return None
    nz,g,nchi=src
    z,chi,Q,P,X,F0=ctx["z"],ctx["chi"],ctx["Q"],ctx["P"],ctx["X"],ctx["F0"]
    F=F0*((1+z)/(1+PIVOT))**eta
    pair={}
    for i in range(4):
      for j in range(i,4):
        gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
        gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*(A*F)/np.maximum(chi,1e-12))[None,:]*X,x=chi,axis=1)
        ii=np.trapezoid((nchi[i]*nchi[j]*(A*F)**2/np.maximum(chi,1e-12)**2)[None,:]*P,x=chi,axis=1)
        pair[(i+1,j+1)]=(gg+gi+ii)*(1+mm[i])*(1+mm[j])
    th=[]
    for r in ctx["des"]["rows"]:
      cl=pair[(min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))]
      ang=np.deg2rad(r["ang"]/60.0)
      wt=b.ELL/(2*np.pi)*jv(0 if r["typ"]=="xip" else 4,b.ELL*ang)
      th.append(np.trapezoid(wt*cl,x=b.ELL))
    return np.asarray(th)

def objective(ctx,p):
    if not (-5<=p[0]<=5 and -5<=p[1]<=5): return 1e30
    th=theory(ctx,p)
    if th is None or np.any(~np.isfinite(th)): return 1e30
    y=np.linalg.solve(ctx["L"],ctx["des"]["data"]-th)
    chi=float(y@y)
    chi+=float(np.sum((np.asarray(p[2:6])/DZ_SIG)**2))
    chi+=float(np.sum(((np.asarray(p[6:10])-M_MEAN)/M_SIG)**2))
    return chi

def fit(ctx):
    x0=np.r_[0.0,0.0,np.zeros(4),M_MEAN]
    bounds=[(-5,5),(-5,5)]+[(-0.08,0.08)]*4+[(-0.08,0.04)]*4
    starts=[x0,x0.copy(),x0.copy()]
    starts[1][0]=1.0; starts[1][1]=1.0
    starts[2][0]=-1.0; starts[2][1]=-1.0
    best=None
    for x in starts:
        r=minimize(lambda q:objective(ctx,q),x,method="L-BFGS-B",bounds=bounds,
                   options={"maxiter":35,"ftol":1e-7,"gtol":3e-5})
        if best is None or r.fun<best.fun: best=r
    return dict(chi2=float(best.fun),success=bool(best.success),message=str(best.message),
                AIA=float(best.x[0]),etaIA=float(best.x[1]),
                dz=[float(v) for v in best.x[2:6]],m=[float(v) for v in best.x[6:10]])

cuts={"c40_150":(40.0,150.0),"c60_180":(60.0,180.0),"c100_250":(100.0,250.0)}
out={"status":"linear native-Weyl DES-Y3 NLA nuisance profile with published dz/m Gaussian priors",
     "priors":{"dz_sigma":DZ_SIG.tolist(),"m_mean":M_MEAN.tolist(),"m_sigma":M_SIG.tolist(),
               "AIA":[-5,5],"etaIA":[-5,5]},"cuts":{}}
rows=[]
for label,(xp,xm) in cuts.items():
    des=subset(official,xp,xm)
    rec={"ndata":len(des["data"]),"xip_min":xp,"xim_min":xm}
    for m in [late,lcdm]:
        ctx=prep_model(m,des); fitr=fit(ctx); fitr["valid_grid_fraction"]=ctx["valid"]
        rec[m["name"]]=fitr
    rec["delta_chi2_late300_minus_local021"]=rec["late300"]["chi2"]-rec["local021"]["chi2"]
    out["cuts"][label]=rec
    rows.append([label,len(des["data"]),rec["late300"]["chi2"],rec["local021"]["chi2"],
                 rec["delta_chi2_late300_minus_local021"],rec["late300"]["AIA"],rec["local021"]["AIA"],
                 rec["late300"]["etaIA"],rec["local021"]["etaIA"]])
    print("WEYL_NUISANCE_FAST",label,json.dumps(rec,sort_keys=True),flush=True)

(OUT/"late300_desy3_nuisance.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
with (OUT/"late300_desy3_nuisance.csv").open("w",newline="") as f:
    w=csv.writer(f); w.writerow(["cuts","ndata","chi2_late300","chi2_local021","delta_chi2",
                                 "AIA_late300","AIA_local021","etaIA_late300","etaIA_local021"]); w.writerows(rows)
