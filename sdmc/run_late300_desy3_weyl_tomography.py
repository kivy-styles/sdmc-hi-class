#!/usr/bin/env python3
"""
late300 DES-Y3 tomographic cosmic-shear pilot using native hi_class Weyl power.

This is deliberately stricter than the earlier compressed S8 screen:
  * loads the public DES-Y3 xi+/xi- data vector, covariance and 4 source n(z);
  * projects the native hi_class Q_W=k^4 P_{(Phi+Psi)/2} output;
  * uses the same pipeline for accepted-action late300 and optimized local021;
  * validates the local021 GR normalization against the matter-power Poisson relation;
  * profiles a one-parameter NLA intrinsic-alignment amplitude A_IA;
  * reports both pure-linear Weyl and a HALOFIT-envelope variant.

The HALOFIT-envelope variant is NOT an exact nonlinear modified-gravity prediction:
Q_W,NL is approximated as Q_W,lin * (P_m,NL/P_m,lin) model by model.
Photo-z shifts, per-bin shear calibration, baryonic-feedback nuisances and eta_IA
are not yet profiled.  Results are therefore a tomographic likelihood pilot,
not the final official DES-Y3 likelihood.
"""
from pathlib import Path
import configparser, json, math, re
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar, minimize
from scipy.special import jv
from astropy.io import fits

C_KMS=299792.458
C1RHO=0.0134
ELL=np.geomspace(2.0,1.0e4,720)

def parse_background(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    a=np.loadtxt(path)
    d={n:a[:,i] for i,n in enumerate(names)}
    z=np.asarray(d["z"])
    hname="H [1/Mpc]"
    cname=next((n for n in names if "comov" in n.lower() and "dist" in n.lower()),None)
    if cname is None:
        raise RuntimeError(f"no comoving-distance column in {names}")
    o=np.argsort(z)
    return dict(z=z[o],H=np.asarray(d[hname])[o],chi=np.asarray(d[cname])[o],raw=d,names=names)

def z_from_header(path):
    for line in Path(path).read_text().splitlines()[:8]:
        m=re.search(r"(?:at\s+)?(?:redshift\s+)?z\s*=\s*([0-9eE+\-.]+)",line,re.I)
        if m: return float(m.group(1))
    raise RuntimeError(f"cannot read redshift from {path}")

def load_spectral_cube(directory, prefix, h, power_kind):
    if power_kind=="weyl":
        files=sorted(Path(directory).glob(f"{prefix}*pk_weyl.dat"))
    else:
        files=sorted(Path(directory).glob(f"{prefix}*pk.dat"))
    rows=[]
    for p in files:
        if p.name.endswith("_ad.dat") or "_ad_" in p.name:
            continue
        z=z_from_header(p)
        a=np.loadtxt(p)
        k=a[:,0]*h
        if power_kind=="weyl":
            y=a[:,1]*h              # h/Mpc -> 1/Mpc
        elif power_kind=="matter":
            y=a[:,1]/h**3           # (Mpc/h)^3 -> Mpc^3
        else:
            raise ValueError(power_kind)
        if np.any(y<=0) or np.any(~np.isfinite(y)):
            raise RuntimeError(f"nonpositive/nonfinite {power_kind} spectrum in {p}")
        rows.append((z,k,y))
    if not rows:
        raise RuntimeError(f"no files for {prefix} / {power_kind} in {directory}")
    rows.sort(key=lambda x:x[0])
    z=np.array([r[0] for r in rows])
    k0=rows[0][1]
    vals=[]
    for zz,k,y in rows:
        if len(k)!=len(k0) or not np.allclose(k,k0,rtol=1e-8,atol=0):
            y=np.exp(np.interp(np.log(k0),np.log(k),np.log(y)))
        vals.append(y)
    return z,k0,np.asarray(vals)

def cube_interp(z,k,val):
    return RegularGridInterpolator((z,np.log(k)),np.log(val),bounds_error=False,fill_value=np.nan)

def eval_cube(interp,zmat,kmat):
    pts=np.column_stack([zmat.ravel(),np.log(kmat.ravel())])
    out=np.exp(interp(pts)).reshape(kmat.shape)
    return out

def tail_integral(z,y):
    zr=z[::-1]; yr=y[::-1]
    cum=cumulative_trapezoid(yr,zr,initial=0.0)
    return (-cum)[::-1]

def read_scale_cuts(path):
    cp=configparser.ConfigParser()
    cp.read(path)
    sec=cp["2pt_like"]
    cuts={}
    for key,val in sec.items():
        m=re.match(r"angle_range_(xip|xim)_(\d+)_(\d+)",key)
        if not m: continue
        cuts[(m.group(1),int(m.group(2)),int(m.group(3)))]=float(val.split()[0])
    return cuts

def load_des(path,scale_path,extra_min=None):
    with fits.open(path) as f:
        xp=f["xip"].data; xm=f["xim"].data
        nz=f["nz_source"].data
        full_cov=np.asarray(f["COVMAT"].data,float)
        zsrc=np.asarray(nz["Z_MID"],float)
        nzbins=np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,5)])
        rows=[]
        for typ,dat,offset in [("xip",xp,0),("xim",xm,len(xp))]:
            for q,r in enumerate(dat):
                rows.append(dict(typ=typ,ang=float(r["ANG"]),b1=int(r["BIN1"]),b2=int(r["BIN2"]),
                                 val=float(r["VALUE"]),idx=offset+q))
    nfull=len(xp)+len(xm)
    if full_cov.shape[0] < nfull:
        raise RuntimeError(f"covariance shape {full_cov.shape} shorter than xi block {nfull}")
    cuts=read_scale_cuts(scale_path)
    keep=[]
    for r in rows:
        key=(r["typ"],min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))
        amin=cuts.get(key,0.0)
        if extra_min is not None:
            amin=max(amin,float(extra_min[0] if r["typ"]=="xip" else extra_min[1]))
        if r["ang"]>=amin:
            keep.append(r)
    idx=np.array([r["idx"] for r in keep],int)
    data=np.array([r["val"] for r in keep])
    cov=full_cov[np.ix_(idx,idx)]
    return dict(rows=keep,data=data,cov=cov,zsrc=zsrc,nz=nzbins,ndata_full=nfull)

def prepare_sources(des,bg,zmax,dz_shift=None):
    zlo=max(0.01,float(np.min(des["zsrc"][des["zsrc"]>0])) if np.any(des["zsrc"]>0) else 0.01)
    zhi=min(zmax,float(np.max(des["zsrc"])))
    z=np.linspace(zlo,zhi,360)
    H=np.interp(z,bg["z"],bg["H"])
    chi=np.interp(z,bg["z"],bg["chi"])
    if dz_shift is None:
        dz_shift=np.zeros(4)
    dz_shift=np.asarray(dz_shift,float)
    nz=[]
    for ib,arr in enumerate(des["nz"]):
        f=interp1d(des["zsrc"],arr,bounds_error=False,fill_value=0.0)
        # Match the CosmoSIS DES-Y3 additive source-photo-z convention:
        # n_biased(z) = n_input(z - bias_i)
        y=np.maximum(f(z-dz_shift[ib]),0.0)
        norm=np.trapezoid(y,z)
        y=y/norm
        nz.append(y)
    nz=np.asarray(nz)
    g=[]
    for y in nz:
        t0=tail_integral(z,y)
        t1=tail_integral(z,y/np.maximum(chi,1e-12))
        g.append(t0-chi*t1)
    g=np.asarray(g)
    nchi=nz*H[None,:]  # n(chi)=n(z) dz/dchi = n(z) H(z), c=1
    return z,H,chi,nz,g,nchi

def model_theory(model,des,variant,shared_boost_model=None,dz_shift=None):
    bg=model["bg"]; h=model["h"]; Om=model["Om"]
    zmax=min(model["zq"].max(),model["zp"].max(),float(np.max(des["zsrc"])))
    z,H,chi,nz,g,nchi=prepare_sources(des,bg,zmax,dz_shift=dz_shift)
    zmat=np.broadcast_to(z[None,:],(len(ELL),len(z)))
    kmat=(ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)

    Q=eval_cube(model["Iq"],zmat,kmat)
    Pl=eval_cube(model["Ip"],zmat,kmat)
    valid=np.isfinite(Q)&np.isfinite(Pl)
    Q=np.where(valid,Q,0.0); Pl=np.where(valid,Pl,0.0)

    if variant=="halofit_envelope":
        Pn=eval_cube(model["Ipnl"],zmat,kmat)
        Pn=np.where(np.isfinite(Pn),Pn,Pl)
        boost=np.divide(Pn,Pl,out=np.ones_like(Pn),where=Pl>0)
        boost=np.clip(boost,0.2,50.0)
        Q=Q*boost
        P=Pl*boost
    elif variant=="shared_lcdm_halofit":
        if shared_boost_model is None:
            raise RuntimeError("shared_lcdm_halofit requires a reference model")
        Pref=eval_cube(shared_boost_model["Ip"],zmat,kmat)
        Pnref=eval_cube(shared_boost_model["Ipnl"],zmat,kmat)
        boost=np.divide(Pnref,Pref,out=np.ones_like(Pnref),where=np.isfinite(Pref)&(Pref>0))
        boost=np.where(np.isfinite(boost),boost,1.0)
        boost=np.clip(boost,0.2,50.0)
        Q=Q*boost
        P=Pl*boost
    else:
        P=Pl

    # NLA growth normalization from a safely linear reference scale.
    kref=0.05
    zref=z
    pref=np.exp(model["Ip"](np.column_stack([zref,np.full_like(zref,np.log(kref))])))
    p0=float(np.exp(model["Ip"]([[0.0,np.log(kref)]])[0]))
    D=np.sqrt(pref/p0)
    D=np.maximum(D,1e-4)
    Funit=-C1RHO*Om/D  # A_IA multiplies this; eta_IA fixed to zero

    # precompute C_ell components for all ten tomographic pairs
    pair_cls={}
    cross=np.sqrt(np.maximum(Q*P,0.0))
    for i in range(4):
        for j in range(i,4):
            gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
            gi_integrand=((g[i]*nchi[j]+g[j]*nchi[i])*Funit/np.maximum(chi,1e-12))[None,:]*cross
            gi=np.trapezoid(gi_integrand,x=chi,axis=1)
            ii_integrand=(nchi[i]*nchi[j]*Funit**2/np.maximum(chi,1e-12)**2)[None,:]*P
            ii=np.trapezoid(ii_integrand,x=chi,axis=1)
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
        zmin=float(z.min()),zmax=float(z.max()),
        valid_grid_fraction=float(np.mean(valid)),
        Sigma_note="Native Q_W projection; no separate Sigma approximation is used."
    )

def chi_profile(des,T):
    d=des["data"]; C=des["cov"]
    # stable inverse action
    L=np.linalg.cholesky(C)
    def chi(A):
        th=T[0]+A*T[1]+A*A*T[2]
        y=np.linalg.solve(L,d-th)
        return float(y@y)
    res=minimize_scalar(chi,bounds=(-5.0,5.0),method="bounded",options=dict(xatol=1e-4))
    return dict(chi2_A0=chi(0.0),chi2_profile=float(res.fun),AIA_best=float(res.x),success=bool(res.success))


def apply_mcal(des,theory,m):
    m=np.asarray(m,float)
    fac=np.array([(1.0+m[r["b1"]-1])*(1.0+m[r["b2"]-1]) for r in des["rows"]])
    return theory*fac

def data_chi2(des,theory):
    L=np.linalg.cholesky(des["cov"])
    y=np.linalg.solve(L,des["data"]-theory)
    return float(y@y)

def chi_profile_des_nuisance(des,model,variant,shared_boost_model):
    """
    Profile the eight Gaussian DES-Y3 source calibration nuisances plus A_IA.
    Source photo-z response is linearized around the public central n(z), then
    validated once with an exact shifted-kernel recomputation at the optimum.
    eta_IA remains fixed at zero and TATT/baryon nuisance sectors remain open.
    """
    dz_sig=np.array([0.018,0.015,0.011,0.017])
    m_mu=np.array([-0.0063,-0.0198,-0.0241,-0.0369])
    m_sig=np.array([0.0091,0.0078,0.0076,0.0076])
    base=model_theory(model,des,variant,shared_boost_model=shared_boost_model,dz_shift=np.zeros(4))
    B=np.vstack(base[:3])  # component, data
    step=0.005
    deriv=[]
    for ib in range(4):
        dp=np.zeros(4); dm=np.zeros(4); dp[ib]=step; dm[ib]=-step
        Tp=model_theory(model,des,variant,shared_boost_model=shared_boost_model,dz_shift=dp)
        Tm=model_theory(model,des,variant,shared_boost_model=shared_boost_model,dz_shift=dm)
        deriv.append((np.vstack(Tp[:3])-np.vstack(Tm[:3]))/(2*step))
    deriv=np.asarray(deriv) # bin, component, data

    simple=chi_profile(des,base[:3])
    x0=np.r_[simple["AIA_best"],np.zeros(4),m_mu]
    bounds=[(-5,5)]+[(-0.08,0.08)]*4+[(-0.1,0.1)]*4
    def predicted(x):
        A=x[0]; dz=x[1:5]; m=x[5:9]
        coeff=np.array([1.0,A,A*A])
        th=coeff@B
        for ib in range(4):
            th += dz[ib]*(coeff@deriv[ib])
        return apply_mcal(des,th,m)
    def objective(x):
        dz=x[1:5]; m=x[5:9]
        prior=float(np.sum((dz/dz_sig)**2)+np.sum(((m-m_mu)/m_sig)**2))
        return data_chi2(des,predicted(x))+prior
    res=minimize(objective,x0,method="L-BFGS-B",bounds=bounds,
                 options=dict(maxiter=450,ftol=1e-10,gtol=1e-6))
    x=np.asarray(res.x)
    # Exact validation of the shifted source kernels at the emulator optimum.
    Tex=model_theory(model,des,variant,shared_boost_model=shared_boost_model,dz_shift=x[1:5])
    th_exact=Tex[0]+x[0]*Tex[1]+x[0]*x[0]*Tex[2]
    th_exact=apply_mcal(des,th_exact,x[5:9])
    prior=float(np.sum((x[1:5]/dz_sig)**2)+np.sum(((x[5:9]-m_mu)/m_sig)**2))
    exact_data=data_chi2(des,th_exact)
    return dict(
      success=bool(res.success),
      chi2_total_emulator=float(res.fun),
      chi2_data_exact=float(exact_data),
      chi2_prior=float(prior),
      chi2_total_exact=float(exact_data+prior),
      emulator_minus_exact=float(res.fun-(exact_data+prior)),
      AIA_best=float(x[0]),
      dz_best=[float(v) for v in x[1:5]],
      m_best=[float(v) for v in x[5:9]],
      eta_IA_fixed=0.0,
      status="A_IA + 4 source-photo-z + 4 shear-calibration nuisances profiled; first-order dz emulator exact-validated"
    )

def gr_identity(model):
    # Compare local021 native Q_W with the GR Poisson prediction A(z)^2 P_delta.
    h=model["h"]; Om=model["Om"]; H0=100*h/C_KMS
    zs=[0.0,0.5,1.0,2.0]
    ratios=[]
    for z in zs:
        ks=np.geomspace(0.005*h,0.2*h,50) # native 1/Mpc
        pts=np.column_stack([np.full_like(ks,z),np.log(ks)])
        Q=np.exp(model["Iq"](pts)); P=np.exp(model["Ip"](pts))
        pred=(1.5*Om*H0*H0*(1+z))**2*P
        r=Q/pred
        ratios.extend(r[np.isfinite(r)])
    ratios=np.asarray(ratios)
    return dict(median=float(np.median(ratios)),max_abs_minus1=float(np.max(np.abs(ratios-1))),
                p05=float(np.percentile(ratios,5)),p95=float(np.percentile(ratios,95)))

def load_model(name,directory,h,Om):
    p00=Path(directory)/f"{name}_00_background.dat"
    pstd=Path(directory)/f"{name}_background.dat"
    bg=parse_background(p00 if p00.exists() else pstd)
    zq,kq,Q=load_spectral_cube(directory,f"{name}_",h,"weyl")
    # isolate the model's linear matter files; exclude Weyl and nonlinear prefix
    files=sorted(Path(directory).glob(f"{name}_z*_pk.dat"))
    if not files:
        # fallback for one-z naming, though workflow uses multi-z
        files=sorted(Path(directory).glob(f"{name}_*pk.dat"))
    # custom load for exact prefix
    rows=[]
    for p in files:
        if "weyl" in p.name or "nl_" in p.name: continue
        z=z_from_header(p); a=np.loadtxt(p); rows.append((z,a[:,0]*h,a[:,1]/h**3))
    rows.sort(key=lambda x:x[0])
    zp=np.array([r[0] for r in rows]); kp=rows[0][1]; P=np.array([r[2] for r in rows])
    # nonlinear files live under name_nl_z*
    rowsn=[]
    for p in sorted(Path(directory).glob(f"{name}_nl_*pk_nl.dat")):
        if "weyl" in p.name: continue
        z=z_from_header(p); a=np.loadtxt(p); rowsn.append((z,a[:,0]*h,a[:,1]/h**3))
    rowsn.sort(key=lambda x:x[0])
    if rowsn:
        zpn=np.array([r[0] for r in rowsn]); kpn=rowsn[0][1]; Pn=np.array([r[2] for r in rowsn])
    else:
        zpn,kpn,Pn=zp,kp,P
    return dict(name=name,h=h,Om=Om,bg=bg,zq=zq,kq=kq,Q=Q,zp=zp,kp=kp,P=P,
                Iq=cube_interp(zq,kq,Q),Ip=cube_interp(zp,kp,P),Ipnl=cube_interp(zpn,kpn,Pn))

def main():
    root=Path("output/weyl_tomography")
    root.mkdir(parents=True,exist_ok=True)
    data=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
    cuts=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
    if not data.exists(): raise RuntimeError(f"missing DES Y3 data {data}")
    if not cuts.exists(): raise RuntimeError(f"missing DES Y3 scale cuts {cuts}")

    h_l=69.71482083084993/100.0
    Om_l=(0.022083219194622913+0.12299536722293603)/h_l**2
    h_c=68.56859/100.0
    Om_c=(0.02240637+0.11824151)/h_c**2
    late=load_model("late300","output/weyl_spectra",h_l,Om_l)
    lcdm=load_model("local021","output/weyl_spectra",h_c,Om_c)

    gri=gr_identity(lcdm)
    print("WEYL_GR_NORMALIZATION",json.dumps(gri,sort_keys=True))

    allres={"scientific_status":{
      "level":"tomographic DES-Y3 xi+/xi- pilot, not final official likelihood",
      "native_weyl":True,
      "scale_cuts":"official DES-Y3 cuts plus a ladder of increasingly conservative angular cuts",
      "IA":"NLA A_IA profiled, eta_IA fixed 0",
      "nonlinear":"reports linear, shared-LCDM HALOFIT-envelope, and model-specific HALOFIT-envelope variants",
      "missing_nuisance":["photo-z shifts","per-bin shear calibration","baryonic feedback","eta_IA"],
    },"gr_normalization":gri,"results":{}}

    cut_ladder=[
      ("official_cuts",None),
      ("large_20_100",(20.0,100.0)),
      ("large_40_150",(40.0,150.0)),
      ("large_80_200",(80.0,200.0)),
      ("large_120_250",(120.0,250.0)),
    ]
    for label,extra in cut_ladder:
        des=load_des(data,cuts,extra_min=extra)
        key={}
        for variant in ["linear","shared_lcdm_halofit","halofit_envelope"]:
            tv={}
            for model in [late,lcdm]:
                T= model_theory(model,des,variant,shared_boost_model=lcdm)
                prof=chi_profile(des,T[:3])
                prof["projection"]=T[3]
                tv[model["name"]]=prof
            tv["delta_chi2_late300_minus_local021"]=tv["late300"]["chi2_profile"]-tv["local021"]["chi2_profile"]
            key[variant]=tv
        key["ndata"]=len(des["data"])
        key["cov_condition"]=float(np.linalg.cond(des["cov"]))
        allres["results"][label]=key


    # More complete DES-Y3 source-systematics profile on the official shear cuts.
    des_off=load_des(data,cuts,extra_min=None)
    nuisance={}
    for variant in ["linear","shared_lcdm_halofit","halofit_envelope"]:
        vv={}
        for model in [late,lcdm]:
            vv[model["name"]]=chi_profile_des_nuisance(des_off,model,variant,lcdm)
        vv["delta_chi2_total_exact_late300_minus_local021"]=(
            vv["late300"]["chi2_total_exact"]-vv["local021"]["chi2_total_exact"])
        vv["delta_chi2_data_exact_late300_minus_local021"]=(
            vv["late300"]["chi2_data_exact"]-vv["local021"]["chi2_data_exact"])
        nuisance[variant]=vv
    allres["official_cuts_source_nuisance_profile"]=nuisance

    (root/"late300_desy3_weyl_tomography.json").write_text(json.dumps(allres,indent=2,sort_keys=True)+"\n")
    rows=[]
    for cut,v in allres["results"].items():
      for variant,r in v.items():
        if variant in ("ndata","cov_condition"): continue
        for model in ("late300","local021"):
          q=r[model]
          rows.append([cut,variant,model,v["ndata"],q["chi2_A0"],q["chi2_profile"],q["AIA_best"],
                       r["delta_chi2_late300_minus_local021"]])
    import csv
    with (root/"late300_desy3_weyl_tomography.csv").open("w",newline="") as f:
      w=csv.writer(f); w.writerow(["cuts","variant","model","ndata","chi2_AIA0","chi2_profile","AIA_best","delta_chi2_late300_minus_local021"])
      w.writerows(rows)
    print("LATE300_DESY3_WEYL_TOMOGRAPHY",json.dumps(allres,sort_keys=True))

if __name__=="__main__":
    main()
