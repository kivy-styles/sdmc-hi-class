#!/usr/bin/env python3
"""
Official DESI DR1 full-shape replay for the exact-covariant geometry-compromise
SDMC leader, using the preserved exact-covariant transfer artifact from run
35574458587.

The covariant artifact was generated at z = 0.295, 0.510, 0.706, 0.919,
1.317, 1.491.  The official DESI effective redshifts differ from those nodes
by <= 3.7e-4.  For each official z_eff we preserve the exact covariant
P(k)-shape from the nearest node and rescale its density and velocity spectra
by sigma8(z)^2 and [f sigma8(z)]^2 interpolated across the exact-covariant
nodes.  Geometry (H, D_A) is evaluated from the exact covariant background at
the exact official z_eff.

This is a deliberately conservative artifact replay: the redshift correction
is many orders smaller than the model-vs-LCDM signal seen in the parameterized
DESI gate, while retaining the exact covariant growth solution.
"""
from pathlib import Path
import sys,re,json,math
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize

sys.path.insert(0,"desi-kp-cosmological-likelihoods/dr1/cobaya")
import lsstypes as types
from cosmoprimo import PowerSpectrumInterpolator1D, PowerSpectrumBAOFilter, Cosmology
from cosmoprimo.fiducial import DESI
from velocileptors.EPT.ept_fullresum_varyDz_nu_fftw import REPT
from desi_fs_bao_all import list_zrange, dataset_fn, get_tracer_label, get_physical_stochastic_settings

ROOT=Path(".")
OUT=Path("output/covariant_geometry_desi")
OUT.mkdir(parents=True,exist_ok=True)
COV=Path("input_cov")
LCDM=Path("output")
DATA=Path("desi_likelihood")
C=299792.458

COVPAR=dict(
    H0=70.5659273196888,
    ob=0.022329086486197414,
    oc=0.12479662145108897,
    ns=0.963943012708798,
)
LCDMPAR=dict(H0=67.36,ob=0.02237,oc=0.1200,ns=0.9649)

def tab(path):
    path=Path(path)
    lines=path.read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def sigmaR_arrays(k,P,R=8.0):
    k=np.asarray(k,float); P=np.asarray(P,float)
    x=k*R
    W=np.ones_like(x); q=x!=0
    W[q]=3*(np.sin(x[q])-x[q]*np.cos(x[q]))/x[q]**3
    y=k**3*P*W**2/(2*np.pi**2)
    return float(np.sqrt(np.trapezoid(y,x=np.log(k))))

# Official likelihood blocks and exact z_eff.
records=[]; blocks=[]
for tracer,iz,zrange in list_zrange:
    if "lya" in tracer.lower(): continue
    fn=dataset_fn(str(DATA),tracer,zrange,observable_name="spectrum-poles-rotated")
    d=types.read(fn)
    sp=d.observable.get("spectrum")
    theory=d.window.theory.get("spectrum")
    kin=np.asarray(theory.get(0).coords("k"),float)
    win=d.window.at.observable.get("spectrum").at.theory.get("spectrum").value()
    ns=f"{get_tracer_label(tracer)}_z{iz}"
    zeff=float(sp.attrs["zeff"])
    records.append(dict(tracer=tracer,iz=iz,zrange=list(zrange),namespace=ns,zeff=zeff,file=str(fn)))
    blocks.append(dict(
        tracer=tracer,namespace=ns,zeff=zeff,
        data=sp.value(),precision=np.linalg.inv(d.covariance.value()),
        window=win,kin=kin,shotnoise=float(np.mean(sp.get(0).values("shotnoise")))
    ))
(OUT/"desi_zeff.json").write_text(json.dumps(records,indent=2))
z_eff=np.array(sorted({float(r["zeff"]) for r in records}),float)
print("COV_GEOM_DESI_ZEFF",z_eff.tolist(),flush=True)

# Exact covariant background and drag ruler from preserved artifact.
cov_bg=tab(COV/"covariant_transfer_00_background.dat").sort_values("z")
ac=pd.read_csv(COV/"covariant_observable_acoustic.csv").iloc[0]
cov_rd=float(ac["rd_Mpc"])

# Reconstruct exact-covariant density and velocity spectra on official z_eff.
z_nodes=np.array([0.295,0.510,0.706,0.919,1.317,1.491],float)
if np.max(np.abs(z_eff-z_nodes))>5e-4:
    raise RuntimeError(f"unexpected DESI z_eff shift: {z_eff-z_nodes}")

h=COVPAR["H0"]/100.
fb=COVPAR["ob"]/(COVPAR["ob"]+COVPAR["oc"])
fc=1.-fb
node_series=[]
node_s8=[]
node_fs8=[]
for i,z in enumerate(z_nodes,1):
    pk=np.loadtxt(COV/f"covariant_transfer_00_z{i}_pk.dat")
    tk=tab(COV/f"covariant_transfer_00_z{i}_tk.dat")
    kh=tk["k (h/Mpc)"].to_numpy()
    dm=tk["d_m"].to_numpy()
    tb=tk["t_b"].to_numpy()
    hp=tk["h_prime"].to_numpy()
    ep=tk["eta_prime"].to_numpy()
    shift=.5*(hp+6.*ep)
    a=1./(1.+z)
    H=float(np.interp(z,cov_bg.z,cov_bg["H [1/Mpc]"]))
    Hc=a*H
    vb=-(tb+shift)/Hc
    vc=-shift/Hc
    vcb=fb*vb+fc*vc
    ratio_v=np.divide(vcb,dm,out=np.zeros_like(vcb),where=np.abs(dm)>1e-300)
    kpk=pk[:,0]
    Pdd=pk[:,1]
    rv=np.interp(kpk,kh,ratio_v)
    Ptt=Pdd*rv**2
    s8=sigmaR_arrays(kpk,Pdd)
    fs8=sigmaR_arrays(kpk,Ptt)
    node_series.append((kpk,Pdd,Ptt))
    node_s8.append(s8); node_fs8.append(fs8)
    print("COV_GEOM_NODE",z,s8,fs8,flush=True)

node_s8=np.array(node_s8); node_fs8=np.array(node_fs8)
s8_eff=np.interp(z_eff,z_nodes,node_s8)
fs8_eff=np.interp(z_eff,z_nodes,node_fs8)

cov_series={}
corr_rows=[]
for j,z in enumerate(z_eff):
    i=int(np.argmin(np.abs(z_nodes-z)))
    zn=z_nodes[i]
    k,Pdd,Ptt=node_series[i]
    sd=node_s8[i]; st=node_fs8[i]
    sd_t=s8_eff[j]; st_t=fs8_eff[j]
    Pdd_t=Pdd*(sd_t/sd)**2
    Ptt_t=Ptt*(st_t/st)**2
    cov_series[float(z)]=(k,Pdd_t,Ptt_t)
    corr_rows.append(dict(
        zeff=float(z),source_z=float(zn),delta_z=float(z-zn),
        sigma8_source=float(sd),sigma8_corrected=float(sd_t),
        fsigma8_source=float(st),fsigma8_corrected=float(st_t),
        density_power_scale=float((sd_t/sd)**2),
        velocity_power_scale=float((st_t/st)**2)
    ))
pd.DataFrame(corr_rows).to_csv(OUT/"covariant_redshift_correction.csv",index=False)
print("COV_GEOM_REDSHIFT_CORRECTION",pd.DataFrame(corr_rows).to_dict("records"),flush=True)

# z=0 density spectrum for REPT base shape from the same covariant artifact.
cov_z0=np.loadtxt(COV/"covariant_observable_00_z1_pk.dat")
cov_k0,cov_P0=cov_z0[:,0],cov_z0[:,1]
f0=float(np.interp(0.0,cov_bg.z,cov_bg["gr.fac. f"]))
cov_series[0.0]=(cov_k0,cov_P0,cov_P0*f0*f0)

# LCDM exact-z outputs generated in the workflow.
lcdm_bg=tab(LCDM/"lcdm_00_background.dat").sort_values("z")
def rd_from(prefix,bg):
    th=tab(LCDM/f"{prefix}_00_thermodynamics.dat").sort_values("z")
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy()
    cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i])
            cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross: raise RuntimeError("no LCDM drag crossing")
    zd=float(cross[0])
    return float(np.interp(zd,bg.z,bg["comov.snd.hrz."])),zd
lcdm_rd,lcdm_zd=rd_from("lcdm",lcdm_bg)

lcdm_zvals=[0.0]+list(z_eff)
lcdm_series={}
hL=LCDMPAR["H0"]/100.
fbL=LCDMPAR["ob"]/(LCDMPAR["ob"]+LCDMPAR["oc"])
fcL=1.-fbL
for i,z in enumerate(lcdm_zvals,1):
    pk=np.loadtxt(LCDM/f"lcdm_00_z{i}_pk.dat")
    tk=tab(LCDM/f"lcdm_00_z{i}_tk.dat")
    kh=tk["k (h/Mpc)"].to_numpy()
    dm=tk["d_m"].to_numpy()
    tb=tk["t_b"].to_numpy()
    hp=tk["h_prime"].to_numpy(); ep=tk["eta_prime"].to_numpy()
    shift=.5*(hp+6.*ep)
    a=1./(1.+z)
    H=float(np.interp(z,lcdm_bg.z,lcdm_bg["H [1/Mpc]"]))
    Hc=a*H
    vb=-(tb+shift)/Hc
    vc=-shift/Hc
    vcb=fbL*vb+fcL*vc
    ratio_v=np.divide(vcb,dm,out=np.zeros_like(vcb),where=np.abs(dm)>1e-300)
    kpk=pk[:,0]; Pdd=pk[:,1]
    rv=np.interp(kpk,kh,ratio_v)
    Ptt=Pdd*rv**2
    lcdm_series[float(z)]=(kpk,Pdd,Ptt)

models={
 "COVARIANT_GEOM":dict(**COVPAR,bg=cov_bg,rd=cov_rd,series=cov_series),
 "LCDM":dict(**LCDMPAR,bg=lcdm_bg,rd=lcdm_rd,series=lcdm_series)
}

kobs=blocks[0]["kin"]
for b in blocks:
    if not np.allclose(b["kin"],kobs): raise RuntimeError("DESI theory k grids differ")

starts={
    "BGS_z0":(1.11348,0.660148,-0.223088),
    "LRG_z0":(1.1351,-0.139101,-0.910397),
    "LRG_z1":(1.22446,-0.671727,-0.22624),
    "LRG_z2":(1.04564,-0.493139,-0.24021),
    "ELG_z1":(0.127936,-0.578801,-1.45356),
    "QSO_z0":(0.780983,0.353843,0.135725),
}

fid=DESI(engine="camb")
all_marg=["alpha0p","alpha2p","alpha4p","alpha6p","sn0p","sn2p","sn4p"]
scales=np.array([12.5]*4+[2.]+[5.]*2)
marg=["alpha0p","alpha2p","sn0p","sn2p"]
gi=np.array([all_marg.index(x) for x in marg])
prior_hess=-np.diag(scales[gi]**-2)

def prepare_basis(m):
    zvals=np.array([0.0]+list(z_eff),float)
    kin=np.geomspace(min(5e-4,kobs[0]/2),max(1.0,kobs[-1]*2),500)
    Pdd=np.stack([np.interp(kin,*m["series"][float(z)][:2]) for z in zvals],axis=-1)
    Ptt=np.stack([np.interp(kin,m["series"][float(z)][0],m["series"][float(z)][2]) for z in zvals],axis=-1)
    pki=PowerSpectrumInterpolator1D(kin,Pdd)
    h=m["H0"]/100.
    cosmo=Cosmology(n_s=m["ns"],Omega_b=m["ob"]/h**2,
                    Omega_cdm=m["oc"]/h**2,Omega_ncdm=0.,H0=m["H0"])
    cosmo.rs_drag=m["rd"]*h
    filt=PowerSpectrumBAOFilter(pki,engine="peakaverage",cosmo=cosmo,cosmo_fid=fid)
    filt(pki,cosmo=cosmo)
    Pnw=filt.smooth_pk_interpolator()(kin)
    sig8=np.asarray(pki.sigma8())
    ptti=PowerSpectrumInterpolator1D(kin,Ptt)
    fs8=np.asarray(ptti.sigma8())
    pt=REPT(kin,Pdd[:,0],pnw=Pnw[:,0],kmin=kobs[0],kmax=kobs[-1],
            nk=200,rbao=110,sbao=None,beyond_gauss=True,one_loop=True,
            shear=True,cutoff=20,jn=5,N=4000,threads=2,
            extrap_min=-4,extrap_max=3,import_wisdom=False)
    extk=np.append(pt.kv,1.)
    def loginterp(arr):
        return 10**interp1d(np.log10(kin),np.log10(np.maximum(arr,1e-300)),
                            kind="cubic",fill_value="extrapolate",axis=0,
                            assume_sorted=True)(np.log10(extk))
    pcb=loginterp(Pdd); pnw=loginterp(Pnw); ptt=loginterp(Ptt)
    tables={}
    for b in blocks:
        z=b["zeff"]; iz=int(np.argmin(np.abs(zvals-z)))
        Hkms=float(np.interp(z,m["bg"].z,m["bg"]["H [1/Mpc]"]))*C
        DA=float(np.interp(z,m["bg"].z,m["bg"]["comov. dist."]))/(1.+z)
        qpar=float(fid.efunc(z)/(Hkms/(100.*h)))
        qper=float(DA*h/fid.angular_diameter_distance(z))
        Dz=float(np.sqrt(pcb[-1,iz]/pcb[-1,0]))
        fk=np.sqrt(np.maximum(ptt[:-1,iz]/pcb[:-1,iz],0.))
        pks=pt.compute_redshift_space_power_multipoles_tables(
            fk,apar=qpar,aperp=qper,ngauss=4,
            pcb=pcb[:-1,iz],pcb_nw=pnw[:-1,iz],Dz=Dz)[1:]
        basis=np.stack([
            interp1d(pt.kv,pks[j],kind="cubic",fill_value="extrapolate",
                     axis=0,assume_sorted=True)(kobs)
            for j in range(3)
        ],axis=0)
        tables[b["namespace"]]=dict(
            basis=basis,sigma8=float(sig8[iz]),fsigma8=float(fs8[iz]),
            qpar=qpar,qper=qper,Dz=Dz)
    return tables

def get_poles(info,pars,settings,sn,return_gradient=False):
    pktable=info["basis"]; sigma8=info["sigma8"]; f=info["fsigma8"]/sigma8
    b1p,b2p,bsp=pars; b3p=0.
    b1L=b1p/sigma8-1.; b2L=b2p/sigma8**2
    bsL=bsp/sigma8**2; b3L=b3p/sigma8**3
    b1=1.+b1L; b2=8./21.*b1L+b2L; bs=bsL-(2./7.)*b1L; b3=3*b3L+b1L
    gradient=np.zeros((7,7))
    gradient[0,0]=b1*b1
    gradient[1,0]=gradient[1,1]=f*b1
    gradient[2,1]=f*f
    gradient[2,2]=f*b1
    gradient[3,2]=f*f
    fsat=settings["fsat"]; sigv=settings["sigv"]
    for jj,pole in enumerate([0,2,4]):
        gradient[jj-3,jj-3]=sn*(fsat if pole>0 else 1.)*sigv**pole
    base=[1,b1,b1*b1,b2,b1*b2,b2*b2,bs,b1*bs,b2*bs,bs*bs,b3,b1*b3]
    nuis=[0.]*7
    mono=np.concatenate([np.asarray(base),gradient.dot(np.asarray(nuis))])
    poles=np.sum(pktable*mono,axis=-1)
    if return_gradient: return poles,pktable[...,-7:].dot(gradient)
    return poles

def profile_model(label,m):
    tables=prepare_basis(m)
    rows=[]; total=0.
    for b in blocks:
        settings=get_physical_stochastic_settings(tracer=b["tracer"].upper()[:3])
        def block_logp(x):
            poles,grad=get_poles(tables[b["namespace"]],x,settings,b["shotnoise"],True)
            theory=b["window"].dot(poles.ravel())
            grad=grad[...,gi].reshape(-1,len(gi))
            grad=b["window"].dot(grad)
            diff=theory-b["data"]
            pgrad=b["precision"].dot(grad)
            postgrad=-pgrad.T.dot(diff)
            likeh=-grad.T.dot(pgrad)
            posth=prior_hess+likeh
            dx=-np.linalg.solve(posth,postgrad)
            logl=-.5*diff.T.dot(b["precision"]).dot(diff)
            logl+=.5*dx.dot(likeh).dot(dx)+postgrad.dot(dx)
            logprior=.5*dx.dot(prior_hess).dot(dx)
            logl+=-.5*np.linalg.slogdet(-posth)[1]
            logbias=-.5*(x[1]/5.)**2-.5*(x[2]/5.)**2
            return float(logl+logprior+logbias)
        seeds=[
            np.array(starts[b["namespace"]],float),
            np.array([1.0,0.0,0.0]),
            np.array([0.55,4.0,-4.0]),
            np.array([1.8,-4.0,4.0]),
            np.array([2.6,8.0,-8.0]),
            np.array([1.4,-10.0,-10.0]),
        ]
        trials=[]
        for iseed,seed in enumerate(seeds):
            opt=minimize(lambda x:-2*block_logp(x),seed,method="L-BFGS-B",
                         bounds=[(0.,3.),(-20.,20.),(-20.,20.)],
                         options={"maxiter":800,"ftol":1e-11,"gtol":3e-7,"maxls":50})
            trials.append(opt)
            print("COV_GEOM_DESI_TRIAL",label,b["namespace"],iseed,float(opt.fun),
                  bool(opt.success),[float(v) for v in opt.x],flush=True)
        finite=[o for o in trials if np.isfinite(o.fun)]
        opt=min(finite,key=lambda o:o.fun)
        vals=np.array([float(o.fun) for o in finite])
        chi=float(opt.fun); total+=chi
        info=tables[b["namespace"]]
        row=dict(
            model=label,tracer=b["tracer"],namespace=b["namespace"],
            zeff=b["zeff"],chi2_profile=chi,success=bool(opt.success),
            start_spread_chi2=float(vals.max()-vals.min()) if len(vals)>1 else 0.,
            n_success=sum(bool(o.success) for o in trials),
            b1p=float(opt.x[0]),b2p=float(opt.x[1]),bsp=float(opt.x[2]),
            sigma8=info["sigma8"],fsigma8=info["fsigma8"],
            qpar=info["qpar"],qper=info["qper"],Dz=info["Dz"])
        rows.append(row)
        print("COV_GEOM_DESI_BLOCK",json.dumps(row,sort_keys=True),flush=True)
    print("COV_GEOM_DESI_TOTAL",label,total,flush=True)
    return rows,total

allrows=[]; totals={}
for label,m in models.items():
    rr,tot=profile_model(label,m)
    allrows.extend(rr); totals[label]=tot

delta=totals["COVARIANT_GEOM"]-totals["LCDM"]
summary=dict(
    covariant_chi2=totals["COVARIANT_GEOM"],
    lcdm_chi2=totals["LCDM"],
    delta_chi2_covariant_minus_lcdm=delta,
    fullplik_delta_covariant_minus_lcdm=5.36115841529,
    combined_planck_plus_desi_delta=5.36115841529+delta,
)
pd.DataFrame(allrows).to_csv(OUT/"covariant_geometry_desi_profile.csv",index=False)
(OUT/"covariant_geometry_desi_summary.json").write_text(json.dumps(summary,indent=2))
print("COV_GEOM_DESI_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
