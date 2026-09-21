#!/usr/bin/env python3
"""
Raw DESI DR1 full-shape promotion for the improved r032 SDMC point.

Models:
  LCDM          : Part-VIII reference control.
  ACOUSTIC_ONLY : improved ordinary cosmology, exact No-Slip (A_lens=0).
  BRAID_BEST    : same ordinary cosmology plus the best stable localized
                  braiding correction from the combined Planck-lite screen.

The official DESI DR1 P0/P2/P4 likelihood blocks are evaluated with the same
REPT/window/nuisance profiling machinery used by the validated Part-VIII gate.
"""
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize

import lsstypes as types
from cosmoprimo import PowerSpectrumInterpolator1D, PowerSpectrumBAOFilter, Cosmology
from cosmoprimo.fiducial import DESI
from velocileptors.EPT.ept_fullresum_varyDz_nu_fftw import REPT

import sys
sys.path.insert(0,"desi-kp-cosmological-likelihoods/dr1/cobaya")
from desi_fs_bao_all import list_zrange, dataset_fn, get_tracer_label, get_physical_stochastic_settings

OUT=Path("output/local002_desi")
OUT.mkdir(parents=True,exist_ok=True)
DATA=Path("desi_likelihood")
C=299792.458

MODELS={
    "LCDM":dict(
        kind="lcdm",
        H0=67.36, ob=0.02237, oc=0.1200,
        As=2.10e-9, ns=0.9649, tau=0.0544,
    ),
    "LOCAL002":dict(
        kind="sdmc",
        H0=70.5653567390982,
        ob=0.022011983189284802,
        oc=0.12404331885203719,
        As=2.132477632355162e-9,
        ns=0.9625227132590487,
        tau=0.055202901571989066,
        AF=0.06017362505197525,
        zc=2.827464461401105,
        width=0.5514493708219379,
    ),
}

# Read the exact effective redshifts and official likelihood blocks.
records=[]
blocks=[]
for tracer,iz,zrange in list_zrange:
    if "lya" in tracer.lower():
        continue
    fn=dataset_fn(str(DATA),tracer,zrange,observable_name="spectrum-poles-rotated")
    d=types.read(fn)
    sp=d.observable.get("spectrum")
    theory=d.window.theory.get("spectrum")
    kin=theory.get(0).coords("k")
    win=d.window.at.observable.get("spectrum").at.theory.get("spectrum").value()
    ns=f"{get_tracer_label(tracer)}_z{iz}"
    zeff=float(sp.attrs["zeff"])
    records.append(dict(tracer=tracer,iz=iz,zrange=list(zrange),namespace=ns,zeff=zeff,file=str(fn)))
    blocks.append(dict(
        tracer=tracer,namespace=ns,zeff=zeff,
        data=sp.value(),
        precision=np.linalg.inv(d.covariance.value()),
        window=win,
        kin=np.asarray(kin,float),
        shotnoise=float(np.mean(sp.get(0).values("shotnoise"))),
    ))
(OUT/"desi_zeff.json").write_text(json.dumps(records,indent=2))
zvals=[0.0]+sorted({float(r["zeff"]) for r in records})
z_to_i={round(z,8):i+1 for i,z in enumerate(zvals)}
print("LOCAL002_DESI_Z",zvals,flush=True)

def ini_lcdm(root,m):
    return textwrap.dedent(f"""\
    H0 = {m['H0']:.15g}
    omega_b = {m['ob']:.15g}
    omega_cdm = {m['oc']:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {m['As']:.17e}
    n_s = {m['ns']:.15g}
    tau_reio = {m['tau']:.15g}
    gauge = synchronous
    modes = s
    output = mPk,dTk,vTk
    extra_metric_transfer_functions = yes
    matter_source_in_current_gauge = no
    P_k_max_h/Mpc = 2.0
    z_pk = {','.join(f'{z:.10g}' for z in zvals)}
    write background = yes
    write thermodynamics = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    output_verbose = 0
    """)

def ini_sdmc(root,m):
    H0=m["H0"]; ob=m["ob"]; oc=m["oc"]
    h=H0/100.
    OR=4.17998772e-5
    OX=1.-(ob+oc+OR)/(h*h)
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {ob:.15g}
    omega_cdm = {oc:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {m['As']:.17e}
    n_s = {m['ns']:.15g}
    tau_reio = {m['tau']:.15g}

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {m['AF']:.17g},{m['zc']:.17g},{m['width']:.17g},0.34231919445927034,1.0,0.045
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},17.925,17.775,0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    a_ini_over_a_today_default = 1.e-8
    a_ini_test_qs_smg = 1.e-8
    pert_ic_ini_z_ref_smg = 1.e7
    a_min_stability_test_smg = 1.e-8

    gauge = synchronous
    modes = s
    output = mPk,dTk,vTk
    extra_metric_transfer_functions = yes
    matter_source_in_current_gauge = no
    P_k_max_h/Mpc = 2.0
    z_pk = {','.join(f'{z:.10g}' for z in zvals)}
    output_background_smg = 3
    write background = yes
    write thermodynamics = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    output_verbose = 0
    """)

# Run CLASS for each model.
for label,m in MODELS.items():
    prefix=label.lower()
    root=str(OUT/(prefix+"_"))
    ip=OUT/(prefix+".ini")
    ip.write_text(ini_lcdm(root,m) if m["kind"]=="lcdm" else ini_sdmc(root,m))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=600)
    print("LOCAL002_DESI_CLASS",label,cp.returncode,flush=True)
    if cp.returncode!=0:
        print(cp.stdout[-4000:],flush=True)
        raise SystemExit(f"CLASS failed for {label}")

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def rd_from(prefix,bg):
    th=tab(OUT/f"{prefix}_00_thermodynamics.dat").sort_values("z")
    zz=th.z.to_numpy()
    kb=th.kappa_b.to_numpy()
    cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i])
            cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross:
        raise RuntimeError(f"no drag crossing for {prefix}")
    zd=float(cross[0])
    return float(np.interp(zd,bg.z,bg["comov.snd.hrz."])),zd

def model_data(label,m):
    prefix=label.lower()
    H0=m["H0"]; ob=m["ob"]; oc=m["oc"]; ns=m["ns"]
    h=H0/100.
    fb=ob/(ob+oc)
    fc=1.-fb
    bg=tab(OUT/f"{prefix}_00_background.dat").sort_values("z")
    rd,zd=rd_from(prefix,bg)
    series={}
    for z in zvals:
        ii=z_to_i[round(z,8)]
        pk=np.loadtxt(OUT/f"{prefix}_00_z{ii}_pk.dat")
        tk=tab(OUT/f"{prefix}_00_z{ii}_tk.dat")
        kh=tk["k (h/Mpc)"].to_numpy()
        dm=tk["d_m"].to_numpy()
        tb=tk["t_b"].to_numpy()
        hp=tk["h_prime"].to_numpy()
        ep=tk["eta_prime"].to_numpy()
        shift=.5*(hp+6.*ep)
        a=1./(1.+z)
        H=float(np.interp(z,bg.z,bg["H [1/Mpc]"]))
        Hc=a*H
        theta_b_gi=tb+shift
        theta_c_gi=shift
        vb=-theta_b_gi/Hc
        vc=-theta_c_gi/Hc
        vcb=fb*vb+fc*vc
        ratio_v=np.divide(vcb,dm,out=np.zeros_like(vcb),where=np.abs(dm)>1e-300)
        kpk=pk[:,0]
        Pdd=pk[:,1]
        rv=np.interp(kpk,kh,ratio_v)
        Ptt=Pdd*rv**2
        series[z]=(kpk,Pdd,Ptt)
    rec=dict(label=label,prefix=prefix,H0=H0,h=h,ob=ob,oc=oc,fb=fb,fc=fc,ns=ns,bg=bg,rd=rd,zd=zd,series=series)
    if m["kind"]=="sdmc":
        rec["min_D"]=float(bg["kin (D)"].min())
        rec["min_cs2"]=float(bg["c_s^2"].min())
        rec["max_cs2"]=float(bg["c_s^2"].max())
        rec["max_abs_slip_driver"]=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"])))
        if not (rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1.):
            raise SystemExit(f"stability failed for {label}: {rec}")
        print("LOCAL002_DESI_STABILITY",json.dumps({k:v for k,v in rec.items() if k not in ("bg","series")},sort_keys=True),flush=True)
    return rec

models={label:model_data(label,m) for label,m in MODELS.items()}

kobs=blocks[0]["kin"]
for b in blocks:
    if not np.allclose(b["kin"],kobs):
        raise RuntimeError("DESI theory k grids differ")

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
    kin=np.geomspace(min(5e-4,kobs[0]/2),max(1.0,kobs[-1]*2),500)
    zs=np.array(zvals)
    Pdd=np.stack([np.interp(kin,*m["series"][z][:2]) for z in zs],axis=-1)
    Ptt=np.stack([np.interp(kin,m["series"][z][0],m["series"][z][2]) for z in zs],axis=-1)
    pki=PowerSpectrumInterpolator1D(kin,Pdd)
    cosmo=Cosmology(n_s=m["ns"],Omega_b=m["ob"]/m["h"]**2,
                    Omega_cdm=m["oc"]/m["h"]**2,Omega_ncdm=0.,H0=m["H0"])
    cosmo.rs_drag=m["rd"]*m["h"]
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
    pcb=loginterp(Pdd)
    pnw=loginterp(Pnw)
    ptt=loginterp(Ptt)

    tables={}
    for b in blocks:
        z=b["zeff"]
        iz=zvals.index(z)
        Hkms=float(np.interp(z,m["bg"].z,m["bg"]["H [1/Mpc]"]))*C
        DA=float(np.interp(z,m["bg"].z,m["bg"]["comov. dist."]))/(1.+z)
        qpar=float(fid.efunc(z)/(Hkms/(100.*m["h"])))
        qper=float(DA*m["h"]/fid.angular_diameter_distance(z))
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
            qpar=qpar,qper=qper,Dz=Dz
        )
    return tables

def get_poles(info,pars,settings,sn,return_gradient=False):
    pktable=info["basis"]
    sigma8=info["sigma8"]
    f=info["fsigma8"]/sigma8
    b1p,b2p,bsp=pars
    b3p=0.
    b1L=b1p/sigma8-1.
    b2L=b2p/sigma8**2
    bsL=bsp/sigma8**2
    b3L=b3p/sigma8**3
    b1=1.+b1L
    b2=8./21.*b1L+b2L
    bs=bsL-(2./7.)*b1L
    b3=3*b3L+b1L

    gradient=np.zeros((7,7))
    gradient[0,0]=b1*b1
    gradient[1,0]=gradient[1,1]=f*b1
    gradient[2,1]=f*f
    gradient[2,2]=f*b1
    gradient[3,2]=f*f
    fsat=settings["fsat"]
    sigv=settings["sigv"]
    for jj,pole in enumerate([0,2,4]):
        gradient[jj-3,jj-3]=sn*(fsat if pole>0 else 1.)*sigv**pole
    base=[1,b1,b1*b1,b2,b1*b2,b2*b2,bs,b1*bs,b2*bs,bs*bs,b3,b1*b3]
    nuis=[0.]*7
    mono=np.concatenate([np.asarray(base),gradient.dot(np.asarray(nuis))])
    poles=np.sum(pktable*mono,axis=-1)
    if return_gradient:
        return poles,pktable[...,-7:].dot(gradient)
    return poles

def profile_model(label,m):
    tables=prepare_basis(m)
    rows=[]
    total=0.
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
            opti=minimize(lambda x:-2*block_logp(x),seed,method="L-BFGS-B",
                          bounds=[(0.,3.),(-20.,20.),(-20.,20.)],
                          options={"maxiter":800,"ftol":1e-11,"gtol":3e-7,"maxls":50})
            trials.append(opti)
            print("LOCAL002_DESI_ROBUST_TRIAL",label,b["namespace"],iseed,float(opti.fun),
                  bool(opti.success),[float(v) for v in opti.x],flush=True)
        finite=[o for o in trials if np.isfinite(o.fun)]
        opt=min(finite,key=lambda o:o.fun)
        vals=np.array([float(o.fun) for o in finite])
        chi=float(opt.fun)
        total+=chi
        info=tables[b["namespace"]]
        row=dict(
            model=label,tracer=b["tracer"],namespace=b["namespace"],
            zeff=b["zeff"],chi2_profile=chi,success=bool(opt.success),
            start_spread_chi2=float(vals.max()-vals.min()) if len(vals)>1 else 0.,
            n_success=sum(bool(o.success) for o in trials),
            b1p=float(opt.x[0]),b2p=float(opt.x[1]),bsp=float(opt.x[2]),
            sigma8=info["sigma8"],fsigma8=info["fsigma8"],
            qpar=info["qpar"],qper=info["qper"],Dz=info["Dz"]
        )
        rows.append(row)
        print("LOCAL002_DESI_BLOCK",json.dumps(row,sort_keys=True),flush=True)
    print("LOCAL002_DESI_TOTAL",label,total,flush=True)
    return rows,total

allrows=[]
totals={}
for label,m in models.items():
    rows,total=profile_model(label,m)
    allrows.extend(rows)
    totals[label]=total

summary={
    "LCDM":totals["LCDM"],
    "SOBOL049":totals["SOBOL049"],
    "SOBOL009":totals["SOBOL009"],
    "delta_sobol049_minus_lcdm":totals["SOBOL049"]-totals["LCDM"],
    "delta_sobol009_minus_lcdm":totals["SOBOL009"]-totals["LCDM"],
}
pd.DataFrame(allrows).to_csv(OUT/"local002_desi_profile.csv",index=False)
(OUT/"local002_desi_summary.json").write_text(json.dumps(summary,indent=2))
print("LOCAL002_DESI_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
