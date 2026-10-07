#!/usr/bin/env python3
"""
KiDS-Legacy native-Weyl nuisance-central COSEBI inputs for late300 vs local021.

Uses:
  * six-bin KiDS-Legacy NZ_SOURCE distributions;
  * published correlated photo-z central shifts;
  * published mass-dependent IA central point (A=5.74, beta=0.44 after
    the KiDS covariance transform);
  * native hi_class Weyl power for GG;
  * native Weyl-matter cross approximation sqrt(Q_W P_m) for GI;
  * linear matter power for II.

This is still a linear/native-Weyl pilot.  The full KiDS COSEBI statistic mixes
nonlinear scales, and there is no calibrated nonlinear SDMC Weyl prescription
in the repository yet.
"""
from pathlib import Path
import importlib.util, json, numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
from astropy.io import fits

S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)

DATA=Path("cosmosis-standard-library/likelihood/KiDS-Legacy/KiDS_Legacy_cosebis.fits")
ROOT=Path("cosmosis-standard-library/likelihood/KiDS-Legacy")
OUT=Path("output/kids_native_weyl_nuisance"); OUT.mkdir(parents=True,exist_ok=True)
ELL=np.geomspace(2.0,2.0e4,800)

hL=0.6971482083084993
OmL=(0.022083219194622913+0.12299536722293603)/hL**2
hC=0.6856859
OmC=(0.02240637+0.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)

with fits.open(DATA) as f:
    nz=f["NZ_SOURCE"].data
    zsrc=np.asarray(nz["Z_MID"],float)
    bins=np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])

# KiDS-Legacy latent Gaussian central values from examples/KiDS-Legacy_values.ini.
unc_dz=np.array([2.5374,-2.44484,-1.37982,-0.494275,1.40234,6.26323])
unc_mass=np.array([
    19.77069603580704893,30.74866003109325874,
    44.30836735636641777,-32.06136674281973,
    143.9212386242852517,-137.8196728635030013,
    -202.8513316870826202,419.8363708097585913
])
Cdz=np.loadtxt(ROOT/"Nz_covariance.txt")
Cmass=np.loadtxt(ROOT/"massdep_cov.txt")
dz=np.linalg.cholesky(Cdz)@unc_dz
massphys=np.linalg.cholesky(Cmass)@unc_mass
A=float(massphys[0]); beta=float(massphys[1]); logM=np.asarray(massphys[2:])
logMpiv=13.5
fr=np.array([0.158,0.198,0.206,0.258,0.207,0.026])
ia_bin_amp=A*fr*10.0**(beta*(logM-logMpiv))

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.0))[::-1]

def make(model):
    zmax=min(model["zq"].max(),model["zp"].max(),float(zsrc.max()))
    z=np.linspace(max(0.01,float(zsrc[zsrc>0].min())),zmax,420)
    H=np.interp(z,model["bg"]["z"],model["bg"]["H"])
    chi=np.interp(z,model["bg"]["z"],model["bg"]["chi"])
    nzs=[]
    for i,arr in enumerate(bins):
        f=interp1d(zsrc,arr,bounds_error=False,fill_value=0.0)
        # KiDS source-photoz-bias convention: positive bias shifts sources up in z.
        y=np.maximum(f(z-dz[i]),0.0)
        norm=np.trapezoid(y,z)
        if norm<=0: raise RuntimeError(f"bad n(z) norm bin {i+1}")
        nzs.append(y/norm)
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
    valid=np.isfinite(Q)&np.isfinite(P)
    Q=np.where(valid,Q,0.0); P=np.where(valid,P,0.0)
    X=np.sqrt(np.maximum(Q*P,0.0))

    kref=0.05
    pref=np.exp(model["Ip"](np.column_stack([z,np.full_like(z,np.log(kref))])))
    p0=float(np.exp(model["Ip"]([[0.0,np.log(kref)]])[0]))
    D=np.maximum(np.sqrt(pref/p0),1e-4)
    Fbase=-b.C1RHO*model["Om"]/D
    Fi=ia_bin_amp[:,None]*Fbase[None,:]

    comp={"gg":{"ell":ELL},"gi":{"ell":ELL},"ii":{"ell":ELL},"total":{"ell":ELL}}
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
            key1=f"bin_{j+1}_{i+1}"; key2=f"bin_{i+1}_{j+1}"
            for key,val in [("gg",gg),("gi",gi),("ii",ii),("total",gg+gi+ii)]:
                comp[key][key1]=val; comp[key][key2]=val
    meta=dict(
        model=model["name"],zmin=float(z.min()),zmax=float(z.max()),
        ell_min=float(ELL.min()),ell_max=float(ELL.max()),
        valid_fraction=float(valid.mean()),
        dz=dz.tolist(),A_IA=A,beta=beta,log10_M_mean=logM.tolist(),
        f_r=fr.tolist(),ia_bin_amplitude=ia_bin_amp.tolist(),
        note="Published KiDS-Legacy nuisance central point; linear native-Weyl GG+GI+II; no nonlinear SDMC Weyl calibration."
    )
    return comp,meta

summary={
  "nuisance_central":{
    "dz":dz.tolist(),"mass_physical":massphys.tolist(),
    "ia_bin_amplitude":ia_bin_amp.tolist()
  },
  "models":{}
}
for model in (late,lcdm):
    comp,meta=make(model)
    for variant,arr in comp.items():
        np.savez(OUT/f"{model['name']}_kids_{variant}_cl.npz",**arr)
    summary["models"][model["name"]]=meta
    print("KIDS_NATIVE_WEYL_NUISANCE_CL",json.dumps(meta,sort_keys=True),flush=True)

(OUT/"kids_native_weyl_nuisance_cl_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
