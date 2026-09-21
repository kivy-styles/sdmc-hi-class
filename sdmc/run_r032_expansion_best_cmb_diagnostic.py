#!/usr/bin/env python3
from pathlib import Path
import math,re,subprocess,json
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

OUT=Path("output/cmb_diagnostic"); OUT.mkdir(parents=True,exist_ok=True)
OR=4.17998772e-5

def run_ini(name,text):
    p=OUT/f"{name}.ini"; p.write_text(text)
    cp=subprocess.run(["./class",str(p)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    if cp.returncode:
        raise RuntimeError(cp.stdout[-3000:])

# exact coordinates of the promoted parameterized expansion-best model
H0=70.5653567390982; ob=0.022011983189284802; oc=0.12404331885203719
ns=0.9625227132590487; tau=0.055202901571989066; lnAs=3.0598696043919773
As=math.exp(lnAs)/1e10; h=H0/100.; OX=1.-(ob+oc+OR)/(h*h)
sdmc=f"""H0={H0}
omega_b={ob}
omega_cdm={oc}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={As:.17e}
n_s={ns}
tau_reio={tau}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg=0.06017362505197525,2.827464461401105,0.5514493708219379,0.34231919445927034,1.0,0.045
expansion_model=sdmc_full
expansion_smg={OX:.17g},17.7,17.1,0.5,0.01105624999,0.25,0.01951933685,1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
output_background_smg=3
modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
write background=yes
write thermodynamics=yes
root={OUT}/sdmc_
format=class
input_verbose=0
background_verbose=0
thermodynamics_verbose=0
perturbations_verbose=0
spectra_verbose=0
lensing_verbose=0
output_verbose=0
"""
lcdm=f"""H0=68.56858744695782
omega_b=0.022406369378007947
omega_cdm=0.11824151052483357
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s=2.05109266435559e-9
n_s=0.9649593164240942
tau_reio=0.044985382026527077
modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
write background=yes
write thermodynamics=yes
root={OUT}/lcdm021_
format=class
input_verbose=0
background_verbose=0
thermodynamics_verbose=0
perturbations_verbose=0
spectra_verbose=0
lensing_verbose=0
output_verbose=0
"""
run_ini("sdmc",sdmc); run_ini("lcdm021",lcdm)

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def acoustic(prefix):
    bg=table(OUT/f"{prefix}_00_background.dat")
    th=table(OUT/f"{prefix}_00_thermodynamics.dat")
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    z=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); zd=None
    for i in range(len(z)-1):
        if (kb[i]-1)*(kb[i+1]-1)<=0 and kb[i]!=kb[i+1]:
            q=(1-kb[i])/(kb[i+1]-kb[i]); zd=float(z[i]+q*(z[i+1]-z[i])); break
    if zd is None: raise RuntimeError("drag crossing missing")
    rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    rd=float(np.interp(zd,bg.z,bg["comov.snd.hrz."]))
    dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return dict(z_star=zstar,z_drag=zd,rs_star=rs,rd=rd,DM_star=dm,ell_A=math.pi*dm/rs)

def cls(prefix,lensed=True):
    suffix="cl_lensed" if lensed else "cl"
    a=np.loadtxt(OUT/f"{prefix}_00_{suffix}.dat")
    return a

def interp(a,ell,col):
    return float(np.interp(ell,a[:,0],a[:,col]))

asdm=acoustic("sdmc"); alcd=acoustic("lcdm021")
print("CMB_DIAG_ACOUSTIC_SDMC",json.dumps(asdm),flush=True)
print("CMB_DIAG_ACOUSTIC_LCDM021",json.dumps(alcd),flush=True)
print("CMB_DIAG_ACOUSTIC_DELTA",json.dumps({k:asdm[k]-alcd[k] for k in asdm}),flush=True)

sl=cls("sdmc",True); ll=cls("lcdm021",True)
su=cls("sdmc",False); lu=cls("lcdm021",False)
# columns: ell, TT, EE, TE, BB, PP, TP, EP in CLASS D_l-style output
ells=[30,100,200,500,800,1000,1200,1500,1800,2000,2200,2500]
rows=[]
for ell in ells:
    r={"ell":ell}
    for tag,a,b in [("lensed",sl,ll),("unlensed",su,lu)]:
        for nm,c in [("TT",1),("EE",2),("TE",3),("PP",5)]:
            sv=interp(a,ell,c); lv=interp(b,ell,c)
            r[f"{tag}_{nm}_sdmc"]=sv; r[f"{tag}_{nm}_lcdm"]=lv
            r[f"{tag}_{nm}_ratio"]=sv/lv if abs(lv)>1e-30 else np.nan
    rows.append(r)
pd.DataFrame(rows).to_csv(OUT/"cmb_ratio_samples.csv",index=False)
print("CMB_DIAG_RATIO_SAMPLES",json.dumps(rows),flush=True)

# RMS fractional difference by broad high-l windows, TT/EE.
windows=[(30,400),(400,800),(800,1200),(1200,1800),(1800,2500)]
wr=[]
for lo,hi in windows:
    for tag,a,b in [("lensed",sl,ll),("unlensed",su,lu)]:
        for nm,c in [("TT",1),("EE",2)]:
            grid=np.arange(lo,hi+1)
            av=np.interp(grid,a[:,0],a[:,c]); bv=np.interp(grid,b[:,0],b[:,c])
            frac=(av-bv)/np.maximum(np.abs(bv),1e-30)
            wr.append(dict(lo=lo,hi=hi,type=tag,spectrum=nm,
                           mean_frac=float(np.mean(frac)),rms_frac=float(np.sqrt(np.mean(frac**2))),
                           max_abs_frac=float(np.max(np.abs(frac)))))
pd.DataFrame(wr).to_csv(OUT/"cmb_window_differences.csv",index=False)
print("CMB_DIAG_WINDOWS",json.dumps(wr),flush=True)

# Peak locations in TT D_l for lensed and unlensed, restricted to acoustic region.
def peaks(a,col=1):
    m=(a[:,0]>=100)&(a[:,0]<=2500)
    x=a[m,0]; y=a[m,col]
    inds,_=find_peaks(y,distance=180,prominence=max(np.ptp(y)*0.01,1e-30))
    # take ordered first 7
    return [{"ell":int(round(x[i])),"Dl":float(y[i])} for i in inds[:7]]
pout={
 "sdmc_lensed_TT":peaks(sl,1),"lcdm_lensed_TT":peaks(ll,1),
 "sdmc_unlensed_TT":peaks(su,1),"lcdm_unlensed_TT":peaks(lu,1),
 "sdmc_lensed_EE":peaks(sl,2),"lcdm_lensed_EE":peaks(ll,2)
}
Path(OUT/"cmb_peaks.json").write_text(json.dumps(pout,indent=2))
print("CMB_DIAG_PEAKS",json.dumps(pout),flush=True)
