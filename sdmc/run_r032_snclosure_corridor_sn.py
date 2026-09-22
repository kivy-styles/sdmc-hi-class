#!/usr/bin/env python3
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np, pandas as pd
from scipy.linalg import cho_factor, cho_solve

OUT=Path("output/snclosure_corridor_sn"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
OR=4.17998772e-5
OB=0.022083219194622913; OC=0.12299536722293603
AS=2.117602412937069e-9; NS=0.9625227132590487; TAU=0.055202901571989066
AF=0.03200255395658314; ZC=4.007449422683567; WIDTH=0.34799921004101636
D0=0.34231919445927034; DF=0.04607192634791136; LAM=18.40625; ZT=16.173189924377947
A0=0.01105624999; B0=0.01951933685; H00=70.79187943335671
A1=0.018589897081255913; B1=0.013815921754576268; H01=69.85635133907199

LOCAL_SN={"pantheonplus":1406.2147033223882,"union3":28.783140002196888,"desy5":1650.6353947147727}

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov mismatch {path}")
    return a.reshape(n,n)

def pack(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0]); u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))
def prof(p,r):
    cf,one,den=p; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(SNROOT/'PantheonPlus/Pantheon+SH0ES.dat',sep=r'\s+',header=0,engine='python')
ppm=pp['m_b_corr'].to_numpy(float); ppz=pp['zHD'].to_numpy(float); ppzh=pp['zHEL'].to_numpy(float)
Cpp=read_cov(SNROOT/'PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov',len(ppm))
m=ppz>0.01; ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=pack(Cpp)
p=SNROOT/'Union3/lcparam_full.txt'; cols=p.read_text().splitlines()[0].removeprefix('#').split()
un=pd.read_csv(p,sep=r'\s+',comment='#',names=cols,engine='python'); cm={c.lower():c for c in un.columns}
unm=un[cm['mb']].to_numpy(float); unz=un[cm['zcmb']].to_numpy(float); unp=pack(read_cov(SNROOT/'Union3/mag_covmat.txt',len(unm)))
de=pd.read_csv(SNROOT/'DESY5/DES-SN5YR_HD.csv'); cm={c.lower():c for c in de.columns}
dem=de[cm['mu']].to_numpy(float); dez=de[cm['zhd']].to_numpy(float); dezh=de[cm['zhel']].to_numpy(float); derr=de[cm['muerr_final']].to_numpy(float)
dep=pack(read_cov(SNROOT/'DESY5/covsys_000.txt',len(dem))+np.diag(derr*derr))

def mu(bg,z,zh):
    zz=bg.z.to_numpy(); dm=bg['comov. dist.'].to_numpy(); o=np.argsort(zz)
    DM=np.interp(z,zz[o],dm[o]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
    return {"pantheonplus":prof(ppp,ppm-mu(bg,ppz,ppzh)),
            "union3":prof(unp,unm-mu(bg,unz,unz)),
            "desy5":prof(dep,dem-mu(bg,dez,dezh))}

def ini(root,A,B,H0):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={AS:.17e}
n_s={NS:.17g}
tau_reio={TAU:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={AF},{ZC},{WIDTH},{D0},1.0,{DF}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},0.5,{A:.17g},0.25,{B:.17g},1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
write background=yes
root={root}
format=class
input_verbose=0
background_verbose=0
thermodynamics_verbose=0
perturbations_verbose=0
spectra_verbose=0
output_verbose=0
""")

rows=[]
for t in [0.0,0.25,0.5,0.75,1.0]:
    A=A0+t*(A1-A0); B=B0+t*(B1-B0); H0=H00+t*(H01-H00)
    tag=f"t{int(round(100*t)):03d}"; root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,A,B,H0))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=180)
    rec={"t":t,"A":A,"B":B,"H0":H0,"status":"FAIL"}
    bgp=Path(root+"00_background.dat")
    if cp.returncode==0 and bgp.exists():
        bg=table(bgp); s=sns(bg)
        rec.update({f"chi2_{k}":v for k,v in s.items()})
        rec.update({f"delta_{k}":v-LOCAL_SN[k] for k,v in s.items()})
        rec["status"]="OK"
    else: rec["error"]=cp.stdout[-800:].replace("\n"," | ")
    rows.append(rec); print("CORRIDOR_SN_POINT",json.dumps(rec,sort_keys=True),flush=True)

pd.DataFrame(rows).to_csv(OUT/"corridor_sn.csv",index=False)
(OUT/"corridor_sn.json").write_text(json.dumps(rows,indent=2))
