#!/usr/bin/env python3
from pathlib import Path
import json,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.linalg import cho_factor,cho_solve

OUT=Path("output/fulljoint_exactbo_sn"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
LOCAL={"pantheonplus":1406.2147033223882,"union3":28.783140002196888,"desy5":1650.6353947147727}
OR=4.17998772e-5
D0=0.34231919445927034
LAM=18.40625; ZT=16.173189924377947; TAU=0.055202901571989066

cands=[
{"id":"exactbo000","A":0.010977504386967048,"B":0.012857573034386149,"H0":69.81677547364065,
 "ob":0.02208511574370519,"oc":0.12298509428428875,"ns":0.962555355281866,"As":2.1197359302086e-9,
 "AF":0.023604633340554453,"ZC":3.4328776987879466,"W":0.36879721635160295,"DF":0.050275813910968366},
{"id":"exactbo001","A":0.010161639647297556,"B":0.013017204296421257,"H0":69.84073210322437,
 "ob":0.022084954853937286,"oc":0.1229859657675159,"ns":0.9625525861642491,"As":2.1197346699560444e-9,
 "AF":0.023604633340554453,"ZC":3.4328776987879466,"W":0.36879721635160295,"DF":0.050275813910968366},
{"id":"exactbo004","A":0.010557210430167612,"B":0.012939807320889084,"H0":69.82911676766862,
 "ob":0.022085032861097482,"oc":0.12298554323019364,"ns":0.96255392876673,"As":2.1197352809875867e-9,
 "AF":0.023632917504307768,"ZC":3.4401668971099166,"W":0.36856316167028097,"DF":0.05010044754926387}
]

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
    if len(a)!=n*n: raise RuntimeError((path,len(a),n*n))
    return a.reshape(n,n)

def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one@u)

def chi(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r@y-(one@y)**2/den)

pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm))
m=ppz>0.01; ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=setup(Cpp)

p=SNROOT/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float)
unp=setup(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

def mu(bg,z,zh):
    zz=bg.z.to_numpy(); dm=bg["comov. dist."].to_numpy(); o=np.argsort(zz)
    DM=np.interp(z,zz[o],dm[o]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)

def ini(c,root):
    h=c["H0"]/100.; ox=1-(c["ob"]+c["oc"]+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={c['H0']:.17g}
omega_b={c['ob']:.17g}
omega_cdm={c['oc']:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={c['As']:.17e}
n_s={c['ns']:.17g}
tau_reio={TAU}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={c['AF']},{c['ZC']},{c['W']},{D0},1.0,{c['DF']}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},0.5,{c['A']},0.25,{c['B']},1.5
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
lensing_verbose=0
output_verbose=0
""")

rows=[]
for c in cands:
    root=str(OUT/(c["id"]+"_")); ip=OUT/(c["id"]+".ini"); ip.write_text(ini(c,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=180)
    if cp.returncode: raise RuntimeError(c["id"]+" "+cp.stdout[-1000:])
    bg=table(root+"00_background.dat")
    scores={
      "pantheonplus":chi(ppp,ppm-mu(bg,ppz,ppzh)),
      "union3":chi(unp,unm-mu(bg,unz,unz)),
      "desy5":chi(dep,dem-mu(bg,dez,dezh))
    }
    r={"id":c["id"],**scores}
    for k,v in scores.items(): r["delta_"+k]=v-LOCAL[k]
    rows.append(r)
    print("FULLJOINT_SN_EXACT",json.dumps(r,sort_keys=True),flush=True)

pd.DataFrame(rows).to_csv(OUT/"sn_exact.csv",index=False)
(OUT/"sn_exact.json").write_text(json.dumps(rows,indent=2))
