#!/usr/bin/env python3
from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve

OUT=Path("output/negative_screen_sn"); OUT.mkdir(parents=True,exist_ok=True)
DATA=Path("sn_data")
OR=4.17998772e-5
D0=0.34231919445927034
DFLOOR=0.045
LAM=17.925
ZT=17.775

CANDIDATES=[
 dict(id="sobol049",H0=70.6661848725751,ob=0.02213832986768335,oc=0.12419066331610083,
      ns=0.9614044621046632,tau=0.05464880801877007,lnAs=3.051477308589965,
      AF=0.05709640212077647,zc=2.8532181778922676,width=0.575869657304138),
 dict(id="sobol009",H0=70.79871074482799,ob=0.02221026752647944,oc=0.12381796840205787,
      ns=0.9633178803175687,tau=0.0596811323207803,lnAs=3.0596968970261513,
      AF=0.05435267057083547,zc=3.1564399330876767,width=0.517175733782351),
]

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def ini(root,c):
    H0=c["H0"]; ob=c["ob"]; oc=c["oc"]; h=H0/100.
    OX=1.-(ob+oc+OR)/(h*h); As=math.exp(c["lnAs"])/1e10
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {ob:.15g}
    omega_cdm = {oc:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {As:.17e}
    n_s = {c["ns"]:.15g}
    tau_reio = {c["tau"]:.15g}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {c["AF"]:.17g},{c["zc"]:.17g},{c["width"]:.17g},{D0:.17g},1.0,{DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {OX:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    write background = yes
    write thermodynamics = no
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    output_verbose = 0
    """)

def lcdm_ini(root):
    return textwrap.dedent(f"""\
    H0 = 67.36
    omega_b = 0.02237
    omega_cdm = 0.1200
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = 2.10e-9
    n_s = 0.9649
    tau_reio = 0.0544
    write background = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    output_verbose = 0
    """)

def run(c):
    root=str(OUT/(c["id"]+"_")); ip=OUT/(c["id"]+".ini"); ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    if cp.returncode:
        print(cp.stdout[-3000:],flush=True); raise RuntimeError(c["id"])
    bg=table(root+"00_background.dat")
    stable=dict(
        min_D=float(bg["kin (D)"].min()),
        min_cs2=float(bg["c_s^2"].min()),
        max_cs2=float(bg["c_s^2"].max()),
        max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
    )
    print("NEG_SN_STABILITY",c["id"],stable,flush=True)
    if not(stable["min_D"]>0 and stable["min_cs2"]>0 and stable["max_cs2"]<=1):
        raise RuntimeError("stability")
    return bg

bgs={c["id"]:run(c) for c in CANDIDATES}
lr=str(OUT/"lcdm_"); li=OUT/"lcdm.ini"; li.write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(li)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-3000:])
bgs["lcdm"]=table(lr+"00_background.dat")

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov size {path}: {len(a)} vs {n*n}")
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))
def chi(pack,resid):
    cf,one,den=pack; y=cho_solve(cf,resid,check_finite=False)
    return float(resid.dot(y)-(one.dot(y))**2/den)
def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1.+z)
    return 5*np.log10((1.+zh)*(1.+z)*DA)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm)); mask=ppz>0.01
ppm,ppz,ppzh,C=ppm[mask],ppz[mask],ppzh[mask],C[np.ix_(mask,mask)]
pppack=setup(C)

p=DATA/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float)
unpack=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
depack=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

def scores(bg):
    return {
      "pantheonplus":chi(pppack,ppm-mu(bg,ppz,ppzh)),
      "union3":chi(unpack,unm-mu(bg,unz,unz)),
      "desy5":chi(depack,dem-mu(bg,dez,dezh)),
    }

ref=scores(bgs["lcdm"])
rows=[]
for key in ["sobol049","sobol009"]:
    s=scores(bgs[key])
    row={"model":key}
    for k in ref:
        row[f"{k}_chi2"]=s[k]; row[f"{k}_lcdm"]=ref[k]; row[f"delta_{k}"]=s[k]-ref[k]
    rows.append(row)
    print("NEG_SN_RESULT",json.dumps(row,sort_keys=True),flush=True)
pd.DataFrame(rows).to_csv(OUT/"negative_screen_sn.csv",index=False)
(OUT/"negative_screen_sn.json").write_text(json.dumps(rows,indent=2))
