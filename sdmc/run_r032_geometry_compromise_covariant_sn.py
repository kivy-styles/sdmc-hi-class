#!/usr/bin/env python3
from pathlib import Path
import json, re
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve

OUT=Path("output/covariant_geometry_sn")
OUT.mkdir(parents=True,exist_ok=True)
COV=Path("input_cov")
DATA=Path("sn_data")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    return a.reshape(n,n)

def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False)
    one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def pchi(pack,r):
    cf,one,den=pack
    y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppmag=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppmag))
m=ppz>0.01
ppmag,ppz,ppzh,C=ppmag[m],ppz[m],ppzh[m],C[np.ix_(m,m)]
pppack=setup(C)

p=DATA/"Union3/lcparam_full.txt"
cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
cm={c.lower():c for c in un.columns}
unmag=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unzh=unz.copy()
unpack=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(unmag)))

de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv")
cm={c.lower():c for c in de.columns}
demag=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
deerr=de[cm["muerr_final"]].to_numpy(float)
depack=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(demag))+np.diag(deerr*deerr))

def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."])
    DA=DM/(1.+z)
    return 5.*np.log10((1.+zh)*(1.+z)*DA)

def score(bg):
    return dict(
      PantheonPlus=pchi(pppack,ppmag-mu(bg,ppz,ppzh)),
      Union3=pchi(unpack,unmag-mu(bg,unz,unzh)),
      DESY5=pchi(depack,demag-mu(bg,dez,dezh))
    )

cov=table(COV/"covariant_observable_00_background.dat")
lcdm=table(OUT/"lcdm_00_background.dat")
cs=score(cov); ls=score(lcdm)
delta={k:cs[k]-ls[k] for k in cs}
planck=5.36115841529
desi=-11.781497184073828
joint_PD=planck+desi

summary={
  "covariant_sn_chi2":cs,
  "lcdm_sn_chi2":ls,
  "delta_sn_covariant_minus_lcdm":delta,
  "delta_planck":planck,
  "delta_desi":desi,
  "delta_planck_plus_desi":joint_PD,
  "delta_joint_with_sn":{k:joint_PD+delta[k] for k in delta}
}
(OUT/"covariant_geometry_sn_summary.json").write_text(json.dumps(summary,indent=2))
pd.DataFrame([{
    "compilation":k,
    "chi2_covariant":cs[k],
    "chi2_lcdm":ls[k],
    "delta_sn":delta[k],
    "delta_planck_plus_desi":joint_PD,
    "delta_joint":joint_PD+delta[k]
} for k in cs]).to_csv(OUT/"covariant_geometry_sn.csv",index=False)
print("COV_GEOM_SN_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
