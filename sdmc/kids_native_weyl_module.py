#!/usr/bin/env python3
from pathlib import Path
import re, numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.integrate import cumulative_trapezoid
from astropy.io import fits
from cosmosis.datablock import option_section

ELL=np.geomspace(2.,1.e4,900)

def zhdr(p):
    for line in Path(p).read_text().splitlines()[:10]:
        m=re.search(r"(?:at\s+)?(?:redshift\s+)?z\s*=\s*([0-9eE+\-.]+)",line,re.I)
        if m:return float(m.group(1))
    raise RuntimeError(str(p))

def cube(directory,prefix,h,kind):
    pats=f"{prefix}*pk_weyl.dat" if kind=="weyl" else f"{prefix}*pk.dat"
    rows=[]
    for p in sorted(Path(directory).glob(pats)):
        if kind!="weyl" and ("weyl" in p.name or "_nl_" in p.name): continue
        a=np.loadtxt(p); z=zhdr(p); k=a[:,0]*h
        y=a[:,1]*h if kind=="weyl" else a[:,1]/h**3
        rows.append((z,k,y))
    rows.sort(key=lambda x:x[0]); z=np.array([x[0] for x in rows]); k=rows[0][1]; v=np.array([x[2] for x in rows])
    return RegularGridInterpolator((z,np.log(k)),np.log(v),bounds_error=False,fill_value=np.nan),z

def tail(z,y):
    return (-cumulative_trapezoid(y[::-1],z[::-1],initial=0.))[::-1]

def setup(options):
    model=options.get_string(option_section,"model")
    specdir=options.get_string(option_section,"spectra_dir")
    data=options.get_string(option_section,"data_file")
    if model=="late300":
        h=0.6971482083084993; Om=(0.022083219194622913+0.12299536722293603)/h**2
    else:
        h=0.6856859; Om=(0.02240637+0.11824151)/h**2
    Iq,zq=cube(specdir,model+"_",h,"weyl")
    Ip,zp=cube(specdir,model+"_",h,"matter")
    # background
    bp=Path(specdir)/(model+"_background.dat")
    lines=bp.read_text().splitlines(); hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    aa=np.loadtxt(bp); d={n:aa[:,i] for i,n in enumerate(names)}; o=np.argsort(d["z"])
    with fits.open(data) as f:
        nz=f["NZ_SOURCE"].data
        zs=np.asarray(nz["Z_MID"],float); ns=np.array([np.asarray(nz[f"BIN{i}"],float) for i in range(1,7)])
    return dict(model=model,h=h,Om=Om,Iq=Iq,Ip=Ip,zq=zq,zp=zp,bgz=np.asarray(d["z"])[o],H=np.asarray(d["H [1/Mpc]"])[o],
                chi=np.asarray(d[next(n for n in names if "comov" in n.lower() and "dist" in n.lower())])[o],zs=zs,ns=ns)

def execute(block,c):
    zmax=min(c["zq"].max(),c["zp"].max(),float(c["zs"].max()))
    z=np.linspace(0.01,zmax,360); H=np.interp(z,c["bgz"],c["H"]); chi=np.interp(z,c["bgz"],c["chi"])
    nz=[]
    for a in c["ns"]:
        y=np.maximum(interp1d(c["zs"],a,bounds_error=False,fill_value=0.)(z),0.); y/=np.trapezoid(y,z); nz.append(y)
    nz=np.asarray(nz); g=[]
    for y in nz:
        g.append(tail(z,y)-chi*tail(z,y/np.maximum(chi,1e-12)))
    g=np.asarray(g)
    zm=np.broadcast_to(z[None,:],(len(ELL),len(z))); km=(ELL[:,None]+.5)/np.maximum(chi[None,:],1e-12)
    pts=np.column_stack([zm.ravel(),np.log(km.ravel())]); Q=np.exp(c["Iq"](pts)).reshape(km.shape); Q=np.where(np.isfinite(Q),Q,0.)
    block["shear_cl","ell"]=ELL
    block["shear_cl","nbin"]=6
    for i in range(6):
        for j in range(i,6):
            cl=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
            block["shear_cl",f"bin_{j+1}_{i+1}"]=cl
    return 0
def cleanup(c): pass
