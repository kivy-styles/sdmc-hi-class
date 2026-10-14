#!/usr/bin/env python3
"""
a010 -> AF11 background/structure corridor with a010 primordial sector frozen.

Purpose:
  Find a background point that preserves the new a010 exact-closure basin while
  moving H0 upward and improving the calibrated-ladder direction.  This is a
  screen only.  Full Plik and official DESI DR1 full shape must promote any
  finalist.

The exact covariance SN gates are included directly in the screen.  BAO, BBN,
age and compressed S8 weak-lensing diagnostics use the same Part IX battery
conventions.  JWST remains an independent nonlinear mapper coefficient g3 and
therefore does not alter the linear screen.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/a010_af11_corridor"); OUT.mkdir(parents=True,exist_ok=True)
DATA=Path("sn_data"); BAO=Path("bao_data")
TCMB=2.7255; CALSIG=.0025; OR=4.17998772e-5
# a010 primordial sector
AS=2.1131951739138815e-9; NS=0.9631153829544783; TAU=0.054600376411341134
OB=0.022083219194622913; OC=0.12299536722293603
D0=0.34231919445927034; POWER=1.0; LAM=18.40625
# endpoint 0 = g022/a010 exact background + structure
E0=dict(H0=69.71482083084993,AF=0.02048146490100771,zc=3.927876388467848,
        width=0.33114133956842123,Dfloor=0.05181542627513409,zt=16.19317622240633,
        Alate=0.024822032167576252,Blate=0.01264704344701022)
# endpoint 1 = AF11 background/structure, but keep a010 primordial sector
E1=dict(H0=69.85635133907199,AF=0.021666666666666667,zc=3.79610941009596,
        width=0.3229390993900597,Dfloor=0.046828035046346486,zt=15.982883202601224,
        Alate=0.018589897081255913,Blate=0.013815921754576268)

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[x for x in lines if x.startswith("#") and re.search(r"1\s*:",x)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        names.append(hdr[m.end():(ms[i+1].start() if i+1<len(ms) else len(hdr))].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:a=a[1:]
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(len(C)); u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one@u)
def prof(pack,r):
    cf,o,d=pack; y=cho_solve(cf,r,check_finite=False); return float(r@y-(o@y)**2/d)
# SN
pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
Cpp=cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm)); m=ppz>0.01
ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=setup(Cpp)
p=DATA/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float); unp=setup(cov(DATA/"Union3/mag_covmat.txt",len(unm)))
de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float); derr=de[cm["muerr_final"]].to_numpy(float)
dep=setup(cov(DATA/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

# BAO
rows=[]
for ln in (BAO/"desi_2024_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
    if ln.strip() and not ln.startswith("#"):
        z,v,q=ln.split(); rows.append((float(z),float(v),q))
bcov=np.loadtxt(BAO/"desi_2024_gaussian_bao_ALL_GCcomb_cov.txt"); binv=np.linalg.inv(bcov); bobs=np.array([x[1] for x in rows])

def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sn(bg):
    return dict(pp=prof(ppp,ppm-mu(bg,ppz,ppzh)),u3=prof(unp,unm-mu(bg,unz,unz)),d5=prof(dep,dem-mu(bg,dez,dezh)))
def rd(bg,th):
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); cr=[]
    for i in range(len(zz)-1):
      if (kb[i]-1)*(kb[i+1]-1)<=0 and kb[i]!=kb[i+1]:
        f=(1-kb[i])/(kb[i+1]-kb[i]); cr.append(zz[i]+f*(zz[i+1]-zz[i]))
    z=float(cr[0]); return float(np.interp(z,bg.z,bg["comov.snd.hrz."]))
def bao(bg,rdrag):
    pred=[]
    for z,_,q in rows:
      DM=float(np.interp(z,bg.z,bg["comov. dist."])); H=float(np.interp(z,bg.z,bg["H [1/Mpc]"])); DH=1/H
      DV=(z*DM*DM*DH)**(1/3); pred.append({"DM_over_rs":DM/rdrag,"DH_over_rs":DH/rdrag,"DV_over_rs":DV/rdrag}[q])
    d=np.array(pred)-bobs; return float(d@binv@d)
def cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n);ee=np.zeros(n);te=np.zeros(n);pp=np.zeros(n);tt[ell]=a[:,1]*conv;ee[ell]=a[:,2]*conv;te[ell]=a[:,3]*conv;L=ell.astype(float);pp[ell]=a[:,5]*L*(L+1)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),dict(tt=tt,ee=ee,te=te,pp=pp)
def planck(path):
    h,s=cls(path)
    def f(A):
      x=float(high.chi_squared(h,A_planck=A))-2*float(lowT.log_likelihood(s["tt"],calib=A))-2*float(lowE.log_likelihood(s["ee"],calib=A))
      lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
      x+=-2*float(lens.log_likelihood(s,**lp))+((A-1)/CALSIG)**2; return x
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded"); return float(op.fun)
def sig8(path):
    a=np.loadtxt(path); k=a[:,0];P=a[:,1];x=8*k;W=np.where(abs(x)<1e-5,1-x*x/10,3*(np.sin(x)-x*np.cos(x))/x**3)
    return float(np.sqrt(np.trapezoid(k**3*P/(2*np.pi**2)*W**2,np.log(k))))

def interp(t):
    return {k:(1-t)*E0[k]+t*E1[k] for k in E0}
def ini(root,c):
    h=c["H0"]/100.; ox=1-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={c['H0']}
omega_b={OB}
omega_cdm={OC}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={AS}
n_s={NS}
tau_reio={TAU}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={c['AF']},{c['zc']},{c['width']},{D0},{POWER},{c['Dfloor']}
expansion_model=sdmc_full
expansion_smg={ox},{LAM},{c['zt']},0.5,{c['Alate']},0.25,{c['Blate']},1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
modes=s
output=tCl,pCl,lCl,mPk
lensing=yes
l_max_scalars=3000
P_k_max_h/Mpc=2.0
z_pk=0
write background=yes
write thermodynamics=yes
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
out=[]
for i,t in enumerate(np.linspace(0,0.8,17)):
    c=interp(float(t)); tag=f"t{i:02d}"; root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,c))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    r=dict(id=tag,t=float(t),**c,status="FAIL")
    try:
      bg=tab(root+"00_background.dat"); th=tab(root+"00_thermodynamics.dat")
      r["Dmin"]=float(bg["kin (D)"].min());r["cs2min"]=float(bg["c_s^2"].min());r["cs2max"]=float(bg["c_s^2"].max())
      if cp.returncode==0 and r["Dmin"]>0 and r["cs2min"]>0 and r["cs2max"]<=1:
        ss=sn(bg); rr=rd(bg,th); bc=bao(bg,rr); pc=planck(root+"00_cl_lensed.dat")
        pk=sorted(OUT.glob(tag+"_*pk.dat"))[0]; s=sig8(pk); Om=float(np.interp(0,bg.z,bg["Omega_m(z)"])); S8=s*math.sqrt(Om/.3)
        age=float(np.interp(0,bg.z,bg["proper time [Gyr]"]))
        r.update(status="OK",planck_lite=pc,sn_pp=ss["pp"],sn_u3=ss["u3"],sn_d5=ss["d5"],rd=rr,bao=bc,sigma8=s,S8=S8,Omega_m0=Om,age=age,
                 shoes_proxy=((c["H0"]-73.49)/1.048)**2)
    except Exception as e:r["error"]=repr(e)
    print("A010_AF11_CORRIDOR",json.dumps(r,sort_keys=True),flush=True);out.append(r)
    for p in OUT.glob(tag+"_*"):
      try:p.unlink()
      except:pass
df=pd.DataFrame(out); df.to_csv(OUT/"corridor.csv",index=False)
ok=df[df.status=="OK"].copy()
base=ok.iloc[(ok.t-0).abs().argmin()]
for x in ["planck_lite","sn_pp","sn_u3","sn_d5","bao","S8","age","shoes_proxy"]:ok["d_"+x]=ok[x]-base[x]
# favor points with little Planck cost, no SN deterioration >0.25, and ladder gain
sel=ok[(ok.d_planck_lite<0.6)&(ok.d_sn_d5<0.35)&(ok.d_sn_pp<0.35)&(ok.d_sn_u3<0.35)].copy()
sel["score"]=sel.d_planck_lite+0.35*np.maximum(sel.d_sn_d5,0)+0.25*np.maximum(sel.d_sn_pp,0)+0.25*np.maximum(sel.d_sn_u3,0)+0.05*sel.d_bao+0.04*sel.d_shoes_proxy
sel=sel.sort_values(["score","shoes_proxy"]).head(8)
sel.to_csv(OUT/"promotion.csv",index=False)
(OUT/"summary.json").write_text(json.dumps({"baseline":base.to_dict(),"promotion":sel.to_dict("records")},indent=2))
print("A010_AF11_SUMMARY",json.dumps({"promotion":sel.to_dict("records")},sort_keys=True),flush=True)
