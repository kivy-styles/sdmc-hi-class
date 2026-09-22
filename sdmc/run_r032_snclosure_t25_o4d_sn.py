#!/usr/bin/env python3
"""
SN-aware ordinary-parameter re-optimization at the t=0.25 late-geometry corridor point.

Fixed:
  late geometry A=0.012939661762813977, B=0.018093483076144067
  edge024 structural action
Re-open:
  omega_b, omega_cdm, n_s, Q=ln(1e10 A_s)-2 tau
For every candidate H0 is re-solved to the acoustic target.  Exact SN
covariances are evaluated on the resulting background together with
Planck TTTEEE-lite + lowT + lowE + lensing.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from scipy.linalg import cho_factor,cho_solve
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure_t25_o4d_sn"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data")
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5
TARGET=301.6798271349195

H00=70.55799740978553
OB0=0.022083219194622913
OC0=0.12299536722293603
NS0=0.9625227132590487
TAU=0.055202901571989066
Q0=2.942463801247999

AF=0.03200255395658314
ZC=4.007449422683567
WIDTH=0.34799921004101636
D0=0.34231919445927034
DF=0.04607192634791136
LAM=18.40625
ZT=16.173189924377947
A_LATE=0.012939661762813977
B_LATE=0.018093483076144067

LOW=np.array([0.02190,0.1214,0.9595,2.9400],float)
HIGH=np.array([0.02232,0.1240,0.9665,2.9455],float)

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def derived(bg,th):
    b=bg.sort_values("z")
    imax=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zstar=float(th.iloc[imax]["z"])
    DM=float(np.interp(zstar,b["z"],b["comov. dist."]))
    rs=float(np.interp(zstar,b["z"],b["comov.snd.hrz."]))
    return dict(z_star=zstar,D_M_star=DM,r_s_star=rs,ell_A=float(np.pi*DM/rs))

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    harr,dls=load_cls(path)
    def pc(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    op=minimize_scalar(lambda A:pc(float(A))[-1],bounds=(.97,1.03),method="bounded",
                       options={"xatol":1e-10})
    A=float(op.x); p=pc(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"cov mismatch {path}")
    return a.reshape(n,n)
def cpack(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))
def prof(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(SNROOT/'PantheonPlus/Pantheon+SH0ES.dat',sep=r'\s+',header=0,engine='python')
ppm=pp['m_b_corr'].to_numpy(float); ppz=pp['zHD'].to_numpy(float); ppzh=pp['zHEL'].to_numpy(float)
Cpp=read_cov(SNROOT/'PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov',len(ppm))
m=ppz>0.01; ppm,ppz,ppzh,Cpp=ppm[m],ppz[m],ppzh[m],Cpp[np.ix_(m,m)]; ppp=cpack(Cpp)

p=SNROOT/'Union3/lcparam_full.txt'; cols=p.read_text().splitlines()[0].removeprefix('#').split()
un=pd.read_csv(p,sep=r'\s+',comment='#',names=cols,engine='python'); cm={c.lower():c for c in un.columns}
unm=un[cm['mb']].to_numpy(float); unz=un[cm['zcmb']].to_numpy(float)
unp=cpack(read_cov(SNROOT/'Union3/mag_covmat.txt',len(unm)))

de=pd.read_csv(SNROOT/'DESY5/DES-SN5YR_HD.csv'); cm={c.lower():c for c in de.columns}
dem=de[cm['mu']].to_numpy(float); dez=de[cm['zhd']].to_numpy(float); dezh=de[cm['zhel']].to_numpy(float)
derr=de[cm['muerr_final']].to_numpy(float)
dep=cpack(read_cov(SNROOT/'DESY5/covsys_000.txt',len(dem))+np.diag(derr*derr))

def mu_bg(bg,z,zh):
    b=bg.sort_values("z"); DM=np.interp(z,b["z"],b["comov. dist."]); DA=DM/(1+z)
    return 5*np.log10((1+zh)*(1+z)*DA)
def sns(bg):
    return dict(pantheonplus=prof(ppp,ppm-mu_bg(bg,ppz,ppzh)),
                union3=prof(unp,unm-mu_bg(bg,unz,unz)),
                desy5=prof(dep,dem-mu_bg(bg,dez,dezh)))

def ini(root,H0,ob,oc,ns,Q,full):
    As=math.exp(Q+2*TAU)/1e10
    h=H0/100.; ox=1.-(ob+oc+OR)/(h*h)
    obs="""modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
""" if full else ""
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={ob:.17g}
omega_cdm={oc:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={As:.17e}
n_s={ns:.17g}
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
expansion_smg={ox:.17g},{LAM},{ZT},0.5,{A_LATE},0.25,{B_LATE},1.5
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
{obs}
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

def run_class(tag,H0,ob,oc,ns,Q,full):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,H0,ob,oc,ns,Q,full))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    if cp.returncode: raise RuntimeError(cp.stdout[-900:])
    bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat")
    d=derived(bg,th)
    st=dict(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
            max_cs2=float(bg["c_s^2"].max()),
            max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    return root,bg,d,st

def candidate(tag,ob,oc,ns,Q):
    vals=[]
    for H in np.linspace(H00-0.8,H00+0.8,7):
        try:
            _,_,d,st=run_class(f"{tag}_r_{H:.4f}".replace(".","p"),float(H),ob,oc,ns,Q,False)
            if st["min_D"]>0 and st["min_cs2"]>0 and st["max_cs2"]<=1:
                vals.append((float(H),float(d["ell_A"]-TARGET)))
        except Exception: pass
    vals=sorted(vals); pair=None
    for a,b in zip(vals[:-1],vals[1:]):
        if a[1]==0 or a[1]*b[1]<=0: pair=(a,b); break
    if pair is None:
        rec=dict(id=tag,omega_b=ob,omega_cdm=oc,n_s=ns,Q=Q,status="NO_ROOT")
        print("T25O4D_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
    (h1,g1),(h2,g2)=pair
    H0=float(h1-g1*(h2-h1)/(g2-g1))
    try:
        root,bg,d,st=run_class(tag,H0,ob,oc,ns,Q,True)
    except Exception as e:
        rec=dict(id=tag,omega_b=ob,omega_cdm=oc,n_s=ns,Q=Q,H0=H0,status="FAIL",error=str(e)[-500:])
        print("T25O4D_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec
    stable=st["min_D"]>0 and st["min_cs2"]>0 and st["max_cs2"]<=1
    rec=dict(id=tag,omega_b=float(ob),omega_cdm=float(oc),n_s=float(ns),Q=float(Q),
             A_s=float(math.exp(Q+2*TAU)/1e10),H0=H0,delta_ellA=d["ell_A"]-TARGET,
             stable_subluminal=bool(stable),**d,**st,status="OK" if stable else "UNSTABLE")
    if stable:
        rec.update(pscore(root+"00_cl_lensed.dat"))
        ss=sns(bg)
        rec.update({f"sn_{k}":float(v) for k,v in ss.items()})
    print("T25O4D_POINT",json.dumps(rec,sort_keys=True),flush=True); return rec

rows=[]
base=candidate("center",OB0,OC0,NS0,Q0)
if base.get("status")!="OK": raise RuntimeError(base)
rows.append(base); P0=float(base["chi2_planck"])

anchors=[
 ("obm",OB0-0.00012,OC0,NS0,Q0),("obp",OB0+0.00012,OC0,NS0,Q0),
 ("ocm",OB0,OC0-0.0007,NS0,Q0),("ocp",OB0,OC0+0.0007,NS0,Q0),
 ("nsm",OB0,OC0,NS0-0.0015,Q0),("nsp",OB0,OC0,NS0+0.0015,Q0),
 ("qm",OB0,OC0,NS0,Q0-0.0012),("qp",OB0,OC0,NS0,Q0+0.0012),
]
for a in anchors: rows.append(candidate(*a))

sam=qmc.Sobol(d=4,scramble=True,seed=250925)
pts=qmc.scale(sam.random_base2(m=5),LOW,HIGH)
for i,p in enumerate(pts): rows.append(candidate(f"sobol{i:03d}",*map(float,p)))

for r in rows:
    if r.get("status")=="OK":
        r["delta_planck_vs_center"]=r["chi2_planck"]-P0
        r["delta_pp_vs_center"]=r["sn_pantheonplus"]-base["sn_pantheonplus"]
        r["delta_u3_vs_center"]=r["sn_union3"]-base["sn_union3"]
        r["delta_d5_vs_center"]=r["sn_desy5"]-base["sn_desy5"]
        r["second_sn_delta_vs_center"]=sorted([r["delta_pp_vs_center"],r["delta_u3_vs_center"],r["delta_d5_vs_center"]])[1]
        r["screen_joint_score"]=r["delta_planck_vs_center"]+r["second_sn_delta_vs_center"]

df=pd.DataFrame(rows); df.to_csv(OUT/"snclosure_t25_o4d_sn.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
bestP=ok.nsmallest(10,"chi2_planck").to_dict("records")
bestJ=ok.nsmallest(10,"screen_joint_score").to_dict("records")
summary={"baseline":base,"n_total":int(len(df)),"n_stable":int(len(ok)),
         "best_planck":bestP,"best_joint_screen":bestJ}
(OUT/"snclosure_t25_o4d_sn_summary.json").write_text(json.dumps(summary,indent=2))
print("T25O4D_SUMMARY",json.dumps(summary,sort_keys=True),flush=True)
