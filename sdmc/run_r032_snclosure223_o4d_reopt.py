#!/usr/bin/env python3
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

OUT=Path("output/snclosure223_o4d"); OUT.mkdir(parents=True,exist_ok=True)
SNROOT=Path("sn_data"); BAOROOT=Path("bao_data")
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5
TARGET=301.6798271349195

# SN-friendly late geometry found by the first closure screen.
AL=0.018589897081255913; BL=0.013815921754576268
TAUA=.25; TAUB=1.5

# Successful edge024 structural anchor.
AF=0.03200255395658314; ZC=4.007449422683567
WIDTH=0.34799921004101636; D0=0.34231919445927034
DF=0.04607192634791136; LAM=18.40625; ZT=16.173189924377947
DN=.5; TAU=0.055202901571989066

OB0=0.022083219194622913; OC0=0.12299536722293603
NS0=0.9625227132590487; Q0=2.942463801247999
H00=69.85635133907199

EDGE_P_LITE=1022.3837030216768
EDGE_PD_FAIR=-1.016724593277559
EDGE_BAO=15.5323571850153
LOCAL_BAO=13.400290718544086
LOCAL_SN={"pantheonplus":1406.2147033223882,
          "union3":28.783140002196888,
          "desy5":1650.6353947147727}

LOW=np.array([.02190,.1215,.9595,2.9395])
HIGH=np.array([.02230,.1245,.9665,2.9458])

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
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False); return cf,one,float(one@u)

def chprof(p,r):
    cf,one,den=p; y=cho_solve(cf,r,check_finite=False)
    return float(r@y-(one@y)**2/den)

# SN covariance products
pp=pd.read_csv(SNROOT/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
ppm=pp["m_b_corr"].to_numpy(float); ppz=pp["zHD"].to_numpy(float); ppzh=pp["zHEL"].to_numpy(float)
C=read_cov(SNROOT/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(ppm))
m=ppz>.01; ppm,ppz,ppzh,C=ppm[m],ppz[m],ppzh[m],C[np.ix_(m,m)]; ppp=pack(C)

p=SNROOT/"Union3/lcparam_full.txt"; cols=p.read_text().splitlines()[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
unm=un[cm["mb"]].to_numpy(float); unz=un[cm["zcmb"]].to_numpy(float)
unp=pack(read_cov(SNROOT/"Union3/mag_covmat.txt",len(unm)))

de=pd.read_csv(SNROOT/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dem=de[cm["mu"]].to_numpy(float); dez=de[cm["zhd"]].to_numpy(float); dezh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dep=pack(read_cov(SNROOT/"DESY5/covsys_000.txt",len(dem))+np.diag(derr*derr))

# DESI Gaussian BAO
br=[]
for ln in (BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
    if not ln.strip() or ln.startswith("#"): continue
    z,v,q=ln.split(); br.append((float(z),float(v),q))
bo=np.array([x[1] for x in br])
bi=np.linalg.inv(np.loadtxt(BAOROOT/"desi_2024_gaussian_bao_ALL_GCcomb_cov.txt"))

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def derived(bg,th):
    b=bg.sort_values("z"); im=int(np.argmax(th["g [Mpc^-1]"].to_numpy()))
    zs=float(th.iloc[im]["z"]); DM=float(np.interp(zs,b.z,b["comov. dist."]))
    rs=float(np.interp(zs,b.z,b["comov.snd.hrz."]))
    return zs,DM,rs,float(np.pi*DM/rs)

def sn(bg):
    b=bg.sort_values("z")
    def mu(z,zh):
        DM=np.interp(z,b.z,b["comov. dist."]); DA=DM/(1+z)
        return 5*np.log10((1+zh)*(1+z)*DA)
    return {"pantheonplus":chprof(ppp,ppm-mu(ppz,ppzh)),
            "union3":chprof(unp,unm-mu(unz,unz)),
            "desy5":chprof(dep,dem-mu(dez,dezh))}

def rd(bg,th):
    z=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); xs=[]
    for i in range(len(z)-1):
        if (kb[i]-1)*(kb[i+1]-1)<=0 and kb[i]!=kb[i+1]:
            t=(1-kb[i])/(kb[i+1]-kb[i]); xs.append(z[i]+t*(z[i+1]-z[i]))
    if not xs: raise RuntimeError("no drag crossing")
    zd=float(xs[0]); b=bg.sort_values("z")
    return float(np.interp(zd,b.z,b["comov.snd.hrz."])),zd

def bao(bg,rdrag):
    b=bg.sort_values("z"); pred=[]
    for z,_,q in br:
        DM=float(np.interp(z,b.z,b["comov. dist."])); H=float(np.interp(z,b.z,b["H [1/Mpc]"]))
        DH=1/H; DV=(z*DM*DM*DH)**(1/3)
        pred.append({"DM_over_rs":DM/rdrag,"DH_over_rs":DH/rdrag,"DV_over_rs":DV/rdrag}[q])
    d=np.array(pred)-bo; return float(d@bi@d)

def cls(path):
    a=np.loadtxt(path); e=a[:,0].astype(int); n=int(e.max())+1; cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[e]=a[:,1]*cv; ee[e]=a[:,2]*cv; te[e]=a[:,3]*cv
    L=e.astype(float); pp[e]=a[:,5]*L*(L+1)
    return np.column_stack([e,tt[e],te[e],ee[e]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    ha,d=cls(path)
    def f(A):
        x=float(high.chi_squared(ha,A_planck=A))
        x+=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        x+=float(-2*lens.log_likelihood(d,**lp))+((A-1)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def ini(root,H0,ob,oc,ns,Q,full):
    As=math.exp(Q+2*TAU)/1e10; h=H0/100; ox=1-(ob+oc+OR)/(h*h)
    obs="modes=s\noutput=tCl,pCl,lCl\nlensing=yes\nl_max_scalars=3000\n" if full else ""
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
expansion_smg={ox:.17g},{LAM},{ZT},{DN},{AL},{TAUA},{BL},{TAUB}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
{obs}write background=yes
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

def run(tag,H0,ob,oc,ns,Q,full):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,H0,ob,oc,ns,Q,full))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    if cp.returncode: raise RuntimeError(cp.stdout[-800:])
    bg=table(root+"00_background.dat"); th=table(root+"00_thermodynamics.dat")
    st=(float(bg["kin (D)"].min()),float(bg["c_s^2"].min()),float(bg["c_s^2"].max()),
        float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
    return root,bg,th,derived(bg,th),st

def cand(tag,ob,oc,ns,Q):
    vals=[]
    for H in (H00-.45,H00+.45):
        try:
            _,_,_,d,st=run(tag+("a" if H<H00 else "b"),H,ob,oc,ns,Q,False)
            if st[0]>0 and st[1]>0 and st[2]<=1: vals.append((H,d[3]))
        except: pass
    if len(vals)<2: return {"id":tag,"status":"ROOT_FAIL"}
    der=(vals[1][1]-vals[0][1])/(vals[1][0]-vals[0][0])
    H=float(np.clip(vals[0][0]+(TARGET-vals[0][1])/der,H00-1.5,H00+1.5))
    try: root,bg,th,d,st=run(tag,H,ob,oc,ns,Q,True)
    except Exception as e: return {"id":tag,"status":"FINAL_FAIL","error":str(e)[-300:]}
    stable=st[0]>0 and st[1]>0 and st[2]<=1
    r={"id":tag,"status":"OK" if stable else "UNSTABLE","omega_b":ob,"omega_cdm":oc,"n_s":ns,"Q":Q,
       "A_s":math.exp(Q+2*TAU)/1e10,"H0":H,"min_D":st[0],"min_cs2":st[1],"max_cs2":st[2],
       "max_abs_noslip":st[3],"ellA":d[3],"delta_ellA":d[3]-TARGET}
    if not stable: return r
    pc,Ap=pscore(root+"00_cl_lensed.dat"); ss=sn(bg); rg,zd=rd(bg,th); bc=bao(bg,rg)
    r.update(chi2_planck=pc,A_planck=Ap,rd=rg,zdrag=zd,bao_chi2=bc)
    r["delta_planck_vs_edge024"]=pc-EDGE_P_LITE
    # Raw-DESI movement is unknown here; use BAO as a conservative 1/4 proxy only for ranking.
    r["pd_proxy"]=EDGE_PD_FAIR+r["delta_planck_vs_edge024"]+.25*(bc-EDGE_BAO)
    js=[]
    for k,v in ss.items():
        dv=v-LOCAL_SN[k]; r["sn_"+k]=v; r["sn_delta_"+k]=dv
        j=r["pd_proxy"]+dv; r["joint_"+k+"_proxy"]=j; js.append(j)
    r["second_joint_proxy"]=sorted(js)[1]
    r["n_sn_closed_proxy"]=sum(x<0 for x in js)
    r["goal_score"]=max(r["pd_proxy"],r["second_joint_proxy"])
    return r

pts=[("center",OB0,OC0,NS0,Q0),
("obm",OB0-.00010,OC0,NS0,Q0),("obp",OB0+.00010,OC0,NS0,Q0),
("ocm",OB0,OC0-.0007,NS0,Q0),("ocp",OB0,OC0+.0007,NS0,Q0),
("nsm",OB0,OC0,NS0-.0012,Q0),("nsp",OB0,OC0,NS0+.0012,Q0),
("qm",OB0,OC0,NS0,Q0-.001),("qp",OB0,OC0,NS0,Q0+.001)]
sam=qmc.Sobol(d=4,scramble=True,seed=22309)
for i,p in enumerate(qmc.scale(sam.random_base2(m=5),LOW,HIGH)):
    pts.append((f"sobol{i:03d}",*map(float,p)))
rows=[]
for p in pts:
    r=cand(*p); rows.append(r); print("SN223_POINT",json.dumps(r,sort_keys=True),flush=True)
df=pd.DataFrame(rows); df.to_csv(OUT/"snclosure223_o4d.csv",index=False)
ok=df[df.status=="OK"].sort_values(["goal_score","second_joint_proxy","pd_proxy"])
best=ok.head(10).to_dict("records")
summary={"late_geometry":{"A":AL,"B":BL,"tauA":TAUA,"tauB":TAUB},
         "n_ok":int(len(ok)),"best":best,"target":"pd_proxy<0 and second_joint_proxy<0"}
(OUT/"snclosure223_o4d_summary.json").write_text(json.dumps(summary,indent=2))
print("SN223_BEST",json.dumps(summary,sort_keys=True),flush=True)
