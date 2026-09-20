#!/usr/bin/env python3
"""
Local 6D response/Hessian map around the frozen r032 exact-leader ordinary
cosmology.  Purpose: detect narrow compensated directions missed by the broad
Sobol screen.

Coordinates y are measured in the step units below.  We evaluate the center,
+/- one step on each axis, and all four corners for every parameter pair,
giving 73 points.  A full quadratic is then fit for Planck and each alternative
Planck+SN objective.  Newton candidates are evaluated only as diagnostics;
the final scientific gate remains exact covariant + full Plik + raw DESI.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap,itertools
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor,cho_solve
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/ordinary6d_local"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
NAMES=["H0","omega_b","omega_cdm","n_s","ln10As","tau_reio"]
X0=np.array([70.8514,0.02239952,0.12444227328918850,0.964,3.076,0.0544],float)
STEP=np.array([0.15,0.00012,0.0008,0.003,0.010,0.004],float)

AF=0.0715; ZC=5.; WIDTH=1.426; D0=0.34231919445927034; DFLOOR=0.045
LAM=17.925; ZT=17.775; OR=4.17998772e-5

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr); names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def planck_score(path):
    harr,dls=load_cls(path)
    def pcs(A):
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp)); cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda a:pcs(float(a))[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x); p=pcs(A)
    return dict(A_planck=A,chi2_high=p[0],chi2_lowT=p[1],chi2_lowE=p[2],
                chi2_lensing=p[3],chi2_cal=p[4],chi2_planck=p[5])

# SN setup
DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n: a=a[1:]
    if len(a)!=n*n: raise RuntimeError(f"bad covariance {path}")
    return a.reshape(n,n)
def setup(C):
    cf=cho_factor(C,lower=True,check_finite=False); one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False); return cf,one,float(one.dot(u))
def proj(pack,r):
    cf,one,den=pack; y=cho_solve(cf,r,check_finite=False)
    return float(r.dot(y)-(one.dot(y))**2/den)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
pm=pp["m_b_corr"].to_numpy(float); pz=pp["zHD"].to_numpy(float); pzh=pp["zHEL"].to_numpy(float)
C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(pm)); mm=pz>0.01
pm,pz,pzh,C=pm[mm],pz[mm],pzh[mm],C[np.ix_(mm,mm)]; ppack=setup(C)
p=DATA/"Union3/lcparam_full.txt"; lines=p.read_text().splitlines(); cols=lines[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python"); cm={c.lower():c for c in un.columns}
um=un[cm["mb"]].to_numpy(float); uz=un[cm["zcmb"]].to_numpy(float); upack=setup(read_cov(DATA/"Union3/mag_covmat.txt",len(um)))
de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv"); cm={c.lower():c for c in de.columns}
dm=de[cm["mu"]].to_numpy(float); dz=de[cm["zhd"]].to_numpy(float); dzh=de[cm["zhel"]].to_numpy(float)
derr=de[cm["muerr_final"]].to_numpy(float)
dpack=setup(read_cov(DATA/"DESY5/covsys_000.txt",len(dm))+np.diag(derr*derr))
def mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."]); DA=DM/(1.+z)
    return 5*np.log10((1.+zh)*(1.+z)*DA)
def sn_scores(bg):
    return proj(ppack,pm-mu(bg,pz,pzh)),proj(upack,um-mu(bg,uz,uz)),proj(dpack,dm-mu(bg,dz,dzh))

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
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
    write background = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    lensing_verbose = 0
    output_verbose = 0
    """)

def sdmc_ini(x,root):
    H0,ob,oc,ns,lnAs,tau=map(float,x); h=H0/100.; Ox=1.-(ob+oc+OR)/(h*h)
    return textwrap.dedent(f"""\
    H0 = {H0:.16g}
    omega_b = {ob:.16g}
    omega_cdm = {oc:.16g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {math.exp(lnAs)/1e10:.17e}
    n_s = {ns:.16g}
    tau_reio = {tau:.16g}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF}, {ZC}, {WIDTH}, {D0:.17g}, 1.0, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {Ox:.17g},{LAM},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl
    lensing = yes
    l_max_scalars = 3000
    write background = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    lensing_verbose = 0
    output_verbose = 0
    """)

# Fixed reference deltas.
lr=str(OUT/"lcdm_"); (OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-2000:])
lps=planck_score(Path(lr+"00_cl_lensed.dat")); lsn=sn_scores(table(Path(lr+"00_background.dat")))
BASE=dict(planck=lps["chi2_planck"],pp=lsn[0],union3=lsn[1],desy5=lsn[2])
print("LOCAL6D_LCDM",BASE,flush=True)

# Build normalized design.
design=[("center",np.zeros(6))]
for i in range(6):
    for s in (-1.,1.):
        y=np.zeros(6); y[i]=s; design.append((f"{NAMES[i]}_{s:+.0f}",y))
for i,j in itertools.combinations(range(6),2):
    for si,sj in itertools.product((-1.,1.),repeat=2):
        y=np.zeros(6); y[i]=si; y[j]=sj
        design.append((f"{NAMES[i]}{si:+.0f}_{NAMES[j]}{sj:+.0f}",y))
assert len(design)==73

rows=[]
for idx,(tag,y) in enumerate(design,1):
    x=X0+STEP*y; root=str(OUT/(f"{idx:03d}_{tag}_")); ip=OUT/(f"{idx:03d}_{tag}.ini")
    ip.write_text(sdmc_ini(x,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"id":tag,"index":idx,**{f"y_{n}":float(v) for n,v in zip(NAMES,y)},
         **{n:float(v) for n,v in zip(NAMES,x)},"returncode":cp.returncode,"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if bgp.exists():
        bg=table(bgp); rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),
                                  max_cs2=float(bg["c_s^2"].max()))
    if cp.returncode==0 and bgp.exists() and clp.exists():
        stable=rec.get("min_D",0)>0 and rec.get("min_cs2",0)>0 and rec.get("max_cs2",2)<=1
        rec["stable_subluminal"]=bool(stable)
        if stable:
            ps=planck_score(clp); sn=sn_scores(table(bgp)); rec.update(ps)
            rec.update(delta_planck=ps["chi2_planck"]-BASE["planck"],
                       delta_pp=sn[0]-BASE["pp"],delta_union3=sn[1]-BASE["union3"],delta_desy5=sn[2]-BASE["desy5"])
            rec["obj_planck"]=rec["delta_planck"]
            rec["obj_pp"]=rec["delta_planck"]+rec["delta_pp"]
            rec["obj_union3"]=rec["delta_planck"]+rec["delta_union3"]
            rec["obj_desy5"]=rec["delta_planck"]+rec["delta_desy5"]
            rec["status"]="OK"
    if rec["status"]!="OK": rec["error"]=cp.stdout[-900:].replace("\n"," | ")
    rows.append(rec); print("LOCAL6D_POINT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows); df.to_csv(OUT/"local6d_response.csv",index=False)
ok=df[df.status=="OK"].copy()
if len(ok)<50: raise SystemExit(f"too few stable local points: {len(ok)}")

# Quadratic feature map q=c+g.y+0.5 Hii yi^2+Hij yi yj
def features(Y):
    Y=np.atleast_2d(Y); cols=[np.ones(len(Y))]
    cols += [Y[:,i] for i in range(6)]
    cols += [0.5*Y[:,i]**2 for i in range(6)]
    cols += [Y[:,i]*Y[:,j] for i,j in itertools.combinations(range(6),2)]
    return np.column_stack(cols)

Y=ok[[f"y_{n}" for n in NAMES]].to_numpy(float)
X=features(Y)
summaries={}
candidates=[]
for obj in ["obj_planck","obj_pp","obj_union3","obj_desy5"]:
    z=ok[obj].to_numpy(float)
    beta=np.linalg.lstsq(X,z,rcond=None)[0]
    c=beta[0]; g=beta[1:7]; diag=beta[7:13]; cross=beta[13:]
    H=np.diag(diag)
    for val,(i,j) in zip(cross,itertools.combinations(range(6),2)):
        H[i,j]=H[j,i]=val
    eig=np.linalg.eigvalsh(H)
    try: ystar=-np.linalg.solve(H,g)
    except np.linalg.LinAlgError: ystar=np.full(6,np.nan)
    pred=float(c+g@ystar+0.5*ystar@H@ystar) if np.all(np.isfinite(ystar)) else np.nan
    summaries[obj]={"center_fit":float(c),"gradient":g.tolist(),"hessian":H.tolist(),
                    "eigenvalues":eig.tolist(),"newton_y":ystar.tolist(),"newton_x":(X0+STEP*ystar).tolist(),
                    "predicted_obj":pred}
    print("LOCAL6D_QUADRATIC",obj,json.dumps(summaries[obj]),flush=True)
    if np.all(np.isfinite(ystar)) and np.max(np.abs(ystar))<=2.0:
        candidates.append((obj,ystar))

Path(OUT/"local6d_quadratic.json").write_text(json.dumps({"names":NAMES,"x0":X0.tolist(),"step":STEP.tolist(),"fits":summaries},indent=2))

# Evaluate unique Newton candidates that stay within 2 local steps.
seen=[]
cres=[]
for obj,y in candidates:
    if any(np.max(np.abs(y-y0))<1e-6 for _,y0 in seen): continue
    seen.append((obj,y.copy()))
    x=X0+STEP*y; tag="newton_"+obj
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(sdmc_ini(x,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
    rec={"id":tag,"source_objective":obj,**{n:float(v) for n,v in zip(NAMES,x)},"returncode":cp.returncode}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp); rec.update(min_D=float(bg["kin (D)"].min()),min_cs2=float(bg["c_s^2"].min()),max_cs2=float(bg["c_s^2"].max()))
        if rec["min_D"]>0 and rec["min_cs2"]>0 and rec["max_cs2"]<=1:
            ps=planck_score(clp); sn=sn_scores(bg); rec.update(ps)
            rec.update(delta_planck=ps["chi2_planck"]-BASE["planck"],delta_pp=sn[0]-BASE["pp"],
                       delta_union3=sn[1]-BASE["union3"],delta_desy5=sn[2]-BASE["desy5"])
            rec["obj_planck"]=rec["delta_planck"]; rec["obj_pp"]=rec["delta_planck"]+rec["delta_pp"]
            rec["obj_union3"]=rec["delta_planck"]+rec["delta_union3"]; rec["obj_desy5"]=rec["delta_planck"]+rec["delta_desy5"]
    cres.append(rec); print("LOCAL6D_NEWTON",json.dumps(rec,sort_keys=True),flush=True)
pd.DataFrame(cres).to_csv(OUT/"local6d_newton_candidates.csv",index=False)
