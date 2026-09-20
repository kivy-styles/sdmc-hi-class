#!/usr/bin/env python3
"""
Stage-1 six-dimensional ordinary-cosmology screen for the frozen r032 SDMC
structural point (lambda_e,z_t,Dfloor)=(17.925,17.775,0.045).

This is a deterministic screening surrogate:
  * SDMC uses the validated independently parameterized kinetic closure.
  * Planck uses Plik-lite + lowT + lowE + lensing with A_planck profiled.
  * Pantheon+, Union3 and DES-Y5 are scored separately (never summed together).
  * Every point is required to satisfy the native global D>0, c_s^2>0,
    c_s^2<=1 stability gate before it can enter a shortlist.
Final candidates must be replayed with the exact covariant reconstruction,
official full Plik nuisance profiling and the raw DESI DR1 full-shape likelihood.
"""

from pathlib import Path
import csv, json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/ordinary6d")
OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255
CAL_SIGMA=0.0025

# Frozen SDMC structural point.
AF=0.0715
ZC=5.00
WIDTH=1.426
D0=0.34231919445927034
POWER=1.0
DFLOOR=0.045
LAMBDA_E=17.925
ZT=17.775
OMEGA_R_PHYS=4.17998772e-5

# Current exact-leader ordinary cosmology and the Planck-control reference.
CURRENT=np.array([
    70.8514,
    0.02239952,
    0.12444227328918850,
    0.964,
    3.076,
    0.0544,
],float)
REFERENCE=np.array([
    67.36,
    0.02237,
    0.1200,
    0.9649,
    math.log(1e10*2.10e-9),
    0.0544,
],float)
NAMES=["H0","omega_b","omega_cdm","n_s","ln10As","tau_reio"]

# Conservative local envelope around the union of current and Planck-control
# coordinates. These are search bounds, not priors or inferred intervals.
MARGIN=np.array([0.55,0.00050,0.0040,0.015,0.040,0.012],float)
LOW=np.minimum(CURRENT,REFERENCE)-MARGIN
HIGH=np.maximum(CURRENT,REFERENCE)+MARGIN

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr))
    names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path)
    ell=a[:,0].astype(int)
    n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv
    ee[ell]=a[:,2]*conv
    te[ell]=a[:,3]*conv
    L=ell.astype(float)
    pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def planck_score(path):
    harr,dls=load_cls(path)
    def pieces(A):
        A=float(A)
        ch=float(high.chi_squared(harr,A_planck=A))
        ct=float(-2*lowT.log_likelihood(dls["tt"],calib=A))
        ce=float(-2*lowE.log_likelihood(dls["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        cl=float(-2*lens.log_likelihood(dls,**lp))
        cp=((A-1.)/CAL_SIGMA)**2
        return ch,ct,ce,cl,cp,ch+ct+ce+cl+cp
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),
                        method="bounded",options={"xatol":1e-10})
    A=float(opt.x)
    return (A,)+pieces(A)

# --- SN data ---
DATA=Path("sn_data")
def read_cov(path,n):
    a=np.asarray(np.loadtxt(path),float).reshape(-1)
    if len(a)==n*n+1 and int(round(a[0]))==n:
        a=a[1:]
    if len(a)!=n*n:
        raise RuntimeError(f"covariance size mismatch {path}: {len(a)} vs {n*n}")
    return a.reshape(n,n)

def project_setup(C):
    cf=cho_factor(C,lower=True,check_finite=False)
    one=np.ones(C.shape[0])
    u=cho_solve(cf,one,check_finite=False)
    return cf,one,float(one.dot(u))

def projected_chi2(pack,resid):
    cf,one,denom=pack
    y=cho_solve(cf,resid,check_finite=False)
    return float(resid.dot(y)-(one.dot(y))**2/denom)

pp=pd.read_csv(DATA/"PantheonPlus/Pantheon+SH0ES.dat",sep=r"\s+",header=0,engine="python")
pp_mag=pp["m_b_corr"].to_numpy(float)
pp_z=pp["zHD"].to_numpy(float)
pp_zh=pp["zHEL"].to_numpy(float)
pp_C=read_cov(DATA/"PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",len(pp_mag))
mask=pp_z>0.01
pp_mag,pp_z,pp_zh=pp_mag[mask],pp_z[mask],pp_zh[mask]
pp_C=pp_C[np.ix_(mask,mask)]
pp_pack=project_setup(pp_C)

p=DATA/"Union3/lcparam_full.txt"
lines=p.read_text().splitlines()
cols=lines[0].removeprefix("#").split()
un=pd.read_csv(p,sep=r"\s+",comment="#",names=cols,engine="python")
cmap={c.lower():c for c in un.columns}
un_mag=un[cmap["mb"]].to_numpy(float)
un_z=un[cmap["zcmb"]].to_numpy(float)
un_zh=un_z.copy()
un_C=read_cov(DATA/"Union3/mag_covmat.txt",len(un_mag))
un_pack=project_setup(un_C)

de=pd.read_csv(DATA/"DESY5/DES-SN5YR_HD.csv")
cmap={c.lower():c for c in de.columns}
de_mag=de[cmap["mu"]].to_numpy(float)
de_z=de[cmap["zhd"]].to_numpy(float)
de_zh=de[cmap["zhel"]].to_numpy(float)
de_err=de[cmap["muerr_final"]].to_numpy(float)
de_C=read_cov(DATA/"DESY5/covsys_000.txt",len(de_mag))+np.diag(de_err*de_err)
de_pack=project_setup(de_C)

def sn_mu(bg,z,zh):
    DM=np.interp(z,bg.z,bg["comov. dist."])
    DA=DM/(1.+z)
    return 5.*np.log10((1.+zh)*(1.+z)*DA)

def sn_scores(bg):
    return (
        projected_chi2(pp_pack,pp_mag-sn_mu(bg,pp_z,pp_zh)),
        projected_chi2(un_pack,un_mag-sn_mu(bg,un_z,un_zh)),
        projected_chi2(de_pack,de_mag-sn_mu(bg,de_z,de_zh)),
    )

def sigma8_from_pk(path,h):
    a=np.loadtxt(path)
    kh=a[:,0]; P=a[:,1]
    q=(kh>1e-4)&(P>0)&np.isfinite(P)
    kh=kh[q]; P=P[q]
    # CLASS P(k) table uses k in h/Mpc and P in (Mpc/h)^3, so R=8 Mpc/h.
    x=kh*8.0
    W=np.ones_like(x)
    m=x!=0
    W[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    y=kh**3*P*W**2/(2*np.pi**2)
    return float(np.sqrt(np.trapezoid(y,x=np.log(kh))))

def acoustic(bg,th):
    zstar=float(th.loc[th["g [Mpc^-1]"].idxmax(),"z"])
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy()
    cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i])
            cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross:
        return {}
    zdrag=float(cross[0])
    rs=float(np.interp(zstar,bg.z,bg["comov.snd.hrz."]))
    rd=float(np.interp(zdrag,bg.z,bg["comov.snd.hrz."]))
    dm=float(np.interp(zstar,bg.z,bg["comov. dist."]))
    return {"z_star":zstar,"z_drag":zdrag,"rs_star":rs,"rd":rd,"ell_A":math.pi*dm/rs}

def sdmc_ini(x,root):
    H0,ob,oc,ns,lnAs,tau=map(float,x)
    h=H0/100.
    Ox=1.-(ob+oc+OMEGA_R_PHYS)/(h*h)
    As=math.exp(lnAs)/1e10
    return textwrap.dedent(f"""\
    H0 = {H0:.15g}
    omega_b = {ob:.15g}
    omega_cdm = {oc:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {As:.17e}
    n_s = {ns:.15g}
    tau_reio = {tau:.15g}

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF}, {ZC}, {WIDTH}, {D0:.17g}, {POWER}, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {Ox:.17g}, {LAMBDA_E}, {ZT}, 0.5, 0.01105624999, 0.25, 0.01951933685, 1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic

    output_background_smg = 3
    modes = s
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
    write background = yes
    write thermodynamics = yes
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
    output = tCl,pCl,lCl,mPk
    lensing = yes
    l_max_scalars = 3000
    P_k_max_h/Mpc = 1.0
    z_pk = 0
    write background = yes
    write thermodynamics = yes
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

# Fixed LCDM reference for screen deltas.
lroot=str(OUT/"lcdm_")
lip=OUT/"lcdm.ini"
lip.write_text(lcdm_ini(lroot))
cp=subprocess.run(["./class",str(lip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode!=0:
    raise RuntimeError("LCDM control failed\n"+cp.stdout[-3000:])
lp=planck_score(Path(lroot+"00_cl_lensed.dat"))
lbg=table(Path(lroot+"00_background.dat"))
lsn=sn_scores(lbg)
lcdm={"planck":lp[-1],"pantheonplus":lsn[0],"union3":lsn[1],"desy5":lsn[2]}
print("ORD6D_LCDM_REFERENCE",lcdm,flush=True)

# Deterministic design: exact current, Planck-reference ordinary coordinates,
# three points on the connecting line, plus a 64-point Sobol space-filling set.
anchors=[
    ("current",CURRENT),
    ("reference_coords",REFERENCE),
    ("line25",CURRENT+0.25*(REFERENCE-CURRENT)),
    ("line50",CURRENT+0.50*(REFERENCE-CURRENT)),
    ("line75",CURRENT+0.75*(REFERENCE-CURRENT)),
]
# Fast correlated line diagnostic only.
design=anchors

meta={
    "names":NAMES,
    "current":CURRENT.tolist(),
    "reference":REFERENCE.tolist(),
    "low":LOW.tolist(),
    "high":HIGH.tolist(),
    "structural":{"AF":AF,"zc":ZC,"width":WIDTH,"D0":D0,"power":POWER,
                  "Dfloor":DFLOOR,"lambda_e":LAMBDA_E,"z_t":ZT},
    "n_design":len(design),
    "note":"screen only; final points require exact covariant replay + full Plik + raw DESI FS",
}
(OUT/"ordinary6d_design.json").write_text(json.dumps(meta,indent=2))

rows=[]
for ic,(tag,x) in enumerate(design,1):
    H0,ob,oc,ns,lnAs,tau=map(float,x)
    h=H0/100.
    Ox=1.-(ob+oc+OMEGA_R_PHYS)/(h*h)
    root=str(OUT/(f"{ic:03d}_{tag}_"))
    ip=OUT/(f"{ic:03d}_{tag}.ini")
    ip.write_text(sdmc_ini(x,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec={
        "id":tag,"index":ic,"H0":H0,"omega_b":ob,"omega_cdm":oc,"n_s":ns,
        "ln10As":lnAs,"A_s":math.exp(lnAs)/1e10,"tau_reio":tau,
        "Omega_X0":Ox,"Omega_m0":(ob+oc)/(h*h),
        "returncode":cp.returncode,"status":"FAIL",
    }
    bgp=Path(root+"00_background.dat")
    clp=Path(root+"00_cl_lensed.dat")
    thp=Path(root+"00_thermodynamics.dat")
    pkp=Path(root+"00_z1_pk.dat")
    if bgp.exists():
        bg=table(bgp)
        rec.update({
            "min_D_all":float(bg["kin (D)"].min()),
            "min_cs2_all":float(bg["c_s^2"].min()),
            "max_cs2_all":float(bg["c_s^2"].max()),
            "max_abs_noslip_all":float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))),
            "F0":float(bg.iloc[np.argmin(np.abs(bg.z.to_numpy()))]["M*^2_smg"]),
        })
    if cp.returncode==0 and bgp.exists() and clp.exists() and thp.exists() and pkp.exists():
        bg=table(bgp)
        th=table(thp)
        stable=(rec["min_D_all"]>0 and rec["min_cs2_all"]>0 and rec["max_cs2_all"]<=1.0)
        rec["stable_subluminal"]=bool(stable)
        if stable:
            ps=planck_score(clp)
            sn=sn_scores(bg)
            rec.update({
                "status":"OK","A_planck":ps[0],
                "chi2_high":ps[1],"chi2_lowT":ps[2],"chi2_lowE":ps[3],
                "chi2_lensing":ps[4],"chi2_cal":ps[5],"chi2_planck":ps[6],
                "chi2_pantheonplus":sn[0],"chi2_union3":sn[1],"chi2_desy5":sn[2],
                "delta_planck":ps[6]-lcdm["planck"],
                "delta_pantheonplus":sn[0]-lcdm["pantheonplus"],
                "delta_union3":sn[1]-lcdm["union3"],
                "delta_desy5":sn[2]-lcdm["desy5"],
                "sigma8":sigma8_from_pk(pkp,h),
            })
            rec.update(acoustic(bg,th))
            rec["screen_planck_pp"]=rec["delta_planck"]+rec["delta_pantheonplus"]
            rec["screen_planck_union3"]=rec["delta_planck"]+rec["delta_union3"]
            rec["screen_planck_desy5"]=rec["delta_planck"]+rec["delta_desy5"]
    if rec["status"]!="OK":
        rec["error"]=cp.stdout[-1400:].replace("\n"," | ")
    print("ORD6D_POINT",json.dumps(rec,sort_keys=True),flush=True)
    rows.append(rec)

df=pd.DataFrame(rows)
df.to_csv(OUT/"ordinary6d_screen.csv",index=False)
ok=df[(df.status=="OK") & (df.stable_subluminal==True)].copy()
if ok.empty:
    raise SystemExit("No stable subluminal candidates in 6D screen")

keep=set(["current","reference_coords","line25","line50","line75"])
for col in ["delta_planck","screen_planck_pp","screen_planck_union3","screen_planck_desy5"]:
    for x in ok.nsmallest(8,col)["id"]:
        keep.add(x)
short=ok[ok.id.isin(sorted(keep))].copy()
# Ranking is only for computational shortlisting. The final scientific result
# is based on full Plik + raw DESI, not these screen columns.
short["screen_best_any"]=short[
    ["delta_planck","screen_planck_pp","screen_planck_union3","screen_planck_desy5"]
].min(axis=1)
short=short.sort_values("screen_best_any")
short.to_csv(OUT/"ordinary6d_shortlist.csv",index=False)

print("ORD6D_BEST_PLANCK",ok.nsmallest(10,"delta_planck")[NAMES+["id","delta_planck","delta_pantheonplus","delta_union3","delta_desy5","sigma8","min_cs2_all","max_cs2_all"]].to_dict("records"),flush=True)
print("ORD6D_BEST_PP",ok.nsmallest(10,"screen_planck_pp")[NAMES+["id","screen_planck_pp","delta_planck","delta_pantheonplus"]].to_dict("records"),flush=True)
print("ORD6D_BEST_UNION3",ok.nsmallest(10,"screen_planck_union3")[NAMES+["id","screen_planck_union3","delta_planck","delta_union3"]].to_dict("records"),flush=True)
print("ORD6D_BEST_DESY5",ok.nsmallest(10,"screen_planck_desy5")[NAMES+["id","screen_planck_desy5","delta_planck","delta_desy5"]].to_dict("records"),flush=True)
print("ORD6D_SHORTLIST",short[NAMES+["id","delta_planck","delta_pantheonplus","delta_union3","delta_desy5","sigma8","rd","ell_A","min_cs2_all","max_cs2_all"]].to_dict("records"),flush=True)
