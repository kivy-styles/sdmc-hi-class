#!/usr/bin/env python3
from pathlib import Path
import math,re,subprocess,textwrap,json
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/planck_gap_diag"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=0.0025
H0=70.8514; ob=0.02239952; oc=0.12444227328918850
ns0=0.964; lnAs0=3.076; tau0=0.0544
OR=4.17998772e-5
AF=0.0715; ZC=5.; WIDTH=1.426; D0=0.34231919445927034
DFLOOR=0.045; LAMBDA=17.925; ZT=17.775

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages")
lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names).sort_values("z")

def load_cls(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1
    conv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[ell]=a[:,1]*conv; ee[ell]=a[:,2]*conv; te[ell]=a[:,3]*conv
    L=ell.astype(float); pp[ell]=a[:,5]*L*(L+1.)
    return np.column_stack([ell,tt[ell],te[ell],ee[ell]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def score(path):
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
    opt=minimize_scalar(lambda A:pieces(A)[-1],bounds=(0.97,1.03),method="bounded",
                        options={"xatol":1e-10})
    A=float(opt.x); p=pieces(A)
    return {"A_planck":A,"chi2_high":p[0],"chi2_lowT":p[1],"chi2_lowE":p[2],
            "chi2_lensing":p[3],"chi2_cal":p[4],"chi2_total":p[5]}

def lcdm_ini(root,extra=""):
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
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    lensing_verbose = 0
    output_verbose = 0
    {extra}
    """)

def sdmc_ini(root,lnAs=lnAs0,tau=tau0,ns=ns0,YHe="0.2453",
             minimal_nu=False,A_L=None):
    h=H0/100.
    onu=0.0
    if minimal_nu:
        onu=0.06/93.14
    Ox=1.-(ob+oc+OR+onu)/(h*h)
    As=math.exp(lnAs)/1e10
    nu=("N_ncdm = 1\nN_ur = 2.0308\nm_ncdm = 0.06" if minimal_nu
        else "N_ncdm = 0\nN_ur = 3.046")
    al="" if A_L is None else f"A_L = {A_L:.12g}"
    return textwrap.dedent(f"""\
    H0 = {H0}
    omega_b = {ob}
    omega_cdm = {oc}
    {nu}
    T_cmb = 2.7255
    YHe = {YHe}
    A_s = {As:.17e}
    n_s = {ns}
    tau_reio = {tau}

    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_independent_kinetic
    parameters_smg = {AF}, {ZC}, {WIDTH}, {D0}, 1.0, {DFLOOR}
    expansion_model = sdmc_full
    expansion_smg = {Ox:.17g},{LAMBDA},{ZT},0.5,0.01105624999,0.25,0.01951933685,1.5
    pert_initial_conditions_smg = zero
    method_qs_smg = fully_dynamic
    output_background_smg = 3

    modes = s
    output = tCl,pCl,lCl,mPk
    lensing = yes
    {al}
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

# Shared LCDM baseline.
lr=str(OUT/"lcdm_"); Path(OUT/"lcdm.ini").write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(OUT/"lcdm.ini")],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
if cp.returncode: raise RuntimeError(cp.stdout[-3000:])
base_lcdm=score(Path(lr+"00_cl_lensed.dat"))
print("PLANCK_GAP_LCDM",base_lcdm,flush=True)

variants=[dict(id="baseline")]
# Diagnostic cube only: this separates lensing smoothing (A_L),
# reionization damping (tau), and scalar amplitude (lnAs).
for al in [1.10,1.15,1.20]:
    for tv in [0.0584,0.0624,0.0664]:
        for lav in [3.056,3.066,3.076]:
            variants.append(dict(
                id=f"comb_AL{al:.2f}_t{tv:.4f}_A{lav:.3f}",
                A_L=al,tau=tv,lnAs=lav
            ))

rows=[]
for i,v in enumerate(variants,1):
    tag=v["id"]; root=str(OUT/(f"{i:02d}_{tag}_"))
    ip=OUT/(f"{i:02d}_{tag}.ini")
    kw=dict(v); kw.pop("id")
    ip.write_text(sdmc_ini(root,**kw))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                      text=True,timeout=240)
    rec={"id":tag,"returncode":cp.returncode,**kw}
    bgp=Path(root+"00_background.dat")
    if bgp.exists():
        bg=tab(bgp)
        rec.update(min_D=float(bg["kin (D)"].min()),
                   min_cs2=float(bg["c_s^2"].min()),
                   max_cs2=float(bg["c_s^2"].max()))
    clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and clp.exists():
        sc=score(clp); rec.update(sc)
        for k in ["chi2_high","chi2_lowT","chi2_lowE","chi2_lensing","chi2_cal","chi2_total"]:
            rec["delta_"+k]=sc[k]-base_lcdm[k]
        rec["stable_subluminal"]=bool(rec.get("min_D",0)>0 and rec.get("min_cs2",0)>0 and rec.get("max_cs2",2)<=1)
    else:
        rec["error"]=cp.stdout[-1200:].replace("\n"," | ")
    rows.append(rec)
    print("PLANCK_GAP_VARIANT",json.dumps(rec,sort_keys=True),flush=True)

df=pd.DataFrame(rows)
df.to_csv(OUT/"planck_gap_diagnostics.csv",index=False)
ok=df[(df.returncode==0)&(df.stable_subluminal==True)].copy()
print("PLANCK_GAP_BEST_TOTAL",ok.nsmallest(12,"chi2_total").to_dict("records"),flush=True)
print("PLANCK_GAP_BEST_HIGH",ok.nsmallest(12,"chi2_high").to_dict("records"),flush=True)
print("PLANCK_GAP_BEST_LENS",ok.nsmallest(12,"chi2_lensing").to_dict("records"),flush=True)
