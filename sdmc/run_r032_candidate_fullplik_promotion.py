#!/usr/bin/env python3
from pathlib import Path
import csv, json, math, subprocess, textwrap
import numpy as np
from iminuit import Minuit
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE import TTTEEE
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/fullplik_candidates"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255
OR=4.17998772e-5
AF=0.0715; ZC=5.; WIDTH=1.426
D0=0.34231919445927034; DFLOOR=.045
LAM=17.925; ZT=17.775

cands=[
 dict(label="CURRENT_NOSLIP",H0=70.8514,ob=0.02239952,oc=0.12444227328918850,
      ns=0.964,tau=0.0544,lnAs=3.076,A=0.0,zl=1.5,sl=.8),
 dict(label="ACOUSTIC_NOSLIP",H0=69.7095092787552,ob=0.02211778594478965,oc=0.1258596659367904,
      ns=0.9637720508351921,tau=0.06456347454525531,lnAs=3.0773358772153037,A=0.0,zl=1.5,sl=.8),
 dict(label="CURRENT_BRAID",H0=70.8514,ob=0.02239952,oc=0.12444227328918850,
      ns=0.964,tau=0.0544,lnAs=3.076,A=0.02,zl=1.5,sl=.8),
 dict(label="ACOUSTIC_BRAID",H0=69.7095092787552,ob=0.02211778594478965,oc=0.1258596659367904,
      ns=0.9637720508351921,tau=0.06456347454525531,lnAs=3.0773358772153037,A=0.02,zl=1.5,sl=.8),
]

def ini(c,root):
    h=c["H0"]/100.; Ox=1.-(c["ob"]+c["oc"]+OR)/(h*h); As=math.exp(c["lnAs"])/1e10
    return textwrap.dedent(f"""\
    H0 = {c['H0']:.15g}
    omega_b = {c['ob']:.15g}
    omega_cdm = {c['oc']:.15g}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {As:.17e}
    n_s = {c['ns']:.15g}
    tau_reio = {c['tau']:.15g}
    Omega_Lambda = 0
    Omega_fld = 3.1443554e-8
    fluid_equation_of_state = SDMC_TRACKER
    cs2_fld = 0.003
    use_ppf = no
    Omega_smg = -1
    gravity_model = sdmc_v3_lens_braiding
    parameters_smg = {AF},{ZC},{WIDTH},{D0:.17g},1.0,{DFLOOR},{c['A']:.17g},{c['zl']:.17g},{c['sl']:.17g}
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

paths={}
for c in cands:
    root=str(OUT/(c["label"].lower()+"_"))
    ip=OUT/(c["label"].lower()+".ini"); ip.write_text(ini(c,root))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    print("FULLPLIK_CLASS",c["label"],cp.returncode,flush=True)
    if cp.returncode:
        print(cp.stdout[-2500:],flush=True); raise SystemExit(2)
    paths[c["label"]]=OUT/(c["label"].lower()+"_00_cl_lensed.dat")

lr=str(OUT/"lcdm_"); li=OUT/"lcdm.ini"; li.write_text(lcdm_ini(lr))
cp=subprocess.run(["./class",str(li)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
if cp.returncode: raise RuntimeError(cp.stdout[-2500:])
paths["LCDM"]=OUT/"lcdm_00_cl_lensed.dat"

high=TTTEEE(packages_path='planck_full')
lowT=TT(packages_path='planck_packages')
lowE=EE(packages_path='planck_packages')
lens=LensingNative(packages_path='planck_packages')

fixed={
 'cib_index':-1.3,'galf_TE_index':-2.4,'galf_EE_index':-2.4,
 'A_sbpx_100_100_TT':1.0,'A_sbpx_143_143_TT':1.0,'A_sbpx_143_217_TT':1.0,'A_sbpx_217_217_TT':1.0,
 'galf_EE_A_100':0.055,'galf_EE_A_100_143':0.040,'galf_EE_A_100_217':0.094,'galf_EE_A_143':0.086,
 'galf_EE_A_143_217':0.21,'galf_EE_A_217':0.70,
 'A_cnoise_e2e_100_100_EE':1.0,'A_cnoise_e2e_143_143_EE':1.0,'A_cnoise_e2e_217_217_EE':1.0,
 'A_sbpx_100_100_EE':1.0,'A_sbpx_100_143_EE':1.0,'A_sbpx_100_217_EE':1.0,
 'A_sbpx_143_143_EE':1.0,'A_sbpx_143_217_EE':1.0,'A_sbpx_217_217_EE':1.0,
 'A_pol':1.0,'calib_100P':1.021,'calib_143P':0.966,'calib_217P':1.040,
}
names=[
 'A_planck','calib_100T','calib_217T',
 'A_cib_217','xi_sz_cib','A_sz','ksz_norm',
 'gal545_A_100','gal545_A_143','gal545_A_143_217','gal545_A_217',
 'ps_A_100_100','ps_A_143_143','ps_A_143_217','ps_A_217_217',
 'galf_TE_A_100','galf_TE_A_100_143','galf_TE_A_100_217',
 'galf_TE_A_143','galf_TE_A_143_217','galf_TE_A_217'
]
x0=np.array([1.0,1.0002,0.99805,67.,0.05,7.,3.,8.6,10.6,23.5,91.9,257.,47.,40.,104.,.130,.130,.46,.207,.69,1.938],float)
seed_cov=np.array([
 1.008245535114707,0.9998750001280258,0.9982119122897034,
 44.22091109885158,0.9997030504321224,6.692287855624427,0.00016370049661516226,
 8.45309812343869,10.60269070310368,20.195272036637586,97.0693107700243,
 252.7505102917758,57.67427029377452,63.72262893550619,129.8627383677383,
 .11563079474456282,.13456076987301357,.4804478052831501,
 .22690668463414007,.6690126964433923,2.1213120205986318
],float)
bounds=[
 (.97,1.03),(.9946,1.0058),(.99285,1.00325),(0,200),(0,1),(0,10),(0,10),
 (0,24.6),(0,26.6),(0,91.5),(0,251.9),(0,400),(0,400),(0,400),(0,400),
 (0,.466),(0,.418),(0,1.18),(0,.783),(0,1.41),(0,6.258)
]
gauss={
 'A_planck':(1.,.0025),'calib_100T':(1.0002,.0007),'calib_217T':(.99805,.00065),
 'gal545_A_100':(8.6,2.),'gal545_A_143':(10.6,2.),'gal545_A_143_217':(23.5,8.5),'gal545_A_217':(91.9,20.),
 'galf_TE_A_100':(.130,.042),'galf_TE_A_100_143':(.130,.036),'galf_TE_A_100_217':(.46,.09),
 'galf_TE_A_143':(.207,.072),'galf_TE_A_143_217':(.69,.09),'galf_TE_A_217':(1.938,.54)
}
expected=set(high.expected_params); supplied=set(fixed)|set(names)
if expected!=supplied: raise RuntimeError(f"Planck parameter mismatch {expected-supplied} {supplied-expected}")

def load(path):
    a=np.loadtxt(path); ell=a[:,0].astype(int); n=int(ell.max())+1; conv=(TCMB*1e6)**2
    dl={k:np.zeros(n) for k in ['tt','ee','bb','te','pp','tp','ep']}
    dl['tt'][ell]=a[:,1]*conv; dl['ee'][ell]=a[:,2]*conv; dl['te'][ell]=a[:,3]*conv; dl['bb'][ell]=a[:,4]*conv
    ll=ell.astype(float)*(ell.astype(float)+1.)
    dl['pp'][ell]=a[:,5]*ll
    dl['tp'][ell]=a[:,6]*np.sqrt(ll)*(TCMB*1e6)
    dl['ep'][ell]=a[:,7]*np.sqrt(ll)*(TCMB*1e6)
    cl={k:np.zeros(n) for k in ['tt','ee','bb','te']}
    fac=np.zeros(n); q=ell>=2; fac[ell[q]]=2*np.pi/ll[q]
    for k in cl: cl[k][ell]=dl[k][ell]*fac[ell]
    return cl,dl

def prior_chi2(p):
    c=0.
    for n,(mu,sig) in gauss.items(): c+=((p[n]-mu)/sig)**2
    c+=((p['ksz_norm']+1.6*p['A_sz']-9.5)/3.)**2
    return float(c)

def make_obj(path):
    cl,dl=load(path)
    def objective(*x):
        p=dict(fixed); p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p)); A=p['A_planck']
        lt=float(lowT.log_likelihood(dl['tt'],calib=A)); le=float(lowE.log_likelihood(dl['ee'],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        if not all(np.isfinite([lh,lt,le,ll])): return 1e100
        return float(-2*(lh+lt+le+ll)+prior_chi2(p))
    def pieces(x):
        p=dict(fixed); p.update(dict(zip(names,map(float,x))))
        lh=float(high.log_likelihood(cl,**p)); A=p['A_planck']
        lt=float(lowT.log_likelihood(dl['tt'],calib=A)); le=float(lowE.log_likelihood(dl['ee'],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,'calibration_param',None) else {}
        ll=float(lens.log_likelihood(dl,**lp))
        return p,lh,lt,le,ll
    return objective,pieces

def profile(label,path,seeds):
    fun,pieces=make_obj(path); best=None
    for iseed,seed in enumerate(seeds):
        m=Minuit(fun,*seed,name=names); m.errordef=1.; m.strategy=1; m.tol=.1
        for n,b in zip(names,bounds): m.limits[n]=b
        for n in names: m.errors[n]=max(1e-5,.02*max(abs(m.values[n]),1.))
        m.migrad(ncall=6000)
        if not m.fmin.is_valid:
            m.simplex(ncall=3000); m.migrad(ncall=6000)
        row=(float(m.fval),np.array([m.values[n] for n in names],float),bool(m.fmin.is_valid),int(m.nfcn))
        print("CAND_FULLPLIK_SEED",label,iseed,*row[::2],row[3],flush=True)
        if best is None or row[0]<best[0]: best=row
    fval,x,valid,nfcn=best; p,lh,lt,le,ll=pieces(x)
    row=dict(model=label,chi2_profile_total=fval,chi2_high_full=-2*lh,chi2_lowT=-2*lt,
             chi2_lowE=-2*le,chi2_lensing=-2*ll,chi2_nuisance_priors=prior_chi2(p),
             A_planck=p['A_planck'],valid=valid,nfcn=nfcn)
    row.update({n:p[n] for n in names if n!='A_planck'})
    print("CAND_FULLPLIK_RESULT",json.dumps(row,sort_keys=True),flush=True)
    return row,x

results=[]; bestseed=seed_cov.copy()
for label in ["CURRENT_NOSLIP","ACOUSTIC_NOSLIP","CURRENT_BRAID","ACOUSTIC_BRAID","LCDM"]:
    seeds=[bestseed,x0,seed_cov]
    row,x=profile(label,paths[label],seeds)
    results.append(row); bestseed=x

lcdm=[r for r in results if r["model"]=="LCDM"][0]
for r in results:
    r["delta_vs_lcdm"]=r["chi2_profile_total"]-lcdm["chi2_profile_total"]
    print("CAND_FULLPLIK_DELTA",r["model"],r["delta_vs_lcdm"],flush=True)

pd_fields=list(results[0].keys())
with (OUT/"candidate_fullplik_profile.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=pd_fields); w.writeheader(); w.writerows(results)
