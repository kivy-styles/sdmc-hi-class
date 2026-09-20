#!/usr/bin/env python3
from pathlib import Path
import sys,re,json,math,subprocess,textwrap
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize

sys.path.insert(0,'desi-kp-cosmological-likelihoods/dr1/cobaya')
import lsstypes as types
from cosmoprimo import PowerSpectrumInterpolator1D,PowerSpectrumBAOFilter,Cosmology
from cosmoprimo.fiducial import DESI
from velocileptors.EPT.ept_fullresum_varyDz_nu_fftw import REPT
from desi_fs_bao_all import list_zrange,dataset_fn,get_tracer_label,get_physical_stochastic_settings

OUT=Path('output/desi_candidates'); OUT.mkdir(parents=True,exist_ok=True)
C=299792.458
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
lcdm=dict(label="LCDM",H0=67.36,ob=0.02237,oc=0.1200,ns=0.9649,tau=0.0544,lnAs=math.log(1e10*2.10e-9))

# Official DESI blocks and effective redshifts
blocks=[]
records=[]
for tracer,iz,zrange in list_zrange:
    if 'lya' in tracer.lower(): continue
    fn=dataset_fn('desi_likelihood',tracer,zrange,observable_name='spectrum-poles-rotated')
    d=types.read(fn)
    sp=d.observable.get('spectrum')
    theory=d.window.theory.get('spectrum')
    kin=theory.get(0).coords('k')
    win=d.window.at.observable.get('spectrum').at.theory.get('spectrum').value()
    rec=dict(tracer=tracer,iz=iz,zrange=list(zrange),
             namespace=f"{get_tracer_label(tracer)}_z{iz}",
             zeff=float(sp.attrs['zeff']),file=str(fn))
    records.append(rec)
    blocks.append(dict(
        tracer=tracer,namespace=rec['namespace'],zeff=rec['zeff'],
        data=sp.value(),precision=np.linalg.inv(d.covariance.value()),
        window=win,kin=np.asarray(kin,float),
        shotnoise=float(np.mean(sp.get(0).values('shotnoise')))
    ))
zvals=[0.0]+sorted({float(r['zeff']) for r in records})
z_to_i={round(z,8):i+1 for i,z in enumerate(zvals)}
(OUT/'desi_zeff.json').write_text(json.dumps(records,indent=2))
print('CAND_DESI_ZVALS',zvals,flush=True)

def sdmc_ini(c,root):
    h=c['H0']/100.; Ox=1.-(c['ob']+c['oc']+OR)/(h*h); As=math.exp(c['lnAs'])/1e10
    zpk=','.join(f'{z:.10g}' for z in zvals)
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
    gauge = synchronous
    modes = s
    output = mPk,dTk,vTk
    extra_metric_transfer_functions = yes
    matter_source_in_current_gauge = no
    P_k_max_h/Mpc = 2.0
    z_pk = {zpk}
    output_background_smg = 3
    write background = yes
    write thermodynamics = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    output_verbose = 0
    """)

def lcdm_ini(c,root):
    As=math.exp(c['lnAs'])/1e10
    zpk=','.join(f'{z:.10g}' for z in zvals)
    return textwrap.dedent(f"""\
    H0 = {c['H0']}
    omega_b = {c['ob']}
    omega_cdm = {c['oc']}
    N_ncdm = 0
    N_ur = 3.046
    T_cmb = 2.7255
    YHe = 0.2453
    A_s = {As:.17e}
    n_s = {c['ns']}
    tau_reio = {c['tau']}
    gauge = synchronous
    modes = s
    output = mPk,dTk,vTk
    extra_metric_transfer_functions = yes
    matter_source_in_current_gauge = no
    P_k_max_h/Mpc = 2.0
    z_pk = {zpk}
    write background = yes
    write thermodynamics = yes
    root = {root}
    format = class
    input_verbose = 0
    background_verbose = 0
    thermodynamics_verbose = 0
    perturbations_verbose = 0
    spectra_verbose = 0
    output_verbose = 0
    """)

def run_class(c,is_lcdm=False):
    prefix=c['label'].lower()
    root=str(OUT/(prefix+'_'))
    ip=OUT/(prefix+'.ini')
    ip.write_text(lcdm_ini(c,root) if is_lcdm else sdmc_ini(c,root))
    cp=subprocess.run(['./class',str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=600)
    print('CAND_DESI_CLASS',c['label'],cp.returncode,flush=True)
    if cp.returncode:
        print(cp.stdout[-3000:],flush=True)
        raise RuntimeError(f"CLASS failed {c['label']}")
    return prefix

for c in cands: c['prefix']=run_class(c)
lcdm['prefix']=run_class(lcdm,True)

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith('#') and re.search(r'1\s*:',l)][-1].lstrip('#').strip()
    ms=list(re.finditer(r'(\d+)\s*:\s*',hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def rd_from(prefix,bg):
    th=tab(OUT/f'{prefix}_00_thermodynamics.dat').sort_values('z')
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i]); cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross: raise RuntimeError(f'no drag crossing {prefix}')
    zd=float(cross[0])
    return float(np.interp(zd,bg.z,bg['comov.snd.hrz.'])),zd

def model_data(c):
    prefix=c['prefix']; H0=c['H0']; ob=c['ob']; oc=c['oc']; ns=c['ns']
    h=H0/100.; fb=ob/(ob+oc); fc=1.-fb
    bg=tab(OUT/f'{prefix}_00_background.dat').sort_values('z')
    rd,zd=rd_from(prefix,bg)
    series={}
    for z in zvals:
        ii=z_to_i[round(z,8)]
        pk=np.loadtxt(OUT/f'{prefix}_00_z{ii}_pk.dat')
        tk=tab(OUT/f'{prefix}_00_z{ii}_tk.dat')
        kh=tk['k (h/Mpc)'].to_numpy()
        dm=tk['d_m'].to_numpy()
        tb=tk['t_b'].to_numpy()
        hp=tk['h_prime'].to_numpy(); ep=tk['eta_prime'].to_numpy()
        shift=.5*(hp+6.*ep)
        a=1./(1.+z)
        H=float(np.interp(z,bg.z,bg['H [1/Mpc]']))
        Hc=a*H
        theta_b_gi=tb+shift
        theta_c_gi=shift
        vb=-theta_b_gi/Hc; vc=-theta_c_gi/Hc
        vcb=fb*vb+fc*vc
        ratio_v=np.divide(vcb,dm,out=np.zeros_like(vcb),where=np.abs(dm)>1e-300)
        kpk=pk[:,0]; Pdd=pk[:,1]
        rv=np.interp(kpk,kh,ratio_v)
        Ptt=Pdd*rv**2
        series[z]=(kpk,Pdd,Ptt)
    return dict(label=c['label'],prefix=prefix,H0=H0,h=h,ob=ob,oc=oc,fb=fb,fc=fc,ns=ns,bg=bg,rd=rd,zd=zd,series=series)

models=[model_data(c) for c in cands]+[model_data(lcdm)]
kobs=blocks[0]['kin']
for b in blocks:
    if not np.allclose(b['kin'],kobs): raise RuntimeError('DESI theory k grids differ')

starts={
 'BGS_z0':(1.11348,0.660148,-0.223088),
 'LRG_z0':(1.1351,-0.139101,-0.910397),
 'LRG_z1':(1.22446,-0.671727,-0.22624),
 'LRG_z2':(1.04564,-0.493139,-0.24021),
 'ELG_z1':(0.127936,-0.578801,-1.45356),
 'QSO_z0':(0.780983,0.353843,0.135725),
}
fid=DESI(engine='camb')
all_marg=['alpha0p','alpha2p','alpha4p','alpha6p','sn0p','sn2p','sn4p']
scales=np.array([12.5]*4+[2.]+[5.]*2)
marg=['alpha0p','alpha2p','sn0p','sn2p']
gi=np.array([all_marg.index(x) for x in marg])
prior_hess=-np.diag(scales[gi]**-2)

def prepare_basis(m):
    kin=np.geomspace(min(5e-4,kobs[0]/2),max(1.0,kobs[-1]*2),500)
    zs=np.array(zvals)
    Pdd=np.stack([np.interp(kin,*m['series'][z][:2]) for z in zs],axis=-1)
    Ptt=np.stack([np.interp(kin,m['series'][z][0],m['series'][z][2]) for z in zs],axis=-1)
    pki=PowerSpectrumInterpolator1D(kin,Pdd)
    cosmo=Cosmology(n_s=m['ns'],Omega_b=m['ob']/m['h']**2,
                    Omega_cdm=m['oc']/m['h']**2,Omega_ncdm=0.,H0=m['H0'])
    cosmo.rs_drag=m['rd']*m['h']
    filt=PowerSpectrumBAOFilter(pki,engine='peakaverage',cosmo=cosmo,cosmo_fid=fid)
    filt(pki,cosmo=cosmo)
    Pnw=filt.smooth_pk_interpolator()(kin)
    sig8=np.asarray(pki.sigma8())
    ptti=PowerSpectrumInterpolator1D(kin,Ptt)
    fs8=np.asarray(ptti.sigma8())
    pt=REPT(kin,Pdd[:,0],pnw=Pnw[:,0],kmin=kobs[0],kmax=kobs[-1],
            nk=200,rbao=110,sbao=None,beyond_gauss=True,one_loop=True,
            shear=True,cutoff=20,jn=5,N=4000,threads=2,
            extrap_min=-4,extrap_max=3,import_wisdom=False)
    extk=np.append(pt.kv,1.)
    def loginterp(arr):
        return 10**interp1d(np.log10(kin),np.log10(np.maximum(arr,1e-300)),
                           kind='cubic',fill_value='extrapolate',axis=0,
                           assume_sorted=True)(np.log10(extk))
    pcb=loginterp(Pdd); pnw=loginterp(Pnw); ptt=loginterp(Ptt)
    tables={}
    for b in blocks:
        z=b['zeff']; iz=zvals.index(z)
        Hkms=float(np.interp(z,m['bg'].z,m['bg']['H [1/Mpc]']))*C
        DA=float(np.interp(z,m['bg'].z,m['bg']['comov. dist.']))/(1.+z)
        qpar=float(fid.efunc(z)/(Hkms/(100.*m['h'])))
        qper=float(DA*m['h']/fid.angular_diameter_distance(z))
        Dz=float(np.sqrt(pcb[-1,iz]/pcb[-1,0]))
        fk=np.sqrt(np.maximum(ptt[:-1,iz]/pcb[:-1,iz],0.))
        pks=pt.compute_redshift_space_power_multipoles_tables(
            fk,apar=qpar,aperp=qperp if False else qper,ngauss=4,
            pcb=pcb[:-1,iz],pcb_nw=pnw[:-1,iz],Dz=Dz)[1:]
        basis=np.stack([
            interp1d(pt.kv,pks[j],kind='cubic',fill_value='extrapolate',
                     axis=0,assume_sorted=True)(kobs) for j in range(3)
        ],axis=0)
        tables[b['namespace']]=dict(basis=basis,sigma8=float(sig8[iz]),
                                     fsigma8=float(fs8[iz]),qpar=qpar,qper=qper,Dz=Dz)
    return tables

def get_poles(info,pars,settings,sn,return_gradient=False):
    pktable=info['basis']; sigma8=info['sigma8']; f=info['fsigma8']/sigma8
    b1p,b2p,bsp=pars; b3p=0.
    b1L=b1p/sigma8-1.; b2L=b2p/sigma8**2
    bsL=bsp/sigma8**2; b3L=b3p/sigma8**3
    b1=1.+b1L; b2=8./21.*b1L+b2L; bs=bsL-(2./7.)*b1L; b3=3*b3L+b1L
    gradient=np.zeros((7,7))
    gradient[0,0]=b1*b1
    gradient[1,0]=gradient[1,1]=f*b1
    gradient[2,1]=f*f; gradient[2,2]=f*b1; gradient[3,2]=f*f
    fsat=settings['fsat']; sigv=settings['sigv']
    for jj,pole in enumerate([0,2,4]):
        gradient[jj-3,jj-3]=sn*(fsat if pole>0 else 1.)*sigv**pole
    base=[1,b1,b1*b1,b2,b1*b2,b2*b2,bs,b1*bs,b2*bs,bs*bs,b3,b1*b3]
    nuis=[0.]*7
    mono=np.concatenate([np.asarray(base),gradient.dot(np.asarray(nuis))])
    poles=np.sum(pktable*mono,axis=-1)
    if return_gradient: return poles,pktable[...,-7:].dot(gradient)
    return poles

def profile_model(m):
    label=m['label']; tables=prepare_basis(m)
    rows=[]; total=0.
    for b in blocks:
        settings=get_physical_stochastic_settings(tracer=b['tracer'].upper()[:3])
        def block_logp(x):
            poles,grad=get_poles(tables[b['namespace']],x,settings,b['shotnoise'],True)
            theory=b['window'].dot(poles.ravel())
            grad=grad[...,gi].reshape(-1,len(gi)); grad=b['window'].dot(grad)
            diff=theory-b['data']
            pgrad=b['precision'].dot(grad)
            postgrad=-pgrad.T.dot(diff)
            likeh=-grad.T.dot(pgrad); posth=prior_hess+likeh
            dx=-np.linalg.solve(posth,postgrad)
            logl=-.5*diff.T.dot(b['precision']).dot(diff)
            logl+=.5*dx.dot(likeh).dot(dx)+postgrad.dot(dx)
            logprior=.5*dx.dot(prior_hess).dot(dx)
            logl+=-.5*np.linalg.slogdet(-posth)[1]
            logbias=-.5*(x[1]/5.)**2-.5*(x[2]/5.)**2
            return float(logl+logprior+logbias)
        seeds=[
          np.array(starts[b['namespace']],float),np.array([1.,0.,0.]),
          np.array([.55,4.,-4.]),np.array([1.8,-4.,4.]),
          np.array([2.6,8.,-8.]),np.array([1.4,-10.,-10.])
        ]
        trials=[]
        for iseed,seed in enumerate(seeds):
            opt=minimize(lambda x:-2*block_logp(x),seed,method='L-BFGS-B',
                         bounds=[(0.,3.),(-20.,20.),(-20.,20.)],
                         options={'maxiter':800,'ftol':1e-11,'gtol':3e-7,'maxls':50})
            trials.append(opt)
            print('CAND_DESI_TRIAL',label,b['namespace'],iseed,float(opt.fun),bool(opt.success),flush=True)
        finite=[o for o in trials if np.isfinite(o.fun)]
        opt=min(finite,key=lambda o:o.fun)
        vals=np.array([float(o.fun) for o in finite])
        chi=float(opt.fun); total+=chi; info=tables[b['namespace']]
        row=dict(model=label,tracer=b['tracer'],namespace=b['namespace'],zeff=b['zeff'],
                 chi2_profile=chi,success=bool(opt.success),
                 start_spread_chi2=float(vals.max()-vals.min()) if len(vals)>1 else 0.,
                 n_success=sum(bool(o.success) for o in trials),
                 b1p=float(opt.x[0]),b2p=float(opt.x[1]),bsp=float(opt.x[2]),
                 sigma8=info['sigma8'],fsigma8=info['fsigma8'],
                 qpar=info['qpar'],qper=info['qper'],Dz=info['Dz'])
        rows.append(row); print('CAND_DESI_BLOCK',json.dumps(row,sort_keys=True),flush=True)
    print('CAND_DESI_TOTAL',label,total,flush=True)
    return rows,total

allrows=[]; sums=[]
for m in models:
    rows,total=profile_model(m)
    allrows.extend(rows); sums.append(dict(model=m['label'],chi2_profile=total))
lcdmchi=[x['chi2_profile'] for x in sums if x['model']=='LCDM'][0]
for s in sums:
    s['delta_vs_lcdm']=s['chi2_profile']-lcdmchi
    print('CAND_DESI_DELTA',s['model'],s['delta_vs_lcdm'],flush=True)

pd.DataFrame(allrows).to_csv(OUT/'candidate_desi_blocks.csv',index=False)
pd.DataFrame(sums).to_csv(OUT/'candidate_desi_summary.csv',index=False)
