"""Exact DESI DR1 full-shape scan for the late266 partial-Kp matter-domain phase map.

This is a candidate observable completion: it leaves the accepted late266
background, Planck spectra, SN distances and growth response unchanged and
maps only the oscillatory late-time matter-transfer phase. It must not be
confused with the retired photon-baryon terminal Kp operator.
"""
from pathlib import Path
import sys,re,json,math,csv
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

OUT=Path('output')
C=299792.458
records=json.loads((OUT/'desi_zeff.json').read_text())
zvals=[0.0]+sorted({float(r['zeff']) for r in records})
z_to_i={round(z,8):i+1 for i,z in enumerate(zvals)}

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
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy()
    cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i])
            cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross: raise RuntimeError(f'no drag crossing for {prefix}')
    zd=float(cross[0])
    return float(np.interp(zd,bg.z,bg['comov.snd.hrz.'])),zd

def model_data(prefix,H0,ob,oc,ns):
    h=H0/100.; fb=ob/(ob+oc); fc=1.-fb
    bg=tab(OUT/f'{prefix}_00_background.dat').sort_values('z')
    th=tab(OUT/f'{prefix}_00_thermodynamics.dat').sort_values('z')
    rd,zd=rd_from(prefix,bg)
    gcol='g [Mpc^-1]'
    jstar=int(np.argmax(th[gcol].to_numpy()))
    zstar=float(th.z.to_numpy()[jstar])
    rsstar=float(np.interp(zstar,bg.z,bg['comov.snd.hrz.']))
    series={}
    for z in zvals:
        ii=z_to_i[round(z,8)]
        pk=np.loadtxt(OUT/f'{prefix}_00_z{ii}_pk.dat')
        tk=tab(OUT/f'{prefix}_00_z{ii}_tk.dat')
        kh=tk['k (h/Mpc)'].to_numpy()
        k=kh*h
        dm=tk['d_m'].to_numpy()
        tb=tk['t_b'].to_numpy()
        hp=tk['h_prime'].to_numpy(); ep=tk['eta_prime'].to_numpy()
        # CLASS default matter_source_in_current_gauge=no already
        # exports d_m as the gauge-invariant matter density source.
        # For the velocity, theta_cb^GI = f_b theta_b^S + k^2 alpha,
        # with k^2 alpha=(h'+6 eta')/2 in synchronous gauge.
        shift=.5*(hp+6.*ep)
        a=1./(1.+z)
        H=float(np.interp(z,bg.z,bg['H [1/Mpc]']))
        Hc=a*H
        theta_b_gi=tb+shift
        theta_c_gi=shift
        vb=-theta_b_gi/Hc
        vc=-theta_c_gi/Hc
        vcb=fb*vb+fc*vc
        ratio_v=np.divide(vcb,dm,out=np.zeros_like(vcb),where=np.abs(dm)>1e-300)

        kpk=pk[:,0]
        Pdd=pk[:,1]
        rv=np.interp(kpk,kh,ratio_v)
        Ptt=Pdd*rv**2
        series[z]=(kpk,Pdd,Ptt)

    return dict(prefix=prefix,H0=H0,h=h,ob=ob,oc=oc,fb=fb,fc=fc,ns=ns,bg=bg,
                rd=rd,zd=zd,zstar=zstar,rsstar=rsstar,series=series)

cov=model_data('desi_cov',69.45160505326464,0.02208511574370519,0.12298509428428875,0.962555355281866)
lcdm=model_data('desi_lcdm',67.36,0.02237,0.1200,0.9649)

# Official likelihood blocks.
blocks=[]
for tracer,iz,zrange in list_zrange:
    if 'lya' in tracer.lower(): continue
    fn=dataset_fn('desi_likelihood',tracer,zrange,observable_name='spectrum-poles-rotated')
    d=types.read(fn)
    sp=d.observable.get('spectrum')
    theory=d.window.theory.get('spectrum')
    kin=theory.get(0).coords('k')
    win=d.window.at.observable.get('spectrum').at.theory.get('spectrum').value()
    blocks.append(dict(
        tracer=tracer,namespace=f'{get_tracer_label(tracer)}_z{iz}',
        zeff=float(sp.attrs['zeff']),data=sp.value(),
        precision=np.linalg.inv(d.covariance.value()),
        window=win,kin=np.asarray(kin,float),
        shotnoise=float(np.mean(sp.get(0).values('shotnoise')))
    ))
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

def prepare_basis(m,eta=0.0):
    kin=np.geomspace(min(5e-4,kobs[0]/2),max(1.0,kobs[-1]*2),500)
    zs=np.array(zvals)
    Pdd=np.stack([np.interp(kin,*m['series'][z][:2]) for z in zs],axis=-1)
    Ptt=np.stack([np.interp(kin,m['series'][z][0],m['series'][z][2]) for z in zs],axis=-1)

    # Smooth the unmodified accepted late266 spectrum first.  Kp acts only
    # on the oscillatory matter-domain component, not on the broadband shape.
    pki_raw=PowerSpectrumInterpolator1D(kin,Pdd)
    cosmo=Cosmology(n_s=m['ns'],Omega_b=m['ob']/m['h']**2,
                    Omega_cdm=m['oc']/m['h']**2,Omega_ncdm=0.,H0=m['H0'])
    cosmo.rs_drag=m['rd']*m['h']
    filt=PowerSpectrumBAOFilter(pki_raw,engine='peakaverage',cosmo=cosmo,cosmo_fid=fid)
    filt(pki_raw,cosmo=cosmo)
    Pnw=filt.smooth_pk_interpolator()(kin)

    # Partial Kp domain conversion:
    # eta=0 -> accepted late266, eta=1 -> full sharp Kp residual-tail map.
    cb=m['fb']**(1./3.)
    ceff=1.-float(eta)*(1.-cb)
    rd_eff=m['rsstar']+ceff*(m['rd']-m['rsstar'])
    phase_scale=rd_eff/m['rd']

    if abs(float(eta))>1e-14:
        wig=np.divide(Pdd,Pnw,out=np.ones_like(Pdd),where=Pnw>0.)-1.
        wig_shift=interp1d(kin,wig,kind='cubic',axis=0,bounds_error=False,
                           fill_value=0.,assume_sorted=True)(kin*phase_scale)
        Pdd_new=Pnw*(1.+wig_shift)
        # Keep the accepted growth/velocity response fixed: the Kp candidate
        # is a matter acoustic phase map, not a broadband growth deformation.
        ratio=np.divide(Ptt,Pdd,out=np.zeros_like(Ptt),where=Pdd>0.)
        Ptt_new=Pdd_new*ratio
        Pdd,Ptt=Pdd_new,Ptt_new

    pki=PowerSpectrumInterpolator1D(kin,Pdd)
    sig8=np.asarray(pki.sigma8())
    ptti=PowerSpectrumInterpolator1D(kin,Ptt)
    fs8=np.asarray(ptti.sigma8())

    # The IR-resummation acoustic scale must match the mapped wiggle phase.
    cosmo.rs_drag=rd_eff*m['h']
    extk=None; pcb=None; pnw=None; ptt=None; pt=None
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
            fk,apar=qpar,aperp=qper,ngauss=4,
            pcb=pcb[:-1,iz],pcb_nw=pnw[:-1,iz],Dz=Dz)[1:]
        basis=np.stack([
            interp1d(pt.kv,pks[j],kind='cubic',fill_value='extrapolate',
                     axis=0,assume_sorted=True)(kobs)
            for j in range(3)
        ],axis=0)
        tables[b['namespace']]=dict(basis=basis,sigma8=float(sig8[iz]),
                                     fsigma8=float(fs8[iz]),qpar=qpar,qper=qper,
                                     Dz=Dz,rd_eff=float(rd_eff),eta=float(eta),
                                     Ceff=float(ceff),phase_scale=float(phase_scale))
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
    gradient[2,1]=f*f
    gradient[2,2]=f*b1
    gradient[3,2]=f*f
    fsat=settings['fsat']; sigv=settings['sigv']
    for jj,pole in enumerate([0,2,4]):
        gradient[jj-3,jj-3]=sn*(fsat if pole>0 else 1.)*sigv**pole
    base=[1,b1,b1*b1,b2,b1*b2,b2*b2,bs,b1*bs,b2*bs,bs*bs,b3,b1*b3]
    nuis=[0.]*7
    mono=np.concatenate([np.asarray(base),gradient.dot(np.asarray(nuis))])
    poles=np.sum(pktable*mono,axis=-1)
    if return_gradient:
        return poles,pktable[...,-7:].dot(gradient)
    return poles

def profile_model(label,m,eta=0.0):
    tables=prepare_basis(m,eta=eta)
    rows=[]; total=0.
    for b in blocks:
        settings=get_physical_stochastic_settings(tracer=b['tracer'].upper()[:3])
        def block_logp(x):
            poles,grad=get_poles(tables[b['namespace']],x,settings,b['shotnoise'],True)
            theory=b['window'].dot(poles.ravel())
            grad=grad[...,gi].reshape(-1,len(gi))
            grad=b['window'].dot(grad)
            diff=theory-b['data']
            pgrad=b['precision'].dot(grad)
            postgrad=-pgrad.T.dot(diff)
            likeh=-grad.T.dot(pgrad)
            posth=prior_hess+likeh
            dx=-np.linalg.solve(posth,postgrad)
            logl=-.5*diff.T.dot(b['precision']).dot(diff)
            logl+=.5*dx.dot(likeh).dot(dx)+postgrad.dot(dx)
            logprior=.5*dx.dot(prior_hess).dot(dx)
            logl+=-.5*np.linalg.slogdet(-posth)[1]
            # Cobaya Gaussian priors on b2p,bsp:
            logbias=-.5*(x[1]/5.)**2-.5*(x[2]/5.)**2
            return float(logl+logprior+logbias)
        seeds=[
            np.array(starts[b['namespace']],float),
            np.array([1.0,0.0,0.0]),
            np.array([0.55,4.0,-4.0]),
            np.array([1.8,-4.0,4.0]),
            np.array([2.6,8.0,-8.0]),
            np.array([1.4,-10.0,-10.0]),
        ]
        trials=[]
        for iseed,seed in enumerate(seeds):
            opti=minimize(lambda x:-2*block_logp(x),seed,method='L-BFGS-B',
                          bounds=[(0.,3.),(-20.,20.),(-20.,20.)],
                          options={'maxiter':800,'ftol':1e-11,'gtol':3e-7,'maxls':50})
            trials.append(opti)
            print('DESI_COV_ROBUST_TRIAL',label,b['namespace'],iseed,float(opti.fun),
                  bool(opti.success),[float(v) for v in opti.x],flush=True)
        finite=[o for o in trials if np.isfinite(o.fun)]
        opt=min(finite,key=lambda o:o.fun)
        vals=np.array([float(o.fun) for o in finite])
        chi=float(opt.fun); total+=chi
        info=tables[b['namespace']]
        row=dict(model=label,tracer=b['tracer'],namespace=b['namespace'],
                 zeff=b['zeff'],chi2_profile=chi,success=bool(opt.success),
                 start_spread_chi2=float(vals.max()-vals.min()) if len(vals)>1 else 0.,
                 n_success=sum(bool(o.success) for o in trials),
                 b1p=float(opt.x[0]),b2p=float(opt.x[1]),bsp=float(opt.x[2]),
                 sigma8=info['sigma8'],fsigma8=info['fsigma8'],
                 qpar=info['qpar'],qper=info['qper'],Dz=info['Dz'])
        rows.append(row)
        print('DESI_FS_BLOCK',row,flush=True)
    print('DESI_FS_PROFILE_TOTAL',label,total,flush=True)
    return rows,total

# Profile the matched LCDM control once.
lcdmrows,lcdmchi=profile_model('LCDM',lcdm,eta=0.0)

# Current audited fair offsets against optimized local021.
PLANCK_FAIR=3.22609
SN_FAIR={'pantheonplus':-0.27895,'union3':-0.19040,'desy5':-0.45328}
# late266 fixed-LCDM DESI delta=-9.83689693657 while its audited fair
# local021 residual is -1.24176.
LOCAL021_DESI_SHIFT=(-1.24176)-(-9.836896936570895)

etas=[0.0,0.20,0.30,0.38,0.44,0.46,0.48,0.52,0.60,0.75,1.0]
summ=[]; allrows=list(lcdmrows)
for eta in etas:
    label=f'LATE266_KP_PHASE_eta{eta:.3f}'
    rows,chi=profile_model(label,cov,eta=eta)
    allrows.extend(rows)
    delta_fixed=chi-lcdmchi
    dfair=delta_fixed+LOCAL021_DESI_SHIFT
    pd_fair=PLANCK_FAIR+dfair
    joints={k:pd_fair+v for k,v in SN_FAIR.items()}
    worst=max(joints.values())
    cb=cov['fb']**(1./3.)
    ceff=1.-eta*(1.-cb)
    rd_eff=cov['rsstar']+ceff*(cov['rd']-cov['rsstar'])
    phase_scale=rd_eff/cov['rd']
    row=dict(eta=eta,Ceff=ceff,rd_eff=rd_eff,
             phase_scale=phase_scale,desi_chi2=chi,lcdm_chi2=lcdmchi,
             delta_desi_vs_fixed=delta_fixed,delta_desi_fair_local021=dfair,
             delta_planck_fair_local021=PLANCK_FAIR,
             delta_PD_fair_local021=pd_fair,
             joint_pantheonplus=joints['pantheonplus'],
             joint_union3=joints['union3'],joint_desy5=joints['desy5'],
             minimax_joint=worst,
             comfortable_all_three=bool(worst<=-0.5))
    summ.append(row)
    print('KP_PHASE_SCAN',json.dumps(row,sort_keys=True),flush=True)

best=min(summ,key=lambda r:r['minimax_joint'])
print('KP_PHASE_BEST',json.dumps(best,sort_keys=True),flush=True)
pd.DataFrame(summ).to_csv(OUT/'late266_kp_phase_scan_summary.csv',index=False)
pd.DataFrame(allrows).to_csv(OUT/'late266_kp_phase_scan_blocks.csv',index=False)
Path(OUT/'late266_kp_phase_best.json').write_text(json.dumps(best,indent=2,sort_keys=True)+'\n')
