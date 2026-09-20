#!/usr/bin/env python3
"""
Stage-1B DESI promotion for the ordinary-cosmology six-dimensional screen.

Reads the shortlist produced by run_r032_ordinary6d_screen.py, evaluates the
official DESI DR1 P0/P2/P4 full-shape likelihood for each stable candidate
using the same independently parameterized SDMC closure, and anchors all
estimated joint totals to the corrected exact-covariant/full-Plik baseline:

  exact Planck+DESI current = +24.539116629523166
  + Pantheon+              = +29.206347490784653
  + Union3                 = +28.731528477197687
  + DES-Y5                 = +34.309916884762176

Only changes relative to the current point are imported from the screen.
This stage is for ranking/promotion, not the final scientific result.
"""

from pathlib import Path
import json, math, re, subprocess, textwrap
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize

import sys
sys.path.insert(0,'desi-kp-cosmological-likelihoods/dr1/cobaya')
import lsstypes as types
from cosmoprimo import PowerSpectrumInterpolator1D, PowerSpectrumBAOFilter, Cosmology
from cosmoprimo.fiducial import DESI
from velocileptors.EPT.ept_fullresum_varyDz_nu_fftw import REPT
from desi_fs_bao_all import list_zrange, dataset_fn, get_tracer_label, get_physical_stochastic_settings

OUT=Path('output/ordinary6d_desi')
OUT.mkdir(parents=True,exist_ok=True)
SCREEN=Path('screen_artifact/ordinary6d_shortlist.csv')
if not SCREEN.exists():
    alt=list(Path('screen_artifact').rglob('ordinary6d_shortlist.csv'))
    if len(alt)!=1:
        raise FileNotFoundError(f'ordinary6d shortlist not found, candidates={alt}')
    SCREEN=alt[0]
screen=pd.read_csv(SCREEN)

# Ensure the exact current ordinary-cosmology point is present.
current_mask=screen['id'].astype(str).eq('current')
if not current_mask.any():
    full=list(Path('screen_artifact').rglob('ordinary6d_screen.csv'))
    if len(full)!=1:
        raise RuntimeError('current point missing and full screen unavailable')
    ff=pd.read_csv(full[0])
    cur=ff[ff.id.astype(str).eq('current')]
    if cur.empty:
        raise RuntimeError('current point missing from screen outputs')
    screen=pd.concat([screen,cur],ignore_index=True)

# Remove duplicates in ordinary coordinates while retaining best screen-ranked row.
coord=['H0','omega_b','omega_cdm','n_s','ln10As','tau_reio']
screen=screen.drop_duplicates(subset=coord).copy()

# Limit the expensive promotion to a robust union: current + 6 best under each
# Planck/SN screen metric. This preserves candidates that trade Planck against
# any one alternative SN compilation.
keep={'current'}
for col in ['delta_planck','screen_planck_pp','screen_planck_union3','screen_planck_desy5']:
    if col in screen:
        keep.update(screen.nsmallest(6,col)['id'].astype(str))
cand=screen[screen.id.astype(str).isin(keep)].copy().reset_index(drop=True)
cand.to_csv(OUT/'desi_promotion_input.csv',index=False)
print('DESI_PROMOTION_NCAND',len(cand),flush=True)
print('DESI_PROMOTION_IDS',cand.id.astype(str).tolist(),flush=True)

C=299792.458
AF=0.0715; ZC=5.0; WIDTH=1.426
D0=0.34231919445927034; POWER=1.0; DFLOOR=0.045
LAMBDA_E=17.925; ZT=17.775
OMEGA_R_PHYS=4.17998772e-5

# Official DESI blocks and effective redshifts.
records=[]
blocks=[]
for tracer,iz,zrange in list_zrange:
    if 'lya' in tracer.lower():
        continue
    fn=dataset_fn('desi_likelihood',tracer,zrange,observable_name='spectrum-poles-rotated')
    d=types.read(fn)
    sp=d.observable.get('spectrum')
    theory=d.window.theory.get('spectrum')
    kin=theory.get(0).coords('k')
    blocks.append(dict(
        tracer=tracer,
        namespace=f'{get_tracer_label(tracer)}_z{iz}',
        zeff=float(sp.attrs['zeff']),
        data=sp.value(),
        precision=np.linalg.inv(d.covariance.value()),
        window=d.window.at.observable.get('spectrum').at.theory.get('spectrum').value(),
        kin=np.asarray(kin,float),
        shotnoise=float(np.mean(sp.get(0).values('shotnoise'))),
    ))
    records.append(dict(tracer=tracer,iz=iz,zrange=list(zrange),
                        namespace=f'{get_tracer_label(tracer)}_z{iz}',
                        zeff=float(sp.attrs['zeff']),file=str(fn)))
Path(OUT/'desi_zeff.json').write_text(json.dumps(records,indent=2))
zvals=[0.0]+sorted({float(r['zeff']) for r in records})
z_to_i={round(z,8):i+1 for i,z in enumerate(zvals)}
zpk=','.join(f'{z:.10g}' for z in zvals)
kobs=blocks[0]['kin']
for b in blocks:
    if not np.allclose(b['kin'],kobs):
        raise RuntimeError('DESI theory k grids differ')

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

def tab(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith('#') and re.search(r'1\s*:',l)][-1].lstrip('#').strip()
    ms=list(re.finditer(r'(\d+)\s*:\s*',hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def rd_from(prefix,bg):
    th=tab(Path(str(prefix)+'00_thermodynamics.dat')).sort_values('z')
    zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy()
    cross=[]
    for i in range(len(zz)-1):
        if (kb[i]-1.)*(kb[i+1]-1.)<=0 and kb[i]!=kb[i+1]:
            q=(1.-kb[i])/(kb[i+1]-kb[i])
            cross.append(zz[i]+q*(zz[i+1]-zz[i]))
    if not cross:
        raise RuntimeError(f'no drag crossing for {prefix}')
    zd=float(cross[0])
    return float(np.interp(zd,bg.z,bg['comov.snd.hrz.'])),zd

def make_ini(row,root):
    H0=float(row.H0); ob=float(row.omega_b); oc=float(row.omega_cdm)
    ns=float(row.n_s); tau=float(row.tau_reio); As=math.exp(float(row.ln10As))/1e10
    h=H0/100.
    Ox=1.-(ob+oc+OMEGA_R_PHYS)/(h*h)
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

def model_data(prefix,row):
    H0=float(row.H0); ob=float(row.omega_b); oc=float(row.omega_cdm); ns=float(row.n_s)
    h=H0/100.; fb=ob/(ob+oc); fc=1.-fb
    bg=tab(Path(str(prefix)+'00_background.dat')).sort_values('z')
    rd,zd=rd_from(prefix,bg)
    series={}
    for z in zvals:
        ii=z_to_i[round(z,8)]
        pk=np.loadtxt(Path(str(prefix)+f'00_z{ii}_pk.dat'))
        tk=tab(Path(str(prefix)+f'00_z{ii}_tk.dat'))
        kh=tk['k (h/Mpc)'].to_numpy()
        dm=tk['d_m'].to_numpy()
        tb=tk['t_b'].to_numpy()
        hp=tk['h_prime'].to_numpy(); ep=tk['eta_prime'].to_numpy()
        shift=.5*(hp+6.*ep)
        a=1./(1.+z)
        H=float(np.interp(z,bg.z,bg['H [1/Mpc]']))
        Hc=a*H
        vb=-(tb+shift)/Hc
        vc=-shift/Hc
        vcb=fb*vb+fc*vc
        ratio_v=np.divide(vcb,dm,out=np.zeros_like(vcb),where=np.abs(dm)>1e-300)
        kpk=pk[:,0]; Pdd=pk[:,1]
        rv=np.interp(kpk,kh,ratio_v)
        Ptt=Pdd*rv**2
        series[z]=(kpk,Pdd,Ptt)
    return dict(prefix=str(prefix),H0=H0,h=h,ob=ob,oc=oc,fb=fb,fc=fc,ns=ns,bg=bg,rd=rd,zd=zd,series=series)

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
            fk,apar=qpar,aperp=qper,ngauss=4,
            pcb=pcb[:-1,iz],pcb_nw=pnw[:-1,iz],Dz=Dz)[1:]
        basis=np.stack([
            interp1d(pt.kv,pks[j],kind='cubic',fill_value='extrapolate',
                     axis=0,assume_sorted=True)(kobs)
            for j in range(3)
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

def profile_model(label,m):
    tables=prepare_basis(m)
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
            logbias=-.5*(x[1]/5.)**2-.5*(x[2]/5.)**2
            return float(logl+logprior+logbias)
        x0=np.array(starts[b['namespace']],float)
        opt=minimize(lambda x:-2*block_logp(x),x0,method='L-BFGS-B',
                     bounds=[(0.,3.),(-20.,20.),(-20.,20.)],
                     options={'maxiter':500,'ftol':1e-10,'gtol':1e-6})
        chi=float(opt.fun); total+=chi
        info=tables[b['namespace']]
        rows.append(dict(candidate=label,tracer=b['tracer'],namespace=b['namespace'],
                         zeff=b['zeff'],chi2_profile=chi,success=bool(opt.success),
                         b1p=float(opt.x[0]),b2p=float(opt.x[1]),bsp=float(opt.x[2]),
                         sigma8=info['sigma8'],fsigma8=info['fsigma8'],
                         qpar=info['qpar'],qper=info['qper'],Dz=info['Dz']))
    return rows,total

allblocks=[]
totals=[]
for i,row in cand.iterrows():
    tag=str(row.id)
    root=OUT/(f'{i+1:03d}_{tag}_')
    ini=OUT/(f'{i+1:03d}_{tag}.ini')
    ini.write_text(make_ini(row,str(root)))
    cp=subprocess.run(['./class',str(ini)],stdout=subprocess.PIPE,
                      stderr=subprocess.STDOUT,text=True,timeout=300)
    if cp.returncode!=0:
        totals.append(dict(id=tag,status='FAIL',returncode=cp.returncode,error=cp.stdout[-1500:].replace('\n',' | ')))
        print('DESI_PROMOTION_FAIL',tag,cp.stdout[-1200:],flush=True)
        continue
    m=model_data(root,row)
    brows,total=profile_model(tag,m)
    allblocks.extend(brows)
    totals.append(dict(id=tag,status='OK',returncode=0,desi_chi2=total))
    print('DESI_PROMOTION_TOTAL',tag,total,flush=True)

pd.DataFrame(allblocks).to_csv(OUT/'desi_promotion_blocks.csv',index=False)
tdf=pd.DataFrame(totals)
merged=cand.merge(tdf,on='id',how='left')

cur=merged[merged.id.astype(str).eq('current')]
if len(cur)!=1 or cur.iloc[0].status!='OK':
    raise RuntimeError('current DESI promotion point did not evaluate cleanly')
cur=cur.iloc[0]

EXACT_PD=24.539116629523166
EXACT_PP=29.206347490784653
EXACT_U3=28.731528477197687
EXACT_DY=34.309916884762176

for c in ['delta_planck','delta_pantheonplus','delta_union3','delta_desy5']:
    if c not in merged:
        raise RuntimeError(f'missing screen column {c}')

merged['dplanck_from_current']=merged['delta_planck']-float(cur.delta_planck)
merged['ddesi_from_current']=merged['desi_chi2']-float(cur.desi_chi2)
merged['dpp_from_current']=merged['delta_pantheonplus']-float(cur.delta_pantheonplus)
merged['du3_from_current']=merged['delta_union3']-float(cur.delta_union3)
merged['ddy_from_current']=merged['delta_desy5']-float(cur.delta_desy5)

merged['estimated_exact_planck_desi']=EXACT_PD+merged.dplanck_from_current+merged.ddesi_from_current
merged['estimated_exact_plus_pantheonplus']=EXACT_PP+merged.dplanck_from_current+merged.ddesi_from_current+merged.dpp_from_current
merged['estimated_exact_plus_union3']=EXACT_U3+merged.dplanck_from_current+merged.ddesi_from_current+merged.du3_from_current
merged['estimated_exact_plus_desy5']=EXACT_DY+merged.dplanck_from_current+merged.ddesi_from_current+merged.ddy_from_current
merged.to_csv(OUT/'desi_promotion_joint_estimates.csv',index=False)

keep2={'current'}
for col in ['estimated_exact_planck_desi','estimated_exact_plus_pantheonplus',
            'estimated_exact_plus_union3','estimated_exact_plus_desy5','desi_chi2']:
    good=merged[merged.status.eq('OK')]
    for x in good.nsmallest(4,col)['id'].astype(str):
        keep2.add(x)
final=merged[merged.id.astype(str).isin(keep2)].copy()
final['promotion_best_any']=final[
    ['estimated_exact_planck_desi','estimated_exact_plus_pantheonplus',
     'estimated_exact_plus_union3','estimated_exact_plus_desy5']
].min(axis=1)
final=final.sort_values('promotion_best_any')
final.to_csv(OUT/'desi_promotion_finalists.csv',index=False)

print('DESI_PROMOTION_CURRENT_CHI2',float(cur.desi_chi2),flush=True)
print('DESI_PROMOTION_BEST_PD',merged[merged.status.eq('OK')].nsmallest(8,'estimated_exact_planck_desi')[
    coord+['id','desi_chi2','estimated_exact_planck_desi']].to_dict('records'),flush=True)
print('DESI_PROMOTION_BEST_PP',merged[merged.status.eq('OK')].nsmallest(8,'estimated_exact_plus_pantheonplus')[
    coord+['id','estimated_exact_plus_pantheonplus']].to_dict('records'),flush=True)
print('DESI_PROMOTION_BEST_U3',merged[merged.status.eq('OK')].nsmallest(8,'estimated_exact_plus_union3')[
    coord+['id','estimated_exact_plus_union3']].to_dict('records'),flush=True)
print('DESI_PROMOTION_BEST_DY',merged[merged.status.eq('OK')].nsmallest(8,'estimated_exact_plus_desy5')[
    coord+['id','estimated_exact_plus_desy5']].to_dict('records'),flush=True)
print('DESI_PROMOTION_FINALISTS',final[
    coord+['id','desi_chi2','estimated_exact_planck_desi','estimated_exact_plus_pantheonplus',
           'estimated_exact_plus_union3','estimated_exact_plus_desy5']].to_dict('records'),flush=True)
