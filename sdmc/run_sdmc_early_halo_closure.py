#!/usr/bin/env python3
"""Reproducible *diagnostic* SDMC early halo / SMBH checks (2026-10-12).

The 3D test uses first-order Zel'dovich/Lagrangian dust, NOT the full SDMC
Einstein-Horndeski+spatial-map nonlinear equations or an N-body solver.
It does not derive g3, alpha_K, or an exact JWST catalog likelihood.
Input numbers trace to the public late300_candidate.json and manuscript.

Requires numpy and scipy. Outputs JSON and CSV in --output directory.
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.constants import c, G, m_p, pi, physical_constants
sigma_T = physical_constants["Thomson cross section"][0]

PARAM = {
    'H0': 69.71482083084993, 'omega_b': 0.022083219194622913,
    'omega_cdm': 0.12299536722293603,
    'Omega_x0': 0.7014079815036496, 'lambda_e': 18.40625,
    'z_t': 16.19317622240633, 'A_late': 0.024822032167576252,
    'B_late': 0.01264704344701022, 'A_F': 0.02048146490100771,
}
C_KMS = c/1000
MPC_M = 3.085677581491367e22
SEC_PER_GYR = 365.25*24*3600*1e9
RHOCRIT_MSUN_MPC3_COEF = 2.77536627e11
# Standard photons + three neutrino species approximately; fixed baseline used
# by Late300's structural derivation script.
OMEGAR_PHYS = 4.17998772e-5


def E_of_z(z):
    z = float(z)
    h = PARAM['H0']/100
    om = (PARAM['omega_b']+PARAM['omega_cdm'])/h**2
    orad = OMEGAR_PHYS/h**2
    ox = 1-om-orad
    a = 1/(1+z)
    x_r = orad*(1+z)**4/(om*(1+z)**3+orad*(1+z)**4)
    # window half-width in ln a = 0.5
    W = .5*(1-np.tanh((np.log(a)+np.log1p(PARAM['z_t']))/.5))
    f = W*(3+x_r)/PARAM['lambda_e']**2
    deltaH = PARAM['A_late']*z*np.exp(-z/.25) - PARAM['B_late']*z*z*np.exp(-z/1.5)
    return float(np.sqrt((om*(1+z)**3+orad*(1+z)**4+ox)/(1-f))*(1+deltaH))


def age_gyr(z):
    # H0^-1 * int_(ln a_min)^(ln a(z)) dln(a)/E. Radiation-era
    # remainder below z~1e8 is negligible (<0.01 yr).
    top=-np.log1p(z)
    def integrand(ln_a):
        return 1.0/E_of_z(np.expm1(-ln_a))
    v,err=quad(integrand, -20, top, epsabs=1e-11, epsrel=1e-10, limit=400)
    H0_s=PARAM['H0']*1000/MPC_M
    return v/H0_s/SEC_PER_GYR


def masses_and_screening():
    h=PARAM['H0']/100
    rho_m = RHOCRIT_MSUN_MPC3_COEF*(PARAM['omega_b']+PARAM['omega_cdm'])
    rows=[]
    kc=3.0 # Mpc^-1, illustrative high-k cutoff from manuscript, not inferred
    for m in [1e8,1e9,1e10,1e11,1e12,3e12,5e12]:
        r=(3*m/(4*np.pi*rho_m))**(1/3)
        k=1/r
        filt=1/(1+(k/kc)**4)
        kc_z20=kc*(1+5.5)/(1+20)  # fixed physical screening wavenumber
        filt20=1/(1+(k/kc_z20)**4)
        am=1+.2123*filt
        am20=1+.2123*filt20
        oldz=20
        newz=(1+oldz)*am-1
        newz20=(1+oldz)*am20-1
        delta_t_Myr=(age_gyr(oldz)-age_gyr(newz))*1000
        delta_t_fixedphys_Myr=(age_gyr(oldz)-age_gyr(newz20))*1000
        rows.append({'halo_mass_Msun':m,'lagrangian_R_Mpc':r,'k_proxy_Mpc_inv':k,
                     'illustrative_filter_kc_3':filt,'effective_AM_if_unfiltered_1p2123':am,
                     'delta_t_Myr_z20_fixed_comoving_kc':delta_t_Myr,
                     'filter_z20_if_fixed_physical_screening':filt20,
                     'AM_z20_if_fixed_physical_screening':am20,
                     'delta_t_Myr_z20_fixed_physical_screening':delta_t_fixedphys_Myr})
    # alpha_K normalized to mass term at a chosen reference k.
    a=1/(1+5.5)
    alpha_ratio=a**4/kc**4
    # A 90% retained halo response at 1e11 Msun and <=10% retained at
    # 1e9 Msun imply a finite, but not uniquely selected, screening band.
    k_h=rows[3]['k_proxy_Mpc_inv']
    k_d=rows[1]['k_proxy_Mpc_inv']
    kc_min=k_h/(1/.9-1)**.25
    kc_max=k_d/(1/.1-1)**.25
    return {'rho_m0_Msun_Mpc3':rho_m,'screening_rows':rows,
            'alpha_K_over_mell2_in_Mpc4_for_kc3_z5p5':alpha_ratio,
            'cutoff_band_kc_Mpc_inverse_90pct_1e11_and_10pct_1e9':[kc_min,kc_max],
            'k_c_z20_fixed_physical_from_k_c3_at_z5p5':kc*(1+5.5)/(1+20),
            'required_early_physical_kc_increase_to_keep_90pct_response_1e11':kc_min/(kc*(1+5.5)/(1+20)),
            'required_early_m2_over_alpha_increase_to_keep_90pct_response_1e11':(kc_min/(kc*(1+5.5)/(1+20)))**4,

            'caveat':'k_c and spatial-mode mass not derived; k~1/R is a diagnostic conversion.'}


def jacobian_3d(N=40, seed=20261012):
    rng=np.random.default_rng(seed)
    noise=rng.normal(size=(N,N,N))
    kx=np.fft.fftfreq(N)*N
    ky=np.fft.fftfreq(N)*N
    kz=np.fft.rfftfreq(N)*N
    k1,k2,k3=np.meshgrid(kx,ky,kz,indexing='ij')
    ksq=k1*k1+k2*k2+k3*k3
    kval=np.sqrt(ksq)
    # Smooth, broad-band Gaussian density field; illustrative, not Late300 P(k).
    # The absolute amplitude is fixed below by a target caustic fraction.
    shape=np.where(kval>0,(kval/2)**.5/(1+(kval/8)**4)**.5,0)
    delta_k=np.fft.rfftn(noise)*shape
    delta_k[0,0,0]=0
    scale=max(float(np.std(np.fft.irfftn(delta_k,s=(N,N,N),axes=(0,1,2)))),1e-30)
    delta_k/=scale
    comps=[k1,k2,k3]
    T=np.empty((N,N,N,3,3),float)
    den=np.where(ksq>0,ksq,1)
    for i in range(3):
        for j in range(3):
            T[...,i,j]=np.fft.irfftn(delta_k*comps[i]*comps[j]/den,s=(N,N,N),axes=(0,1,2))
    eigen=np.linalg.eigvalsh(T)
    maxi=float(eigen[...,-1].max())
    D=.85/maxi  # pre-shell-crossing; max(D lambda)=.85
    dlam=D*eigen
    J=np.prod(1-dlam,axis=-1)
    density=1/J
    trace=dlam.sum(axis=-1)
    # Exactly J=1-I1+I2-I3 in the scaled deformation tensor.
    invariant=1-trace+(dlam[...,0]*dlam[...,1]+dlam[...,0]*dlam[...,2]
                      +dlam[...,1]*dlam[...,2])-np.prod(dlam,axis=-1)
    iso=(1-trace/3)**3
    # Zel'dovich directional expansion rates at matter-era f_D=1
    haxes=1-dlam/(1-dlam)
    hloc=haxes.mean(axis=-1)
    shear_sq_by_Hsq=.5*np.sum((haxes-hloc[...,None])**2,axis=-1)
    delta=density-1
    # Proxies: exact anisotropic Hamiltonian includes shear, but the Zel'dovich
    # field is not itself a solution of the full SDMC constraints.
    curvature_proxy=1+delta-hloc*hloc+shear_sq_by_Hsq/3
    u = -PARAM['A_F']*np.maximum(hloc,0)*np.maximum(curvature_proxy,0)**3
    # Avoid exponentiating highly nonlinear invalid patch responses.
    active=(delta>1)&(hloc>0)
    field_summary={
        'grid_cells':int(N**3),'seed':seed,'D_selected':float(D),
        'max_D_eigenvalue':float(np.max(dlam)),
        'min_J':float(J.min()),'mean_J':float(J.mean()),
        'min_positive_J':float(np.min(J[J>0])),
        'max_det_invariant_error':float(np.max(np.abs(J-invariant))),
        'mass_conservation_cellwise_max_error':float(np.max(np.abs(J*density-1))),
        'max_AMGM_violation':float(np.max(J-iso)),
        'mean_aniso_to_iso_J_ratio':float(np.mean(J/iso)),
        'fraction_cells_delta_above_1':float(np.mean(delta>1)),
        'fraction_cells_delta_above_3':float(np.mean(delta>3)),
        'fraction_cells_turnaround_hloc_nonpositive':float(np.mean(hloc<=0)),
        'active_fraction':float(np.mean(active)),
        'max_negative_u':float(np.min(u)),
        'active_u_median':float(np.median(u[active])) if np.any(active) else None,
        'active_u_p95_amplitude':float(np.percentile(-u[active],95)) if np.any(active) else None,
        'fraction_active_u_below_minus_0p1':float(np.mean(u[active]<-.1)) if np.any(active) else None,
        'physical_assumption':'Pre-shell-crossing 3D Zel\'dovich and anisotropic local-curvature proxy; not full SDMC Einstein/Horndeski N-body.'
    }
    # Explicit anisotropy benchmark at fixed trace.
    Dtest=.5
    eig_iso=np.array([.5,.5,.5]); eig_an=np.array([1.2,.3,0])
    ji=np.prod(1-Dtest*eig_iso); ja=np.prod(1-Dtest*eig_an)
    field_summary['fixed_trace_case']={'D':Dtest,'isotropic_eigenvalues':eig_iso.tolist(),
                                       'anisotropic_eigenvalues':eig_an.tolist(),
                                       'J_iso':float(ji),'J_aniso':float(ja),
                                       'rho_iso_over_background':float(1/ji),
                                       'rho_aniso_over_background':float(1/ja),
                                       'first_caustic_D_iso':float(1/np.max(eig_iso)),
                                       'first_caustic_D_aniso':float(1/np.max(eig_an))}
    return field_summary


def smbhs():
    tE = sigma_T*c/(4*np.pi*G*m_p)/SEC_PER_GYR
    # BH accretion efficiency epsilon, Eddington ratio lambda, active fraction duty.
    rows=[]
    for zobs in [6,7,8,10,12]:
        age_obs=age_gyr(zobs)
        for zseed in [20,23.26,30,34.82]:
            age_seed=age_gyr(zseed)
            dt=age_obs-age_seed
            for eps,edd,duty in [(0.1,1.,1.),(0.1,1.,0.5),(.15,1.,1.),(0.1,2.,1.)]:
                ef=tE*eps/(1-eps)/(edd*duty)
                growth=np.exp(np.clip(dt/ef,-650,650)) if dt>=0 else np.nan
                rows.append({'z_observed':zobs,'z_seed_formation':zseed,'age_observed_Gyr':age_obs,
                             'age_seed_Gyr':age_seed,'available_Gyr':dt,
                             'radiative_efficiency':eps,'Eddington_ratio':edd,'duty_cycle':duty,
                             'e_fold_time_Myr':ef*1000,'mass_growth_factor':growth,
                             'min_seed_for_1e9_Msun':1e9/growth if growth>0 else np.nan})
    bb=[]
    fb=PARAM['omega_b']/(PARAM['omega_b']+PARAM['omega_cdm'])
    for mstar in [1.096e11,1.514e11,2.344e11]:
        for eps in [.20,.315,.63]:
            bb.append({'Mstar_Msun':mstar,'epsilon':eps,'cosmic_baryon_fraction':fb,
                       'minimum_host_Mhalo_Msun':mstar/(eps*fb)})
    return {'tE_Gyr':tE,'ages_Gyr':{str(z):age_gyr(z) for z in [0,5.5,6,7,8,10,12,20,23.26,30,34.82]},
            'growth_rows':rows,'stellar_host_bounds':bb,
            'warning':'Illustrative continuous exponential accretion, not prediction of a seed mass, fuel supply, mergers, or gas fragmentation.'}


def data_audit():
    return {
     'repository_jwst_script':'sdmc/run_g022a034_jwst_likelihood_equivalent.py',
     'published_author_posterior':'Comini, Vagnozzi and Loeb, arXiv:2604.13866 (2026)',
     'raw_fresco_data':'MAST JWST program 1895; DOI:10.17909/gdyc-7g80',
     'advanced_study_dataset':'Available from corresponding author on reasonable request per Xiao et al., Nature 635 (2024)',
     'exact_likelihood_reproduced':False,
     'reason':'Individual catalog completeness/selection and joint likelihood/covariance plus model-derived nonlinear HMF are not supplied by accepted Late300 outputs or repository.',
     'late300_published_conditional':{'g3':.02048146490100771,'A_eff':1.20743461646488,'epsilon95_equiv':.31905230682809604},
     'source_sanity':'One-sided likelihood posterior derived from published constraints is not a sample-by-sample catalog posterior.'
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='sdmc_early_halo_results')
    parser.add_argument('--grid',type=int,default=40)
    args=parser.parse_args()
    output=Path(args.output);output.mkdir(parents=True,exist_ok=True)
    screening=masses_and_screening()
    geometry=jacobian_3d(args.grid)
    bh=smbhs()
    data=data_audit()
    out={'program':'SDMC early halo/SMBH complementary derivation diagnostics',
         'input_params':PARAM,'geometry':geometry,'screening':screening,
         'black_holes':bh,'observational_audit':data}
    (output/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False))
    for name,rows in [('screening_mass_table.csv',screening['screening_rows']),
                      ('smbh_growth_table.csv',bh['growth_rows']),
                      ('stellar_host_mass_table.csv',bh['stellar_host_bounds'])]:
        with (output/name).open('w',newline='') as f:
            dw=csv.DictWriter(f,fieldnames=list(rows[0])); dw.writeheader();dw.writerows(rows)
    print(json.dumps({'geometry':geometry, 'alpha_ratio':screening['alpha_K_over_mell2_in_Mpc4_for_kc3_z5p5'],
      'ages_Gyr':bh['ages_Gyr'], 'published_JWST':data['late300_published_conditional']},indent=2))

if __name__=='__main__':main()