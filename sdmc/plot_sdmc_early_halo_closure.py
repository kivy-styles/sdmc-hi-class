#!/usr/bin/env python3
"""Plot screening and seed-mass diagnostics from run results."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent/'results'
r=json.loads((p/'results.json').read_text())
sc=r['screening']['screening_rows']
mass=np.array([x['halo_mass_Msun'] for x in sc]); f=np.array([x['illustrative_filter_kc_3'] for x in sc]);
f20=np.array([x['filter_z20_if_fixed_physical_screening'] for x in sc]);
fig,ax=plt.subplots(figsize=(9,5))
ax.plot(mass,f,'o-',label=r'Fixed comoving $k_c=3$ Mpc$^{-1}$')
ax.plot(mass,f20,'s--',label=r'Fixed physical cutoff (normalized at $z=5.5$), evaluated at $z=20$')
ax.axhline(.9,ls=':',lw=1,label='90% retained response')
ax.axhline(.1,ls=':',lw=1,label='10% retained response')
ax.set(xscale='log',xlabel=r'Lagrangian halo mass ($M_\odot$)',ylabel='Response factor',ylim=(-.03,1.06),
       title='SDMC candidate nonlinear mapper: mass–redshift selectivity')
ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(p/'screening_mass_redshift.png',dpi=180);plt.close(fig)
rows=r['black_holes']['growth_rows']
fig,ax=plt.subplots(figsize=(9,5))
for zseed in [20,23.26]:
 for duty in [1.,.5]:
    subset=sorted([x for x in rows if x['z_seed_formation']==zseed and x['radiative_efficiency']==.1
                   and x['Eddington_ratio']==1 and x['duty_cycle']==duty],key=lambda x:x['z_observed'])
    z=[x['z_observed'] for x in subset]; seed=[x['min_seed_for_1e9_Msun'] for x in subset]
    ax.plot(z,seed,marker='o',label=f'seed z={zseed}; duty={duty}')
ax.set(yscale='log',xlabel='Observation redshift',ylabel=r'Minimum seed for $10^9\,M_\odot$ black hole ($M_\odot$)',
       title='Late300 cosmic ages and idealized Eddington-limited SMBH growth')
ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(p/'smbh_seed_constraints.png',dpi=180);plt.close(fig)
print('Plots written to',p)