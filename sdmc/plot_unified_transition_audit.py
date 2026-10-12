#!/usr/bin/env python3
"""Create plain diagnostic plots from the audited accepted target profiles."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from run_unified_transition_audit import audit

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='output/unified_transition_plots');args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    params=json.loads((Path(__file__).resolve().parent/'late300_candidate.json').read_text())['parameters']
    res=audit(params)
    z=np.expm1(np.linspace(0,np.log1p(200),1200));n=-np.log1p(z)
    AF=params['A_F'];w=params['width'];zc=params['z_c'];Df=params['D_floor'];D0=params['D0'];zt=params['z_t']
    nc=-np.log1p(zc);U=.5*(1+np.tanh((n-nc)/(2*w)))
    D=Df+D0*U;am=AF/w*U*(1-U);K=D+6*am-6*am**2
    fig,ax=plt.subplots(figsize=(9,4.5))
    ax.plot(z,U,label='Common gravitational/kinetic transition U')
    ax.axvline(zt,linestyle='--',label='Early tracker shutoff zₜ=16.19')
    ax.axvline(zc,linestyle=':',label='Gravitational centre zc=3.93')
    ax.scatter([zt,zc],[res['exact_late300_transition']['U_at_early_handoff'],.5])
    ax.set_xscale('symlog',linthresh=1)
    ax.set_xlim(0,110)
    ax.set_ylim(0,1.05)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Transition fraction U(z)')
    ax.set_title('Late300 has distinct tracker and gravitational transition epochs')
    ax.legend(fontsize=8,loc='upper right')
    fig.tight_layout();fig.savefig(out/'transition_separation.png',dpi=165);plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,4.5))
    ax.plot(z,D,label='D = Dfloor + D0 U')
    ax.plot(z,K,label='Required kinetic composite Kcov = D + 6αM − 6αM²')
    ax.plot(z,6*am**2,label='Direct braiding squared 6αM²')
    ax.axvline(zt,linestyle='--')
    ax.axvline(zc,linestyle=':')
    ax.set_xscale('symlog',linthresh=1);ax.set_xlim(0,100)
    ax.set_xlabel('Redshift z');ax.set_ylabel('Dimensionless kinetic response')
    ax.set_title('The No-Slip kinetic composite required by the accepted Late300 target')
    ax.legend(fontsize=8);fig.tight_layout()
    fig.savefig(out/'required_kinetic_composite.png',dpi=165);plt.close(fig)
    print('PLOTS',json.dumps([str(p) for p in out.glob('*.png')]))
if __name__=='__main__':main()