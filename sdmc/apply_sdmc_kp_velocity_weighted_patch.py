#!/usr/bin/env python3
"""Patch CLASS with the Part-VI velocity-weighted reciprocal Kp common acceleration.

This implements Eqs. (134)--(137) of the frozen SDMC Part VI construction,
rather than the earlier diagnostic that modified the photon pressure term
alone.  In the near-tight-coupling terminal era the manuscript target is

  d/deta [(1+R_Hu)/C Theta0'] + C k^2/3 Theta0 = S_K,

and a phenomenological Euler-system completion is

  F_K = (C'/C) theta + (C^2-1) k^2 Theta0/(1+R_Hu),

with the *same* acceleration F_K added to photon and baryon Euler equations.
The common force cancels from the leading photon-baryon slip equation.

R_Hu=3 rho_b/(4 rho_gamma) is the Hu-Sugiyama convention. CLASS uses its
inverse R_CLASS=4 rho_gamma/(3 rho_b), so R_Hu=1/R_CLASS.  Outside exact TCA
we use the momentum-weighted common velocity

  theta = (theta_gamma + R_Hu theta_b)/(1+R_Hu).

The terminal C(z) profile is the same fixed, non-Planck-fitted two-edge window
used in the distance-budget scan. A guard aborts if a material part of the
window overlaps CLASS tight coupling: this patch is intentionally an auditable
post-TCA implementation of the manuscript's phenomenological common-mode
completion.
"""
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument('--dz',type=float,required=True)
ap.add_argument('--zon',type=float,required=True)
ap.add_argument('--zoff',type=float,default=1059.0)
ap.add_argument('--cb',type=float,default=0.53395806)
a=ap.parse_args()
if not (0.0<a.cb<=1.0): raise SystemExit('cb must be in (0,1]')
if a.dz<=0 or a.zon<=a.zoff: raise SystemExit('require dz>0 and zon>zoff')

p=Path('source/perturbations.c')
text=p.read_text(encoding='utf-8')
marker='SDMC_KP_VELOCITY_WEIGHTED_RECIPROCAL'
if marker in text: raise SystemExit('Kp reciprocal common force already patched')

def rep(old,new):
    global text
    n=text.count(old)
    if n!=1: raise RuntimeError(f'expected one anchor, found {n}:\n{old[:300]}')
    text=text.replace(old,new,1)

# Local variables in perturbations_derivs.
old="""  double delta_g=0.,theta_g=0.,shear_g=0.;
  double delta_b,theta_b;
  double delta_idm = 0., theta_idm = 0.;
"""
new="""  double delta_g=0.,theta_g=0.,shear_g=0.;
  double delta_b,theta_b;
  /* SDMC_KP_VELOCITY_WEIGHTED_RECIPROCAL */
  double sdmc_kp_C=1.,sdmc_kp_C2=1.,sdmc_kp_dlnC_dtau=0.,sdmc_kp_W=0.;
  double sdmc_kp_F=0.,sdmc_kp_RHu=0.,sdmc_kp_theta=0.;
  double delta_idm = 0., theta_idm = 0.;
"""
rep(old,new)

# Build the fixed terminal window after a,H and CLASS R are available and after
# the photon/baryon short-cut variables have been loaded.  Insert immediately
# before the actual evolution equations.
old="""    /** - --> (e) BEGINNING OF ACTUAL SYSTEM OF EQUATIONS OF EVOLUTION */

    /* start with idm as it might be needed during (normal) tca  */
"""
new=f"""    /* SDMC reciprocal common-mode terminal force. */
    {{
      const double kp_Cb={a.cb:.15g};
      const double kp_zoff={a.zoff:.15g};
      const double kp_zon={a.zon:.15g};
      const double kp_dz={a.dz:.15g};
      double kp_z=1./a-1.;
      /* Part-VI velocity weight W_v(k), interpolated from the retained
       * density/velocity decomposition.  k is in 1/Mpc in CLASS. */
      static const double kp_kv[10]={{0.02,0.03,0.04,0.05,0.07,0.10,0.15,0.20,0.25,0.30}};
      static const double kp_wvtab[10]={{0.0678,0.2558,0.5032,0.6873,0.8559,0.9285,0.9547,0.9533,0.9350,0.8899}};
      double kp_wv=0.;
      if (k<=kp_kv[0]) {{
        kp_wv=kp_wvtab[0]*fmax(k,0.)/kp_kv[0];
      }}
      else if (k>=kp_kv[9]) {{
        kp_wv=kp_wvtab[9];
      }}
      else {{
        int ik;
        for (ik=0;ik<9;ik++) {{
          if (k>=kp_kv[ik] && k<kp_kv[ik+1]) {{
            double f=(k-kp_kv[ik])/(kp_kv[ik+1]-kp_kv[ik]);
            kp_wv=kp_wvtab[ik]+f*(kp_wvtab[ik+1]-kp_wvtab[ik]);
            break;
          }}
        }}
      }}
      /* Only the velocity-generated acoustic contribution receives the
       * full direct-domain map: C_p(k)=1-(1-C_b)W_v(k). */
      double kp_Ceff=1.-(1.-kp_Cb)*kp_wv;
      double u_on=(kp_zon-kp_z)/kp_dz;
      double u_off=(kp_z-kp_zoff)/kp_dz;
      double t_on=tanh(u_on),t_off=tanh(u_off);
      double Awin=.5*(1.+t_on),Bwin=.5*(1.+t_off);
      double dA_dz=-.5*(1.-t_on*t_on)/kp_dz;
      double dB_dz= .5*(1.-t_off*t_off)/kp_dz;
      double dW_dz;
      sdmc_kp_W=Awin*Bwin;
      if (sdmc_kp_W<1.e-12) sdmc_kp_W=0.;
      dW_dz=dA_dz*Bwin+Awin*dB_dz;
      if (sdmc_kp_W==0.) dW_dz=0.;
      sdmc_kp_C=1.-(1.-kp_Ceff)*sdmc_kp_W;
      sdmc_kp_C2=sdmc_kp_C*sdmc_kp_C;
      /* dz/deta=-H in the CLASS background convention. */
      sdmc_kp_dlnC_dtau=pvecback[pba->index_bg_H]
                         *(1.-kp_Ceff)*dW_dz/sdmc_kp_C;

      if ((sdmc_kp_W>1.e-3) &&
          (ppw->approx[ppw->index_ap_tca]==(int)tca_on)) {{
        class_stop(error_message,
                   \"SDMC reciprocal Kp window overlaps tight coupling at z=%g (W=%g); TCA completion required\",
                   kp_z,sdmc_kp_W);
      }}

      if ((sdmc_kp_W>0.) &&
          (ppw->approx[ppw->index_ap_tca]==(int)tca_off) &&
          (ppw->approx[ppw->index_ap_rsa]==(int)rsa_off)) {{
        /* CLASS R is 4 rho_gamma/(3 rho_b); manuscript R_Hu is inverse. */
        sdmc_kp_RHu=1./R;
        sdmc_kp_theta=(theta_g+sdmc_kp_RHu*theta_b)/(1.+sdmc_kp_RHu);
        sdmc_kp_F=sdmc_kp_dlnC_dtau*sdmc_kp_theta
          +(sdmc_kp_C2-1.)*k2*(delta_g/4.)/(1.+sdmc_kp_RHu);
      }}
    }}

    /** - --> (e) BEGINNING OF ACTUAL SYSTEM OF EQUATIONS OF EVOLUTION */

    /* start with idm as it might be needed during (normal) tca  */
"""
rep(old,new)

# Add the same acceleration to baryons in the exact (TCA-off) Euler equation.
old="""      dy[pv->index_pt_theta_b] =
        - a_prime_over_a*theta_b
        + metric_euler
        + k2*delta_p_b_over_rho_b
        + R*pvecthermo[pth->index_th_dkappa]*(theta_g-theta_b);
"""
new="""      dy[pv->index_pt_theta_b] =
        - a_prime_over_a*theta_b
        + metric_euler
        + k2*delta_p_b_over_rho_b
        + R*pvecthermo[pth->index_th_dkappa]*(theta_g-theta_b)
        + sdmc_kp_F;
"""
rep(old,new)

# Add the identical acceleration to photons.  We do NOT rescale the bare
# photon pressure term here; the (C^2-1)/(1+R_Hu) term in F_K is precisely the
# reciprocal common-mode correction derived in Eq. (136).
old="""        dy[pv->index_pt_theta_g] =
          k2*(delta_g/4.-s2_squared*y[pv->index_pt_shear_g])
          + metric_euler
          + photon_scattering_rate*(theta_b-theta_g);
"""
new="""        dy[pv->index_pt_theta_g] =
          k2*(delta_g/4.-s2_squared*y[pv->index_pt_shear_g])
          + metric_euler
          + photon_scattering_rate*(theta_b-theta_g)
          + sdmc_kp_F;
"""
rep(old,new)

p.write_text(text,encoding='utf-8')
print(f'Applied velocity-weighted reciprocal Kp operator: Cb={a.cb:.9g}, zoff={a.zoff:.9g}, zon={a.zon:.9g}, dz={a.dz:.9g}')
