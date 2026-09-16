#!/usr/bin/env python3
"""Smoothly connect the frozen early tracker to the NKp-v2 B0 No-Slip sector.

The first endpoint-clean implementation switched the B0 shooting ODE on with a
hard C0 boundary at z_start.  The analytic profile itself remained finite, but
hi_class obtains alpha_B' numerically from the background table.  The kink in
alpha_M' at the hard activation was therefore large enough to make the native
scalar-gradient numerator briefly negative.

Keep the original No-Slip stability convention

  N_s = -2(1+alpha_M)(h+alpha_M) - 2 alpha_M'
        - (3 Omega_m + 4 Omega_r)/F,

but replace the hard activation by a smooth target-numerator handoff.  Before
B0 activation, alpha_M=0 and F=1 are exact fixed points if the target numerator
is the natural branch value

  N_early = -2 h - (3 Omega_m + 4 Omega_r).

We blend N_early into the late positive B0 target with a tanh window centered
at z_start and use the already frozen DeltaN_handoff as its width.  Thus no new
width parameter is introduced.  Re-shooting alpha_M(0)=0 with DeltaN=0.5 gives
z_start=21.45883112 on the frozen NKp-v2 background and leaves the recovered B0
profile essentially unchanged (F0 ~ 1.02795, alpha_M peak ~ 0.0234 near z~4.39).
"""
from pathlib import Path

p = Path("gravity_smg/gravity_models_smg.c")
text = p.read_text(encoding="utf-8")

old_sig = '''static void sdmc_v2_b0_rhs(struct background *pba,
                           double N,
                           double Ns,
                           double F,
                           double am,
                           double *Fp,
                           double *amp) {'''
new_sig = '''static void sdmc_v2_b0_rhs(struct background *pba,
                           double N,
                           double Ns,
                           double zact,
                           double dNact,
                           double F,
                           double am,
                           double *Fp,
                           double *amp) {'''

old_rhs = '''  double h, Om, Or;
  sdmc_v2_background_for_b0(pba,N,&h,&Om,&Or);
  *Fp = am*F;
  *amp = -.5*(Ns + 2.*(1.+am)*(h+am) + (3.*Om+4.*Or)/F);
}'''
new_rhs = '''  double h, Om, Or;
  sdmc_v2_background_for_b0(pba,N,&h,&Om,&Or);

  /* Smoothly hand the scalar-gradient numerator from the exact early
   * alpha_M=0 fixed point to the late B0 target. */
  double Q = 3.*Om+4.*Or;
  double Ns_early = -2.*h-Q;
  double Nact = -log(1.+zact);
  double Wact = .5*(1.-tanh((N-Nact)/dNact));
  double Ns_eff = Wact*Ns_early+(1.-Wact)*Ns;

  *Fp = am*F;
  *amp = -.5*(Ns_eff + 2.*(1.+am)*(h+am) + Q/F);
}'''

old_profile_sig = '''static void sdmc_v2_b0_profile(struct background *pba,
                               double N,
                               double Ns,
                               double zstart,
                               double *F,
                               double *am,
                               double *amp) {'''
new_profile_sig = '''static void sdmc_v2_b0_profile(struct background *pba,
                               double N,
                               double Ns,
                               double zstart,
                               double dNact,
                               double *F,
                               double *am,
                               double *amp) {'''

repls = [
    (old_sig, new_sig),
    (old_rhs, new_rhs),
    (old_profile_sig, new_profile_sig),
    ("  double N0 = -log(1.+zstart);", "  /* At z=1e4 the activation window is unity to machine precision, so\n   * F=1 and alpha_M=0 are already the exact early fixed point. */\n  double N0 = -log(1.+1.e4);"),
    ("sdmc_v2_b0_rhs(pba,x,Ns,f,m,&k1f,&k1m);", "sdmc_v2_b0_rhs(pba,x,Ns,zstart,dNact,f,m,&k1f,&k1m);"),
    ("sdmc_v2_b0_rhs(pba,x+.5*dN,Ns,f+.5*dN*k1f,m+.5*dN*k1m,&k2f,&k2m);", "sdmc_v2_b0_rhs(pba,x+.5*dN,Ns,zstart,dNact,f+.5*dN*k1f,m+.5*dN*k1m,&k2f,&k2m);"),
    ("sdmc_v2_b0_rhs(pba,x+.5*dN,Ns,f+.5*dN*k2f,m+.5*dN*k2m,&k3f,&k3m);", "sdmc_v2_b0_rhs(pba,x+.5*dN,Ns,zstart,dNact,f+.5*dN*k2f,m+.5*dN*k2m,&k3f,&k3m);"),
    ("sdmc_v2_b0_rhs(pba,x+dN,Ns,f+dN*k3f,m+dN*k3m,&k4f,&k4m);", "sdmc_v2_b0_rhs(pba,x+dN,Ns,zstart,dNact,f+dN*k3f,m+dN*k3m,&k4f,&k4m);"),
    ("sdmc_v2_b0_rhs(pba,N,Ns,f,m,&dummy,amp);", "sdmc_v2_b0_rhs(pba,N,Ns,zstart,dNact,f,m,&dummy,amp);"),
    ("sdmc_v2_b0_profile(pba,N,Ns_target,zstart,&F,&am,&amp);", "sdmc_v2_b0_profile(pba,N,Ns_target,zstart,dNhand,&F,&am,&amp);"),
]

for old, new in repls:
    if new in text:
        continue
    if text.count(old) != 1:
        raise RuntimeError(f"expected one smooth No-Slip patch anchor:\n{old}")
    text = text.replace(old, new)

p.write_text(text, encoding="utf-8")
print("SDMC NKp-v2 smooth No-Slip activation patch complete.")
