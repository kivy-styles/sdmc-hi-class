#!/usr/bin/env python3
"""Add the SDMC split-sector Boltzmann completion.

This diagnostic deliberately separates the two pieces that the consolidated
manuscript treats as conceptually distinct:

  1. the early lambda_e=20 tracker is represented as an ordinary conserved
     fluid with the exact frozen homogeneous density history and rest-frame
     sound speed c_phi^2=0.003;
  2. the late structural release is represented by a separate Horndeski
     degree of freedom with alpha_B=-2 alpha_M.

The motivation is to test whether the very large growth found when one
hi_class scalar is forced to carry both sectors is an inherited dynamical-mode
artifact.  The expansion history is NOT changed: sdmc_full still prescribes
exactly the frozen background.  The tracker fluid merely moves the early
tracker energy/perturbations out of rho_smg and into the standard fluid sector.

Two late gravity closures are registered:

  sdmc_v2_late_noslip  - the recovered NKp-v2 B0 F(a) trajectory;
  sdmc_kp_late_noslip  - the frozen analytic Kp release
                         (A_M=0.02586, z_c=4.38, w=0.4467).

For both, the scalar is a spectator before the late release.  A tiny positive
D floor regularises the decoupled scalar, while the active late kinetic term is
chosen from the *native hi_class* cs2 numerator so that the requested late
sound speed is obtained wherever the release is active.  Native hi_class
stability checks remain enabled.

This is a split-sector completion test, not a claim that the auxiliary tracker
fluid is a unique microscopic Lagrangian.
"""
from pathlib import Path


def replace_once(path: str, old: str, new: str) -> None:
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if new in text:
        print(f"{path}: split-sector patch already present")
        return
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{path}: expected one patch anchor, found {n}:\n{old}")
    p.write_text(text.replace(old, new), encoding="utf-8")
    print(f"{path}: patched")


# ---------------------------------------------------------------------------
# 1) Standard-fluid carrier for the frozen early tracker.
# ---------------------------------------------------------------------------
replace_once(
    "include/background.h",
    "enum equation_of_state {CLP,EDE};",
    "enum equation_of_state {CLP,EDE,SDMC_TRACKER};",
)

# Accept fluid_equation_of_state = SDMC_TRACKER.
replace_once(
    "source/input.c",
    '''      else if ((strstr(string1,"EDE") != NULL) || (strstr(string1,"ede") != NULL)) {
        pba->fluid_equation_of_state = EDE;
      }
      else {
        class_stop(errmsg,"incomprehensible input '%s' for the field 'fluid_equation_of_state'",string1);
      }''',
    '''      else if ((strstr(string1,"EDE") != NULL) || (strstr(string1,"ede") != NULL)) {
        pba->fluid_equation_of_state = EDE;
      }
      else if ((strstr(string1,"SDMC_TRACKER") != NULL) || (strstr(string1,"sdmc_tracker") != NULL)) {
        pba->fluid_equation_of_state = SDMC_TRACKER;
      }
      else {
        class_stop(errmsg,"incomprehensible input '%s' for the field 'fluid_equation_of_state'",string1);
      }'''
)

replace_once(
    "source/input.c",
    '''    if (pba->fluid_equation_of_state == EDE) {
      /** 8.a.2.3) Equation of state of the fluid in 'EDE' case */
      /* Read */
      class_read_double("w0_fld",pba->w0_fld);
      class_read_double("Omega_EDE",pba->Omega_EDE);
      class_read_double("cs2_fld",pba->cs2_fld);
    }''',
    '''    if (pba->fluid_equation_of_state == EDE) {
      /** 8.a.2.3) Equation of state of the fluid in 'EDE' case */
      /* Read */
      class_read_double("w0_fld",pba->w0_fld);
      class_read_double("Omega_EDE",pba->Omega_EDE);
      class_read_double("cs2_fld",pba->cs2_fld);
    }
    if (pba->fluid_equation_of_state == SDMC_TRACKER) {
      /** SDMC split-sector tracker: only the rest-frame sound speed is free
          here; the homogeneous density and w(a) are fixed by sdmc_full. */
      class_read_double("cs2_fld",pba->cs2_fld);
    }'''
)

# Insert exact tracker-fluid density/w helper in background.c.
p = Path("source/background.c")
text = p.read_text(encoding="utf-8")
helper_marker = "static void sdmc_tracker_fluid_target"
if helper_marker not in text:
    anchor = '#include "hi_class.h"\n'
    if text.count(anchor) != 1:
        raise RuntimeError("source/background.c: include anchor not unique")
    helper = r'''

/* -------------------------------------------------------------------------
 * SDMC split-sector early tracker fluid.
 *
 * sdmc_full prescribes
 *   H^2/H0^2 = rho_late / (1-f_trk),
 * with f_trk=W*3(1+w_b)/lambda_e^2.  Therefore the separately conserved
 * carrier density is rho_trk = f_trk/(1-f_trk) rho_late.  Its equation of
 * state follows exactly from conservation, w=-1-(1/3)d ln rho_trk/dN.
 * ------------------------------------------------------------------------- */
static void sdmc_tracker_fluid_target(struct background *pba,
                                      double a,
                                      double *rho,
                                      double *w) {
  double z = 1./a-1.;
  double Om0 = pba->Omega0_b + pba->Omega0_cdm;
  double Or0 = pba->Omega0_g + pba->Omega0_ur;
  double Ox = pba->parameters_smg[0];
  double lambda_e = pba->parameters_smg[1];
  double zt = pba->parameters_smg[2];
  double dNt = pba->parameters_smg[3];
  double A = pba->parameters_smg[4];
  double tauA = pba->parameters_smg[5];
  double B = pba->parameters_smg[6];
  double tauB = pba->parameters_smg[7];

  double rm = Om0*pow(a,-3.);
  double rr = Or0*pow(a,-4.);
  double base = rm+rr+Ox;
  double base_N = -3.*rm-4.*rr;

  double eA = exp(-z/tauA);
  double eB = exp(-z/tauB);
  double delta = A*z*eA-B*z*z*eB;
  double ddelta_dz = A*eA*(1.-z/tauA)-B*eB*(2.*z-z*z/tauB);
  double delta_N = -(1.+z)*ddelta_dz;
  double g = 1.+delta;
  double rho_late_dimless = base*g*g;
  double rho_late_N = base_N*g*g+2.*base*g*delta_N;

  double xr = rr/(rm+rr);
  double f0 = (3.+xr)/(lambda_e*lambda_e);
  double f0_N = -xr*(1.-xr)/(lambda_e*lambda_e);
  double N = log(a);
  double Nt = -log(1.+zt);
  double u = (N-Nt)/dNt;
  double th = tanh(u);
  double W = .5*(1.-th);
  double W_N = -.5*(1.-th*th)/dNt;
  double f = W*f0;
  double f_N = W_N*f0+W*f0_N;

  /* The frozen handoff leaves f positive (though tiny) at a=1. */
  if (f < 1.e-40) f = 1.e-40;
  if (f > 1.-1.e-12) f = 1.-1.e-12;

  *rho = pba->H0*pba->H0 * f/(1.-f)*rho_late_dimless;
  {
    double lnrho_N = f_N/(f*(1.-f)) + rho_late_N/rho_late_dimless;
    *w = -1. - lnrho_N/3.;
  }
}

static double sdmc_tracker_fluid_w_only(struct background *pba,double a) {
  double rho,w;
  sdmc_tracker_fluid_target(pba,a,&rho,&w);
  return w;
}
'''
    text = text.replace(anchor, anchor+helper)
    p.write_text(text, encoding="utf-8")
    print("source/background.c: inserted split tracker helper")

# Add SDMC_TRACKER to the three background_w_fld switches.
replace_once(
    "source/background.c",
    '''  case CLP:
    *w_fld = pba->w0_fld + pba->wa_fld * (1. - a);
    break;
  case EDE:''',
    '''  case CLP:
    *w_fld = pba->w0_fld + pba->wa_fld * (1. - a);
    break;
  case SDMC_TRACKER: {
    double rho_dummy;
    sdmc_tracker_fluid_target(pba,a,&rho_dummy,w_fld);
    break;
  }
  case EDE:'''
)

replace_once(
    "source/background.c",
    '''  case CLP:
    *dw_over_da_fld = - pba->wa_fld;
    break;
  case EDE:''',
    '''  case CLP:
    *dw_over_da_fld = - pba->wa_fld;
    break;
  case SDMC_TRACKER: {
    const double eps = 1.e-5;
    double wp = sdmc_tracker_fluid_w_only(pba,a*exp(eps));
    double wm = sdmc_tracker_fluid_w_only(pba,a*exp(-eps));
    *dw_over_da_fld = (wp-wm)/(2.*eps*a);
    break;
  }
  case EDE:'''
)

replace_once(
    "source/background.c",
    '''  case CLP:
    *integral_fld = 3.*((1.+pba->w0_fld+pba->wa_fld)*log(1./a) + pba->wa_fld*(a-1.));
    break;
  case EDE:''',
    '''  case CLP:
    *integral_fld = 3.*((1.+pba->w0_fld+pba->wa_fld)*log(1./a) + pba->wa_fld*(a-1.));
    break;
  case SDMC_TRACKER: {
    double rho_a,w_a,rho_0,w_0;
    sdmc_tracker_fluid_target(pba,a,&rho_a,&w_a);
    sdmc_tracker_fluid_target(pba,1.,&rho_0,&w_0);
    *integral_fld = log(rho_a/rho_0);
    break;
  }
  case EDE:'''
)

# ---------------------------------------------------------------------------
# 2) Late-only No-Slip scalar.  Runtime order: this patch is applied after
#    apply_sdmc_v2_noslip_patch.py and apply_sdmc_v2_smooth_noslip_patch.py.
# ---------------------------------------------------------------------------
replace_once(
    "include/background.h",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip,\n",
    "    constant_alphas, sdmc_early_cs, sdmc_v2_noslip, sdmc_v2_late_noslip, sdmc_kp_late_noslip,\n",
)

# Parse the two new models immediately before the stock EFT model parser.
replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''  if (strcmp(string1,"eft_alphas_power_law") == 0) {''',
    '''  if (strcmp(string1,"sdmc_v2_late_noslip") == 0) {
     pba->gravity_model_smg = sdmc_v2_late_noslip;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 5;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0. ||
                pba->parameters_2_smg[2] <= 0. || pba->parameters_2_smg[3] <= 0. ||
                pba->parameters_2_smg[4] <= 0., errmsg,
                "sdmc_v2_late_noslip requires positive Ns, c_s^2, zact, width and D floor");
   }

  if (strcmp(string1,"sdmc_kp_late_noslip") == 0) {
     pba->gravity_model_smg = sdmc_kp_late_noslip;
     pba->field_evolution_smg = _FALSE_;
     pba->M2_evolution_smg = _FALSE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 2;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
     class_test(pba->parameters_2_smg[0] <= 0. || pba->parameters_2_smg[1] <= 0., errmsg,
                "sdmc_kp_late_noslip requires positive late c_s^2 and D floor");
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {'''
)

# Add both late closures before the stock EFT alpha parametrisation.
replace_once(
    "gravity_smg/gravity_models_smg.c",
    '''  else if (pba->gravity_model_smg == eft_alphas_power_law) {''',
    r'''  else if (pba->gravity_model_smg == sdmc_v2_late_noslip) {

    double Ns_target = pba->parameters_2_smg[0];
    double cs2_late = pba->parameters_2_smg[1];
    double zact = pba->parameters_2_smg[2];
    double dNact = pba->parameters_2_smg[3];
    double Dfloor = pba->parameters_2_smg[4];
    double N = log(a);
    double F, am, amp;
    sdmc_v2_b0_profile(pba,N,Ns_target,zact,dNact,&F,&am,&amp);

    double bra = -2.*am;
    double bra_N = -2.*amp;
    double H = pvecback[pba->index_bg_H];
    double rho_o = pvecback[pba->index_bg_rho_tot_wo_smg];
    double p_o = pvecback[pba->index_bg_p_tot_wo_smg];
    double rho_s = pvecback[pba->index_bg_rho_smg];
    double p_s = pvecback[pba->index_bg_p_smg];

    /* Native hi_class Horndeski cs2 numerator (alpha_H=alpha_T=0). */
    double cs2num =
        .5*(2.-bra)*(bra+2.*am)
      + 1.5*(2.-bra)*(rho_s+p_s)/(H*H)
      - 1.5*(2.-2.*F+bra*F)/F*(rho_o+p_o)/(H*H)
      + bra_N;

    double Nact = -log(1.+zact);
    double Wlate = .5*(1.+tanh((N-Nact)/dNact));
    double Dlate = cs2num/cs2_late;
    double D = (1.-Wlate)*Dfloor + Wlate*Dlate;
    if (D < Dfloor) D = Dfloor;
    double ak = D-1.5*bra*bra;

    pvecback[pba->index_bg_kineticity_smg] = ak;
    pvecback[pba->index_bg_braiding_smg] = bra;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = am;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = F-1.;
    pvecback[pba->index_bg_M2_smg] = F;
  }

  else if (pba->gravity_model_smg == sdmc_kp_late_noslip) {

    const double AM = 0.02586;
    const double zc = 4.38;
    const double width = 0.4467;
    double cs2_late = pba->parameters_2_smg[0];
    double Dfloor = pba->parameters_2_smg[1];
    double N = log(a);
    double Nc = -log(1.+zc);
    double x = (N-Nc)/(2.*width);
    double th = tanh(x);
    double sech2 = 1./(cosh(x)*cosh(x));
    double S = .5*(1.+th);
    double F = exp(AM*S);
    double am = AM/(4.*width)*sech2;
    double amp = -AM/(4.*width*width)*sech2*th;
    double bra = -2.*am;
    double bra_N = -2.*amp;

    double H = pvecback[pba->index_bg_H];
    double rho_o = pvecback[pba->index_bg_rho_tot_wo_smg];
    double p_o = pvecback[pba->index_bg_p_tot_wo_smg];
    double rho_s = pvecback[pba->index_bg_rho_smg];
    double p_s = pvecback[pba->index_bg_p_smg];
    double cs2num =
        .5*(2.-bra)*(bra+2.*am)
      + 1.5*(2.-bra)*(rho_s+p_s)/(H*H)
      - 1.5*(2.-2.*F+bra*F)/F*(rho_o+p_o)/(H*H)
      + bra_N;

    double Dlate = cs2num/cs2_late;
    double D = (1.-S)*Dfloor + S*Dlate;
    if (D < Dfloor) D = Dfloor;
    double ak = D-1.5*bra*bra;

    pvecback[pba->index_bg_kineticity_smg] = ak;
    pvecback[pba->index_bg_braiding_smg] = bra;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_running_smg] = am;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_delta_M2_smg] = F-1.;
    pvecback[pba->index_bg_M2_smg] = F;
  }

  else if (pba->gravity_model_smg == eft_alphas_power_law) {'''
)

print("SDMC split-sector patch complete.")
