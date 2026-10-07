#!/usr/bin/env python3
from pathlib import Path

p=Path("source/perturbations.c")
s=p.read_text()

def repl(old,new,label):
    global s
    if new in s:
        print(label,"already installed")
        return
    if s.count(old)!=1:
        raise RuntimeError(f"{label}: expected one anchor, found {s.count(old)}")
    s=s.replace(old,new,1)
    print(label,"installed")

old1='''      if (pba->use_ppf == _FALSE_) {
        ppw->delta_rho_fld = ppw->pvecback[pba->index_bg_rho_fld]*y[ppw->pv->index_pt_delta_fld];
        ppw->rho_plus_p_theta_fld = (1.+w_fld)*ppw->pvecback[pba->index_bg_rho_fld]*y[ppw->pv->index_pt_theta_fld];
        ca2_fld = w_fld - w_prime_fld / 3. / (1.+w_fld) / a_prime_over_a;
        /** We must gauge transform the pressure perturbation from the fluid rest-frame to the gauge we are working in */
        ppw->delta_p_fld = pba->cs2_fld * ppw->delta_rho_fld + (pba->cs2_fld-ca2_fld)*(3*a_prime_over_a*ppw->rho_plus_p_theta_fld/k/k);
      }'''
new1='''      if (pba->use_ppf == _FALSE_) {
        /* SDMC FUTURE AUX-FLD ISOLATION TEST:
           keep the accepted equations exactly unchanged through ln(a)<=2.
           Beyond that point the tiny bookkeeping fluid is excluded from
           perturbation stress-energy before its future rho_fld~0 pathology.
           This does not alter the SMG action or background. */
        if ((pba->has_smg == _TRUE_) && (a > exp(2.0))) {
          ppw->delta_rho_fld = 0.;
          ppw->rho_plus_p_theta_fld = 0.;
          ppw->delta_p_fld = 0.;
        }
        else {
          ppw->delta_rho_fld = ppw->pvecback[pba->index_bg_rho_fld]*y[ppw->pv->index_pt_delta_fld];
          ppw->rho_plus_p_theta_fld = (1.+w_fld)*ppw->pvecback[pba->index_bg_rho_fld]*y[ppw->pv->index_pt_theta_fld];
          ca2_fld = w_fld - w_prime_fld / 3. / (1.+w_fld) / a_prime_over_a;
          /** We must gauge transform the pressure perturbation from the fluid rest-frame to the gauge we are working in */
          ppw->delta_p_fld = pba->cs2_fld * ppw->delta_rho_fld + (pba->cs2_fld-ca2_fld)*(3*a_prime_over_a*ppw->rho_plus_p_theta_fld/k/k);
        }
      }'''
repl(old1,new1,"future fld stress-energy gate")

old2='''      if (pba->use_ppf == _FALSE_){

        /** - ----> factors w, w_prime, adiabatic sound speed ca2 (all three background-related),
            plus actual sound speed in the fluid rest frame cs2 */

        class_call(background_w_fld(pba,a,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, ppt->error_message);
        w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

        ca2 = w_fld - w_prime_fld / 3. / (1.+w_fld) / a_prime_over_a;
        cs2 = pba->cs2_fld;

        /** - ----> fluid density */

        dy[pv->index_pt_delta_fld] =
          -(1+w_fld)*(y[pv->index_pt_theta_fld]+metric_continuity)
          -3.*(cs2-w_fld)*a_prime_over_a*y[pv->index_pt_delta_fld]
          -9.*(1+w_fld)*(cs2-ca2)*a_prime_over_a*a_prime_over_a*y[pv->index_pt_theta_fld]/k2;

        /** - ----> fluid velocity */

        dy[pv->index_pt_theta_fld] = /* fluid velocity */
          -(1.-3.*cs2)*a_prime_over_a*y[pv->index_pt_theta_fld]
          +cs2*k2/(1.+w_fld)*y[pv->index_pt_delta_fld]
          +metric_euler;
      }'''
new2='''      if (pba->use_ppf == _FALSE_){

        if ((pba->has_smg == _TRUE_) && (a > exp(2.0))) {
          /* Future-only auxiliary-fluid isolation.  Freeze the two unused
             bookkeeping perturbation coordinates after their stress-energy
             has been removed.  The accepted a<=1 solution is untouched. */
          dy[pv->index_pt_delta_fld] = 0.;
          dy[pv->index_pt_theta_fld] = 0.;
        }
        else {
          /** - ----> factors w, w_prime, adiabatic sound speed ca2 (all three background-related),
              plus actual sound speed in the fluid rest frame cs2 */

          class_call(background_w_fld(pba,a,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, ppt->error_message);
          w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

          ca2 = w_fld - w_prime_fld / 3. / (1.+w_fld) / a_prime_over_a;
          cs2 = pba->cs2_fld;

          /** - ----> fluid density */

          dy[pv->index_pt_delta_fld] =
            -(1+w_fld)*(y[pv->index_pt_theta_fld]+metric_continuity)
            -3.*(cs2-w_fld)*a_prime_over_a*y[pv->index_pt_delta_fld]
            -9.*(1+w_fld)*(cs2-ca2)*a_prime_over_a*a_prime_over_a*y[pv->index_pt_theta_fld]/k2;

          /** - ----> fluid velocity */

          dy[pv->index_pt_theta_fld] =
            -(1.-3.*cs2)*a_prime_over_a*y[pv->index_pt_theta_fld]
            +cs2*k2/(1.+w_fld)*y[pv->index_pt_delta_fld]
            +metric_euler;
        }
      }'''
repl(old2,new2,"future fld derivative gate")

old3='''    if (ppt->has_source_delta_fld == _TRUE_) {
      _set_source_(ppt->index_tp_delta_fld) = ppw->delta_rho_fld/pvecback[pba->index_bg_rho_fld]
        + 3.*a_prime_over_a*(1.+pvecback[pba->index_bg_w_fld])*theta_over_k2; // N-body gauge correction
    }'''
new3='''    if (ppt->has_source_delta_fld == _TRUE_) {
      if ((pba->has_smg == _TRUE_) && (a > exp(2.0)))
        _set_source_(ppt->index_tp_delta_fld) = 0.;
      else
        _set_source_(ppt->index_tp_delta_fld) = ppw->delta_rho_fld/pvecback[pba->index_bg_rho_fld]
          + 3.*a_prime_over_a*(1.+pvecback[pba->index_bg_w_fld])*theta_over_k2;
    }'''
repl(old3,new3,"future fld density source gate")

old4='''    if (ppt->has_source_theta_fld == _TRUE_) {

      class_call(background_w_fld(pba,a,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, ppt->error_message);

      _set_source_(ppt->index_tp_theta_fld) = ppw->rho_plus_p_theta_fld/(1.+w_fld)/pvecback[pba->index_bg_rho_fld]
        + theta_shift; // N-body gauge correction
    }'''
new4='''    if (ppt->has_source_theta_fld == _TRUE_) {

      if ((pba->has_smg == _TRUE_) && (a > exp(2.0))) {
        _set_source_(ppt->index_tp_theta_fld) = 0.;
      }
      else {
        class_call(background_w_fld(pba,a,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, ppt->error_message);
        _set_source_(ppt->index_tp_theta_fld) = ppw->rho_plus_p_theta_fld/(1.+w_fld)/pvecback[pba->index_bg_rho_fld]
          + theta_shift;
      }
    }'''
repl(old4,new4,"future fld velocity source gate")

p.write_text(s)
print("FUTURE_AUX_FLD_DECOUPLE ln_a_switch=2.0")
