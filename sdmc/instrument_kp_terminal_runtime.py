#!/usr/bin/env python3
from pathlib import Path

p=Path("source/perturbations.c")
s=p.read_text()

old=r'''  double drs = pth->rs_d - pth->rs_star;
  if (!(drs > 0.)) return;

  double x = (pvecback[pba->index_bg_rs] - pth->rs_star)/drs;'''
new=r'''  double drs = pth->rs_d - pth->rs_star;
  {
    static int sdmc_kp_diag_pre = 0;
    if (sdmc_kp_diag_pre < 12) {
      fprintf(stdout,
              "KPTERM_RUNTIME_PRE drs=%.12e rs_star=%.12e rs_d=%.12e rs_now=%.12e\\n",
              drs,pth->rs_star,pth->rs_d,pvecback[pba->index_bg_rs]);
      fflush(stdout);
      sdmc_kp_diag_pre++;
    }
  }
  if (!(drs > 0.)) return;

  double x = (pvecback[pba->index_bg_rs] - pth->rs_star)/drs;
  {
    static int sdmc_kp_diag_x = 0;
    if ((x > -5.) && (x < 5.) && (sdmc_kp_diag_x < 80)) {
      fprintf(stdout,
              "KPTERM_RUNTIME_X x=%.12e rs=%.12e star=%.12e drag=%.12e\\n",
              x,pvecback[pba->index_bg_rs],pth->rs_star,pth->rs_d);
      fflush(stdout);
      sdmc_kp_diag_x++;
    }
  }'''
  {
    static int sdmc_kp_diag_terms = 0;
    if ((fabs(F) > 1.e-10) && (sdmc_kp_diag_terms < 60)) {
      fprintf(stdout,
              "KPTERM_RUNTIME_TERMS x=%.12e F=%.12e C=%.12e dlnC=%.12e Fk=%.12e drs=%.12e rs=%.12e k2=%.12e dg=%.12e tg=%.12e tb=%.12e\n",
              x,F,*C,*dlnC,*Fk,drs,pvecback[pba->index_bg_rs],k2,delta_g,theta_g,theta_b);
      fflush(stdout);
      sdmc_kp_diag_terms++;
    }
  }
}'''
if old not in s:
    raise SystemExit("helper diagnostic anchor not found")
s=s.replace(old,new,1)

old2=r'''      sdmc_kp_terminal_terms(pba,pth,ppt,pvecback,R,k2,
                             delta_g,theta_g,theta_b,
                             &kp_C,&kp_dlnC,&kp_Fk);
      dy[pv->index_pt_theta_b] += kp_Fk;
      dy[pv->index_pt_theta_g] += kp_Fk;'''
new2=r'''      sdmc_kp_terminal_terms(pba,pth,ppt,pvecback,R,k2,
                             delta_g,theta_g,theta_b,
                             &kp_C,&kp_dlnC,&kp_Fk);
      {
        static int sdmc_kp_diag_calls = 0;
        if ((fabs(kp_C-1.) > 1.e-10 || fabs(kp_Fk) > 1.e-30) &&
            (sdmc_kp_diag_calls < 60)) {
          fprintf(stdout,
                  "KPTERM_RUNTIME_CALL tca=%d rsa=%d C=%.12e dlnC=%.12e Fk=%.12e dyb_before=%.12e dyg_before=%.12e\n",
                  ppw->approx[ppw->index_ap_tca],
                  ppw->approx[ppw->index_ap_rsa],
                  kp_C,kp_dlnC,kp_Fk,
                  dy[pv->index_pt_theta_b],
                  dy[pv->index_pt_theta_g]);
          fflush(stdout);
          sdmc_kp_diag_calls++;
        }
      }
      dy[pv->index_pt_theta_b] += kp_Fk;
      dy[pv->index_pt_theta_g] += kp_Fk;'''
if old2 not in s:
    raise SystemExit("call diagnostic anchor not found")
s=s.replace(old2,new2,1)

old_fk=r'''  *Fk = (*dlnC)*theta_common
    + ((*C)*(*C)-1.)*k2*R_class/(1.+R_class)*(delta_g/4.);
}'''
new_fk=r'''  *Fk = (*dlnC)*theta_common
    + ((*C)*(*C)-1.)*k2*R_class/(1.+R_class)*(delta_g/4.);

  /* Runtime audit only: prove whether the terminal profile is actually active. */
  {
    static int sdmc_kp_diag_terms = 0;
    if ((fabs(F) > 1.e-10) && (sdmc_kp_diag_terms < 60)) {
      fprintf(stdout,
              "KPTERM_RUNTIME_TERMS x=%.12e F=%.12e C=%.12e dlnC=%.12e Fk=%.12e drs=%.12e rs=%.12e k2=%.12e dg=%.12e tg=%.12e tb=%.12e\\n",
              x,F,*C,*dlnC,*Fk,drs,pvecback[pba->index_bg_rs],k2,delta_g,theta_g,theta_b);
      fflush(stdout);
      sdmc_kp_diag_terms++;
    }
  }
}'''
if old_fk not in s:
    raise SystemExit("Fk diagnostic anchor not found")
s=s.replace(old_fk,new_fk,1)

p.write_text(s)
print("instrumented Kp terminal runtime audit")
