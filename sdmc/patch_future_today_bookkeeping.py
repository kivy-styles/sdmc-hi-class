#!/usr/bin/env python3
"""
Patch CLASS future-background bookkeeping so quantities named "today" remain
anchored at log(a/a0)=0 even when background_solve integrates to loga_final>0.

The future audit extends the interpolation table beyond a=1. Upstream CLASS
normally assumes the last table row *is* today, so without this patch it would
silently redefine age, conformal age, growth normalization and Omega0_m/r at
the future endpoint.  This patch evaluates those quantities at log(a/a0)=0
using a one-sided linear extrapolation from the final two accepted-side rows.
That avoids allowing a deliberately fast future-only action release to feed
back into quantities whose semantics are explicitly present-day.

Applied only in the isolated future section of the experimental workflow.
"""
from pathlib import Path

p=Path("source/background.c")
s=p.read_text()

def replace_once(old,new,label):
    global s
    if new in s:
        return
    n=s.count(old)
    if n!=1:
        raise RuntimeError(f"{label}: expected one anchor, found {n}")
    s=s.replace(old,new,1)

replace_once(
"""  /* growth factor today */
  double D_today;
""",
"""  /* growth factor today */
  double D_today;
  /* SDMC future audit: bracket and interpolation weight for log(a/a0)=0. */
  int index_today_future_audit = -1;
  double weight_today_future_audit = 0.;
""",
"declaration")

anchor="""  /** - recover some quantities today */
"""
insert="""  /* SDMC FUTURE AUDIT: when the table extends beyond a=1, the last row is
     no longer today.  Use the final two accepted-side rows only.  The second
     row index is the last log(a/a0)<=0 point, and the interpolation weight is
     allowed to exceed unity slightly because this is a one-sided extrapolation
     to the exact boundary. */
  {
    int jj_today_future_audit;
    int last_past_today_future_audit = -1;
    for (jj_today_future_audit=0;
         jj_today_future_audit<pba->bt_size;
         jj_today_future_audit++) {
      if (pba->loga_table[jj_today_future_audit] <= 0.)
        last_past_today_future_audit = jj_today_future_audit;
      else
        break;
    }
    class_test(last_past_today_future_audit < 1,
               pba->error_message,
               "SDMC future audit could not find two accepted-side rows for today");
    index_today_future_audit = last_past_today_future_audit-1;
    weight_today_future_audit =
      (0.-pba->loga_table[index_today_future_audit])/
      (pba->loga_table[index_today_future_audit+1]
       -pba->loga_table[index_today_future_audit]);
  }

  /** - recover some quantities today */
"""
replace_once(anchor,insert,"today-index insertion")

replace_once(
"""  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_tau];
""",
"""  pba->age =
    ((1.-weight_today_future_audit)
      *pba->background_table[index_today_future_audit*pba->bg_size
                             +pba->index_bg_time]
     +weight_today_future_audit
      *pba->background_table[(index_today_future_audit+1)*pba->bg_size
                             +pba->index_bg_time])/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age =
    (1.-weight_today_future_audit)*pba->tau_table[index_today_future_audit]
    +weight_today_future_audit*pba->tau_table[index_today_future_audit+1];
""",
"age bookkeeping")

replace_once(
"""    pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm]/pba->H0/pba->H0;
""",
"""    pba->Omega0_dcdm =
      ((1.-weight_today_future_audit)
        *pba->background_table[index_today_future_audit*pba->bg_size
                               +pba->index_bg_rho_dcdm]
       +weight_today_future_audit
        *pba->background_table[(index_today_future_audit+1)*pba->bg_size
                               +pba->index_bg_rho_dcdm])/pba->H0/pba->H0;
""",
"dcdm today")

replace_once(
"""    pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr]/pba->H0/pba->H0;
""",
"""    pba->Omega0_dr =
      ((1.-weight_today_future_audit)
        *pba->background_table[index_today_future_audit*pba->bg_size
                               +pba->index_bg_rho_dr]
       +weight_today_future_audit
        *pba->background_table[(index_today_future_audit+1)*pba->bg_size
                               +pba->index_bg_rho_dr])/pba->H0/pba->H0;
""",
"dr today")

replace_once(
"""  D_today = pvecback_integration[pba->index_bi_D];
""",
"""  D_today =
    (1.-weight_today_future_audit)
      *pba->background_table[index_today_future_audit*pba->bg_size
                             +pba->index_bg_D]
    +weight_today_future_audit
      *pba->background_table[(index_today_future_audit+1)*pba->bg_size
                             +pba->index_bg_D];
""",
"growth today")

replace_once(
"""  pba->Omega0_m = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_m];
  pba->Omega0_r = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_r];
""",
"""  pba->Omega0_m =
    (1.-weight_today_future_audit)
      *pba->background_table[index_today_future_audit*pba->bg_size
                             +pba->index_bg_Omega_m]
    +weight_today_future_audit
      *pba->background_table[(index_today_future_audit+1)*pba->bg_size
                             +pba->index_bg_Omega_m];
  pba->Omega0_r =
    (1.-weight_today_future_audit)
      *pba->background_table[index_today_future_audit*pba->bg_size
                             +pba->index_bg_Omega_r]
    +weight_today_future_audit
      *pba->background_table[(index_today_future_audit+1)*pba->bg_size
                             +pba->index_bg_Omega_r];
""",
"Omega0 bookkeeping")

p.write_text(s)
print("FUTURE_TODAY_BOOKKEEPING installed")
