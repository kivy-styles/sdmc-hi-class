#!/usr/bin/env python3
"""
Patch CLASS future-background bookkeeping so quantities named "today" remain
anchored at log(a/a0)=0 even when background_solve integrates to loga_final>0.

The future audit extends the interpolation table beyond a=1. Upstream CLASS
normally assumes the last table row *is* today, so without this patch it would
silently redefine age, conformal age, growth normalization and Omega0_m/r at
the future endpoint. This patch keeps those quantities tied to the row nearest
loga=0 while retaining all future rows for the SDMC health audit.

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
  /* SDMC future audit: row representing the actual present epoch log(a/a0)=0 */
  int index_today_future_audit = 0;
""",
"declaration")

anchor="""  evolver_ndf15_abstol = 1e-15;

  /** - recover some quantities today */
"""
insert="""  evolver_ndf15_abstol = 1e-15;

  /* SDMC FUTURE AUDIT: when the table extends beyond a=1, the last row is
     no longer today. Find the row closest to log(a/a0)=0 once and use it for
     all quantities whose semantics are explicitly present-day. */
  {
    double abs_loga_today_future_audit = fabs(pba->loga_table[0]);
    int jj_today_future_audit;
    for (jj_today_future_audit=1;
         jj_today_future_audit<pba->bt_size;
         jj_today_future_audit++) {
      if (fabs(pba->loga_table[jj_today_future_audit])
          < abs_loga_today_future_audit) {
        abs_loga_today_future_audit =
          fabs(pba->loga_table[jj_today_future_audit]);
        index_today_future_audit = jj_today_future_audit;
      }
    }
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
    pba->background_table[index_today_future_audit*pba->bg_size
                          +pba->index_bg_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pba->tau_table[index_today_future_audit];
""",
"age bookkeeping")

replace_once(
"""    pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm]/pba->H0/pba->H0;
""",
"""    pba->Omega0_dcdm =
      pba->background_table[index_today_future_audit*pba->bg_size
                            +pba->index_bg_rho_dcdm]/pba->H0/pba->H0;
""",
"dcdm today")

replace_once(
"""    pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr]/pba->H0/pba->H0;
""",
"""    pba->Omega0_dr =
      pba->background_table[index_today_future_audit*pba->bg_size
                            +pba->index_bg_rho_dr]/pba->H0/pba->H0;
""",
"dr today")

replace_once(
"""  D_today = pvecback_integration[pba->index_bi_D];
""",
"""  D_today =
    pba->background_table[index_today_future_audit*pba->bg_size
                          +pba->index_bg_D];
""",
"growth today")

replace_once(
"""  pba->Omega0_m = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_m];
  pba->Omega0_r = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_r];
""",
"""  pba->Omega0_m =
    pba->background_table[index_today_future_audit*pba->bg_size
                          +pba->index_bg_Omega_m];
  pba->Omega0_r =
    pba->background_table[index_today_future_audit*pba->bg_size
                          +pba->index_bg_Omega_r];
""",
"Omega0 bookkeeping")

p.write_text(s)
print("FUTURE_TODAY_BOOKKEEPING installed")
