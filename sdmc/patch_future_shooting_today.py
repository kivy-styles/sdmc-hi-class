#!/usr/bin/env python3
"""
Patch hi_class shooting targets so a dedicated future-background run still
interprets "today" at log(a/a0)=0 after background_solve is extended to
loga_final>0.

Only the Omega_smg and M2_today_smg target readout is changed.  Ordinary
accepted runs are completed before this patch is applied in the experimental
workflow.
"""
from pathlib import Path

p=Path("source/input.c")
s=p.read_text()

old=r'''    case Omega_smg:
      output[i] = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_smg]/pow(ba.H0,2) - ba.Omega0_smg;
      if (input_verbose > 2)
        printf(" param[%i] = %e, Omega_smg = %.3e, %.3e, target = %.2e \n",ba.tuning_index_smg,
          ba.parameters_smg[ba.tuning_index_smg],
          ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_smg]
              /ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_crit], ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_rho_smg]
              /pow(ba.H0,2),output[i]);
      break;
    case M2_today_smg:
      output[i] = ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_M2_smg] - ba.M2_today_smg;
      if (input_verbose > 2)
        printf("M2 = %e, want %e, param=%e\n",
         ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_M2_smg],
         ba.M2_today_smg,
         ba.parameters_smg[ba.tuning_index_2_smg]
        );
      break;
    }
'''
new=r'''    case Omega_smg: {
      int index_today_smg = 0;
      double abs_loga_today_smg = fabs(ba.loga_table[0]);
      int jj_today_smg;
      for (jj_today_smg=1; jj_today_smg<ba.bt_size; jj_today_smg++) {
        if (fabs(ba.loga_table[jj_today_smg]) < abs_loga_today_smg) {
          abs_loga_today_smg = fabs(ba.loga_table[jj_today_smg]);
          index_today_smg = jj_today_smg;
        }
      }
      output[i] = ba.background_table[index_today_smg*ba.bg_size+ba.index_bg_rho_smg]/pow(ba.H0,2) - ba.Omega0_smg;
      if (input_verbose > 2)
        printf(" param[%i] = %e, Omega_smg(a=1) = %.3e, %.3e, target = %.2e \n",ba.tuning_index_smg,
          ba.parameters_smg[ba.tuning_index_smg],
          ba.background_table[index_today_smg*ba.bg_size+ba.index_bg_rho_smg]
              /ba.background_table[index_today_smg*ba.bg_size+ba.index_bg_rho_crit],
          ba.background_table[index_today_smg*ba.bg_size+ba.index_bg_rho_smg]
              /pow(ba.H0,2),output[i]);
      break;
    }
    case M2_today_smg: {
      int index_today_m2 = 0;
      double abs_loga_today_m2 = fabs(ba.loga_table[0]);
      int jj_today_m2;
      for (jj_today_m2=1; jj_today_m2<ba.bt_size; jj_today_m2++) {
        if (fabs(ba.loga_table[jj_today_m2]) < abs_loga_today_m2) {
          abs_loga_today_m2 = fabs(ba.loga_table[jj_today_m2]);
          index_today_m2 = jj_today_m2;
        }
      }
      output[i] = ba.background_table[index_today_m2*ba.bg_size+ba.index_bg_M2_smg] - ba.M2_today_smg;
      if (input_verbose > 2)
        printf("M2(a=1) = %e, want %e, param=%e\n",
         ba.background_table[index_today_m2*ba.bg_size+ba.index_bg_M2_smg],
         ba.M2_today_smg,
         ba.parameters_smg[ba.tuning_index_2_smg]
        );
'''
if old not in s:
    if "Omega_smg(a=1)" in s:
        print("FUTURE_SHOOTING_TODAY already installed")
        raise SystemExit(0)
    raise RuntimeError("future shooting target anchor not found")
if s.count(old)!=1:
    raise RuntimeError(f"expected one shooting target anchor, found {s.count(old)}")
p.write_text(s.replace(old,new,1))
print("FUTURE_SHOOTING_TODAY installed")
