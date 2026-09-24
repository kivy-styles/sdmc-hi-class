#!/usr/bin/env python3
"""
Patch hi_class shooting targets so a dedicated future-background run still
interprets "today" at log(a/a0)=0 after background_solve is extended to
loga_final>0.

Only the Omega_smg and M2_today_smg target readout is changed. The target is linearly interpolated at exact log(a/a0)=0 rather than taken from the nearest future-grid row.  Ordinary
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
'''
new=r'''    case Omega_smg: {
      int jlo_today_smg = -1;
      int jj_today_smg;
      double w_today_smg, rho_today_smg, crit_today_smg;
      for (jj_today_smg=0; jj_today_smg<ba.bt_size-1; jj_today_smg++) {
        if (ba.loga_table[jj_today_smg] <= 0. &&
            ba.loga_table[jj_today_smg+1] >= 0.) {
          jlo_today_smg = jj_today_smg;
          break;
        }
      }
      class_test(jlo_today_smg < 0,
                 errmsg,
                 "future audit could not bracket log(a/a0)=0 for Omega_smg shooting");
      w_today_smg =
        (0.-ba.loga_table[jlo_today_smg])/
        (ba.loga_table[jlo_today_smg+1]-ba.loga_table[jlo_today_smg]);
      rho_today_smg =
        (1.-w_today_smg)*ba.background_table[jlo_today_smg*ba.bg_size+ba.index_bg_rho_smg]
        +w_today_smg*ba.background_table[(jlo_today_smg+1)*ba.bg_size+ba.index_bg_rho_smg];
      crit_today_smg =
        (1.-w_today_smg)*ba.background_table[jlo_today_smg*ba.bg_size+ba.index_bg_rho_crit]
        +w_today_smg*ba.background_table[(jlo_today_smg+1)*ba.bg_size+ba.index_bg_rho_crit];
      output[i] = rho_today_smg/pow(ba.H0,2) - ba.Omega0_smg;
      if (input_verbose > 2)
        printf(" param[%i] = %e, Omega_smg(a=1 exact interp) = %.3e, %.3e, target = %.2e \\n",
          ba.tuning_index_smg,
          ba.parameters_smg[ba.tuning_index_smg],
          rho_today_smg/crit_today_smg,
          rho_today_smg/pow(ba.H0,2),
          output[i]);
      break;
    }
    case M2_today_smg: {
      int jlo_today_m2 = -1;
      int jj_today_m2;
      double w_today_m2, M2_today_interp;
      for (jj_today_m2=0; jj_today_m2<ba.bt_size-1; jj_today_m2++) {
        if (ba.loga_table[jj_today_m2] <= 0. &&
            ba.loga_table[jj_today_m2+1] >= 0.) {
          jlo_today_m2 = jj_today_m2;
          break;
        }
      }
      class_test(jlo_today_m2 < 0,
                 errmsg,
                 "future audit could not bracket log(a/a0)=0 for M2 shooting");
      w_today_m2 =
        (0.-ba.loga_table[jlo_today_m2])/
        (ba.loga_table[jlo_today_m2+1]-ba.loga_table[jlo_today_m2]);
      M2_today_interp =
        (1.-w_today_m2)*ba.background_table[jlo_today_m2*ba.bg_size+ba.index_bg_M2_smg]
        +w_today_m2*ba.background_table[(jlo_today_m2+1)*ba.bg_size+ba.index_bg_M2_smg];
      output[i] = M2_today_interp - ba.M2_today_smg;
      if (input_verbose > 2)
        printf("M2(a=1 exact interp) = %e, want %e, param=%e\\n",
          M2_today_interp,
          ba.M2_today_smg,
          ba.parameters_smg[ba.tuning_index_2_smg]
        );
      break;
    }
'''
if old not in s:
    if "Omega_smg(a=1 exact interp)" in s:
        print("FUTURE_SHOOTING_TODAY already installed")
        raise SystemExit(0)
    raise RuntimeError("future shooting target anchor not found")
if s.count(old)!=1:
    raise RuntimeError(f"expected one shooting target anchor, found {s.count(old)}")
p.write_text(s.replace(old,new,1))
print("FUTURE_SHOOTING_TODAY installed")
