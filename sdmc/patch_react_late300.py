#!/usr/bin/env python3
from pathlib import Path
import re, shutil
src=Path("reactv2/reactions/src/BeyondLCDM.cpp")
s=src.read_text()
if '#include "late300_table.h"' not in s:
    s=s.replace('#include "SpecialFunctions.h"','#include "SpecialFunctions.h"\n#include "late300_table.h"')
# Exact accepted-action alpha histories for this dedicated build.
patterns=[
(r'inline double alphai_eft\(double a, double omega0, double alpha0, int model\)\{.*?\n\}', '''inline double alphai_eft(double a, double omega0, double alpha0, int model){
  switch(model){
    case 1: return late300_ak(a);
    case 2: return late300_ab(a);
    case 3: return late300_am(a);
    case 4: return 0.;
    case 5: return late300_m2(a);
    default: return 0.;
  }
}'''),
(r'inline double dalphai_eft\(double a, double omega0, double alpha0, int model\)\{.*?\n\}', '''inline double dalphai_eft(double a, double omega0, double alpha0, int model){
  switch(model){
    case 1: return late300_dak(a);
    case 2: return late300_dab(a);
    case 3: return late300_dam(a);
    case 4: return 0.;
    case 5: return 0.;
    default: return 0.;
  }
}'''),
(r'inline double ddalphai_eft\(double a, double omega0, double alpha0, int model\)\{.*?\n\}', '''inline double ddalphai_eft(double a, double omega0, double alpha0, int model){
  switch(model){
    case 1: return late300_d2ak(a);
    case 2: return late300_d2ab(a);
    case 3: return late300_d2am(a);
    case 4: return 0.;
    case 5: return 0.;
    default: return 0.;
  }
}''')
]
for pat,rep in patterns:
    s2,n=re.subn(pat,rep,s,count=1,flags=re.S)
    if n!=1: raise RuntimeError(f"failed alpha patch {pat}: {n}")
    s=s2
# Dedicated build uses stock EFT IDs 8/9 but exact late300 background.
s=s.replace('''case 8:
			/* EFTofDE: unscreened approx */
				return  HA( a, omega0);''','''case 8:
            /* late300 accepted-action background, unscreened nonlinear bracket */
                return late300_E(a);''')
s=s.replace('''case 9:
			/* EFTofDE: superscreened approx */
				return  HA( a, omega0);''','''case 9:
            /* late300 accepted-action background, superscreened nonlinear bracket */
                return late300_E(a);''')
s=s.replace('''case 8:
		/* EFTofDE: unscreened approximation */
			return HA1( a, omega0);''','''case 8:
        /* late300 accepted-action background derivative */
            return late300_HA1(a);''')
s=s.replace('''case 9:
		/* EFTofDE: superscreened approximation */
			return HA1( a, omega0);''','''case 9:
        /* late300 accepted-action background derivative */
            return late300_HA1(a);''')
src.write_text(s)
shutil.copy("output/react_tables/late300_table.h","reactv2/reactions/src/late300_table.h")
print("PATCH_REACT_LATE300_OK")
