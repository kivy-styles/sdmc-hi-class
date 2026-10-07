#!/usr/bin/env python3
from pathlib import Path
import shutil,re

src=Path("reactv2/reactions/src/BeyondLCDM.cpp")
s=src.read_text()

if '#include "late300_table.h"' not in s:
    s=s.replace('#include "SpecialFunctions.h"','#include "SpecialFunctions.h"\n#include "late300_table.h"')

def function_span(text, signature):
    start=text.find(signature)
    if start<0: raise RuntimeError(f"signature not found: {signature}")
    brace=text.find("{",start)
    if brace<0: raise RuntimeError(f"opening brace not found: {signature}")
    depth=0
    for i in range(brace,len(text)):
        if text[i]=="{": depth+=1
        elif text[i]=="}":
            depth-=1
            if depth==0: return start,i+1
    raise RuntimeError(f"unclosed function: {signature}")

def replace_function(text,signature,newbody):
    a,b=function_span(text,signature)
    return text[:a]+newbody+text[b:]

s=replace_function(s,
    "inline double alphai_eft(double a, double omega0, double alpha0, int model)",
'''inline double alphai_eft(double a, double omega0, double alpha0, int model){
  switch(model){
    case 1: return late300_ak(a);
    case 2: return late300_ab(a);
    case 3: return late300_am(a);
    case 4: return 0.;
    case 5: return late300_m2(a);
    default: return 0.;
  }
}''')
s=replace_function(s,
    "inline double dalphai_eft(double a, double omega0, double alpha0, int model)",
'''inline double dalphai_eft(double a, double omega0, double alpha0, int model){
  switch(model){
    case 1: return late300_dak(a);
    case 2: return late300_dab(a);
    case 3: return late300_dam(a);
    case 4: return 0.;
    case 5: return 0.;
    default: return 0.;
  }
}''')
s=replace_function(s,
    "inline double ddalphai_eft(double a, double omega0, double alpha0, int model)",
'''inline double ddalphai_eft(double a, double omega0, double alpha0, int model){
  switch(model){
    case 1: return late300_d2ak(a);
    case 2: return late300_d2ab(a);
    case 3: return late300_d2am(a);
    case 4: return 0.;
    case 5: return 0.;
    default: return 0.;
  }
}''')

# Patch only the HAg/HA1g function bodies, preserving every other ReACT model.
def patch_cases(text,signature,replacements):
    a,b=function_span(text,signature)
    q=text[a:b]
    for case,repl in replacements.items():
        pat=rf'(case\s+{case}\s*:.*?)(?=\n\s*case\s+\d+\s*:|\n\s*default\s*:)'
        m=re.search(pat,q,flags=re.S)
        if not m: raise RuntimeError(f"case {case} not found in {signature}")
        q=q[:m.start()]+repl+q[m.end():]
    return text[:a]+q+text[b:]

s=patch_cases(s,"double HAg(double a, double omega0, double extpars[], int model)",{
  8:'''case 8:
            /* late300 exact accepted-action background: unscreened nonlinear bracket */
            return late300_E(a);
''',
  9:'''case 9:
            /* late300 exact accepted-action background: superscreened nonlinear bracket */
            return late300_E(a);
'''
})
s=patch_cases(s,"double HA1g(double a, double omega0, double extpars[], int model)",{
  8:'''case 8:
            /* late300 exact accepted-action time derivative */
            return late300_HA1(a);
''',
  9:'''case 9:
            /* late300 exact accepted-action time derivative */
            return late300_HA1(a);
'''
})

src.write_text(s)
shutil.copy("output/react_tables/late300_table.h","reactv2/reactions/src/late300_table.h")

# Sanity checks on patched source before compilation.
t=src.read_text()
checks=[
  '#include "late300_table.h"',
  'case 1: return late300_ak(a);',
  'case 2: return late300_ab(a);',
  'case 3: return late300_am(a);',
  'return late300_E(a);',
  'return late300_HA1(a);',
]
for x in checks:
    if x not in t: raise RuntimeError(f"missing patch marker: {x}")
print("PATCH_REACT_LATE300_OK",len(t))
