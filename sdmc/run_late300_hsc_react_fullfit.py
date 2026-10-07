#!/usr/bin/env python3
import sys, json, importlib.util
from pathlib import Path

S=importlib.util.spec_from_file_location("rx","sdmc/run_late300_hsc_react_cutscan.py")
rx=importlib.util.module_from_spec(S); S.loader.exec_module(rx)

def main():
    model_name=sys.argv[1]
    variant=sys.argv[2]
    if variant not in {"us_native","ss_native","us_screen03","ss_screen03"}:
        raise SystemExit("unsupported variant")
    m={"late300":rx.h.late,"local021":rx.h.lcdm}[model_name]
    print("HSC_REACT_FULL_FIT_START",model_name,variant,flush=True)
    obj=rx.HSCReact(m,variant)
    fit=obj.fit()
    out={
      "status":"HSC-Y3 corrected ReACT full-band NLA-z nuisance refit",
      "model":model_name,"variant":variant,"ndata":len(rx.h.DATA),
      "fit":fit,
      "missing":["full TATT A2/alpha2/bias_ta sector","baryonic-feedback nuisance"],
      "nonlinear_limits":{"ReACT_zmax":2.5,"high_z":"boost tapered to unity z=2.5..3",
                          "high_k":"boost frozen at validated ReACT kmax"}
    }
    p=Path("output/hsc_react_full");p.mkdir(parents=True,exist_ok=True)
    (p/f"{model_name}_{variant}.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("HSC_REACT_FULL_FIT_DONE",json.dumps(out,sort_keys=True),flush=True)

if __name__=="__main__": main()
