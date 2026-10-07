from cosmosis.datablock import option_section
from pathlib import Path
import json

def setup(options):
    return options.get_string(option_section,"output")

def execute(block,out):
    sec="cosebis"
    n=[int(x) for x in block[sec,"n"]]
    vals=[]; pairs=[]
    for i in range(1,7):
      for j in range(i,7):
        x=[float(v) for v in block[sec,f"bin_{j}_{i}"]]
        vals.extend(x); pairs.append([j,i])
    Path(out).write_text(json.dumps({"n":n,"pairs":pairs,"theory":vals})+"\n")
    print("COSEBI_VECTOR_READY",out,len(vals),flush=True)
    return 0
def cleanup(config): return
