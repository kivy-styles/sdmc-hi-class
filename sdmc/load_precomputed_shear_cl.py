import numpy as np
from cosmosis.datablock import option_section

def setup(options):
    path=options.get_string(option_section,"npz")
    x=np.load(path)
    return {k:x[k] for k in x.files}

def execute(block,cfg):
    sec="shear_cl"
    block[sec,"ell"]=cfg["ell"]
    block[sec,"nbin"]=6
    block[sec,"nbin_a"]=6
    block[sec,"nbin_b"]=6
    for i in range(1,7):
      for j in range(i,7):
        k=f"bin_{j}_{i}"
        block[sec,k]=cfg[k]
    return 0

def cleanup(cfg):
    return
