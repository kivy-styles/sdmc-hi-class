from cosmosis.datablock import option_section
import numpy as np

def setup(options):
    return options.get_string(option_section,"npz")

def execute(block,config):
    a=np.load(config)
    ell=a["ell"]
    block["shear_cl","ell"]=ell
    block["shear_cl","nbin"]=6
    block["shear_cl","nbin_a"]=6
    block["shear_cl","nbin_b"]=6
    for i in range(1,7):
        for j in range(1,7):
            key=f"bin_{j}_{i}"
            if key in a:
                block["shear_cl",key]=a[key]
    return 0

def cleanup(config): pass
