from cosmosis.datablock import names,option_section
from pathlib import Path
import json

def setup(options):
    return options.get_string(option_section,"output")

def execute(block,out):
    chi=float(block[names.data_vector,"cosebis_CHI2"])
    like=float(block[names.likelihoods,"cosebis_LIKE"])
    n=[int(x) for x in block["cosebis","n"]]
    vec=[]
    pairs=[]
    for i in range(1,7):
        for j in range(i,7):
            key=f"bin_{j}_{i}"
            vals=[float(x) for x in block["cosebis",key]]
            vec.extend(vals)
            pairs.append([j,i])
    payload={"chi2":chi,"loglike":like,"n":n,"pairs":pairs,"theory_vector":vec}
    Path(out).write_text(json.dumps(payload)+"\n")
    print("KIDS_COSEBI_RESULT",out,chi,like,"NTH",len(vec),flush=True)
    return 0

def cleanup(config):
    return
