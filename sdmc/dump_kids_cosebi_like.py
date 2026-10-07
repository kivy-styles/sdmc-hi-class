from cosmosis.datablock import names,option_section
from pathlib import Path
import json

def setup(options):
    return options.get_string(option_section,"output")

def execute(block,out):
    chi=float(block[names.data_vector,"cosebis_CHI2"])
    like=float(block[names.likelihoods,"cosebis_LIKE"])
    Path(out).write_text(json.dumps({"chi2":chi,"loglike":like})+"\n")
    print("KIDS_COSEBI_RESULT",out,chi,like,flush=True)
    return 0

def cleanup(config):
    return
