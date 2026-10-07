from cosmosis.datablock import option_section
import json

def setup(options):
    return options.get_string(option_section,"label",default="model")

def execute(block,label):
    chi=float(block["data_vector","cosebis_CHI2"])
    like=float(block["likelihoods","cosebis_LIKE"])
    print("KIDS_NATIVE_WEYL_LIKELIHOOD",json.dumps({"label":label,"chi2":chi,"loglike":like},sort_keys=True),flush=True)
    return 0

def cleanup(config): pass
