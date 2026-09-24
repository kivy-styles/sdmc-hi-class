#!/usr/bin/env python3
"""
Extend the installed sdmc_v3_covariant_logstruct_audit gravity model with
controlled background-scalar initial-condition perturbations.

The underlying covariant functions are unchanged.  Only the initial value and
initial conformal-time derivative of the structural-clock scalar are varied:

    psi_ini      -> psi_ini + delta_psi
    psi'_ini     -> psi'_ini (1 + delta_v)

parameters_smg = V_offset, delta_psi, delta_v

The potential offset remains tuning parameter index 0.  This makes it possible
to test whether the chronological structural-clock trajectory is a dynamical
attractor of the already reconstructed action rather than merely one exact
solution selected by its original initial conditions.
"""
from pathlib import Path

p = Path("gravity_smg/gravity_models_smg.c")
s = p.read_text()

# Expand only the log-structural model's parameter list.
tag = 'if (strcmp(string1,"sdmc_v3_covariant_logstruct_audit") == 0) {'
i = s.index(tag)
j = s.index('if (strcmp(string1,"sdmc_v3_covariant_linear_audit") == 0)', i)
block = s[i:j]
old = "pba->parameters_size_smg = 1;"
new = "pba->parameters_size_smg = 3;"
if old not in block:
    if new not in block:
        raise RuntimeError("logstruct parser parameter-size anchor not found")
else:
    block = block.replace(old, new, 1)
s = s[:i] + block + s[j:]

# Perturb only the log-structural initial conditions.
ic_tag = "case sdmc_v3_covariant_logstruct_audit:"
i = s.index(ic_tag, s.index("int gravity_models_initial_conditions_smg("))
j = s.index("break;", i)
block = s[i:j]

old_phi = "pvecback_integration[pba->index_bi_phi_smg]=ps;"
new_phi = (
    "pvecback_integration[pba->index_bi_phi_smg]="
    "ps+pba->parameters_smg[1];"
)
old_v = (
    "pvecback_integration[pba->index_bi_phi_prime_smg]="
    "a/(SDMC_LS_T0*exp(ps));"
)
new_v = (
    "pvecback_integration[pba->index_bi_phi_prime_smg]="
    "a/(SDMC_LS_T0*exp(ps))*(1.+pba->parameters_smg[2]);"
)

if old_phi in block:
    block = block.replace(old_phi, new_phi, 1)
elif new_phi not in block:
    raise RuntimeError("logstruct phi IC anchor not found")

if old_v in block:
    block = block.replace(old_v, new_v, 1)
elif new_v not in block:
    raise RuntimeError("logstruct velocity IC anchor not found")

s = s[:i] + block + s[j:]
p.write_text(s)

print("LOGSTRUCT_ATTRACTOR_IC_PATCH")
print("parameters_smg = V_offset, delta_psi, delta_v")
