#!/usr/bin/env python3
"""
Audit the provenance of the accepted late300 linear-G3 replay.

Purpose
-------
The accepted covariant action is genuinely free-evolved after reconstruction,
but the reconstruction itself is conditioned on the late300 target.  This
script records that distinction mechanically from the repository files so that
replay freedom is not confused with model-family parameter counting.
"""
from pathlib import Path
import json, re

WF = Path(".github/workflows/sdmc_late300_linear_covariant_audit.yml")
RC = Path("sdmc/run_late300_linear_covariant.py")

wf = WF.read_text()
rc = RC.read_text()

def heredoc_after(marker):
    i = wf.index(marker)
    start = wf.index("\n", i) + 1
    m = re.search(r"^\s*EOF\s*$", wf[start:], re.M)
    if not m:
        raise RuntimeError(f"missing heredoc terminator after {marker}")
    return wf[start:start+m.start()]

target = heredoc_after("cat > output/linear_cov_target.ini")
free = heredoc_after("cat > output/linear_cov_free.ini")

def ini_value(block, key):
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(.+?)\s*$", block, re.M)
    return m.group(1).strip() if m else None

def vector(s):
    return [float(x.strip()) for x in s.split(",")] if s else None

ordinary_keys = ["H0","omega_b","omega_cdm","A_s","n_s","tau_reio","YHe"]
ordinary_target = {k: ini_value(target,k) for k in ordinary_keys}
ordinary_free = {k: ini_value(free,k) for k in ordinary_keys}
kin = vector(ini_value(target,"parameters_smg"))
exp = vector(ini_value(target,"expansion_smg"))

checks = {
    "target_is_explicitly_generated": "Generate late300 parameterized target" in wf,
    "reconstruction_reads_target_ini": 'TARGET_INI = Path("output/linear_cov_target.ini")' in rc,
    "reconstruction_reads_target_background": 'output/linear_cov_target_00_background.dat' in rc,
    "reconstruction_builds_g_k1_k2_V": all(x in rc for x in ["g=-F1/(H*H)","k2=","k1=","V="]),
    "reconstruction_freezes_expansion_constants": all(x in rc for x in [
        "double Ox =", "double lambda_e =", "double zt =", "double A =", "double B ="
    ]),
    "free_run_uses_reconstructed_action": ini_value(free,"gravity_model") == "sdmc_v3_covariant_linear_audit",
    "free_run_structural_offset_zero": ini_value(free,"parameters_smg") == "0.0",
    "shared_ordinary_values_unchanged": ordinary_target == ordinary_free,
    "target_has_kinetic_vector": kin is not None and len(kin) == 6,
    "target_has_expansion_vector": exp is not None and len(exp) == 8,
}

if not all(checks.values()):
    raise RuntimeError("provenance audit failed: " + json.dumps(checks, sort_keys=True))

kin_names = ["AF","zc","width","D0","power","Dfloor"]
exp_names = ["Omega_x0","lambda_e","z_t","dN_t","A_late","tau_A","B_late","tau_B"]

out = {
    "status": "accepted-action replay provenance audit",
    "checks": checks,
    "ordinary_cosmology_carried_into_free_replay": ordinary_target,
    "target_kinetic_coordinates": dict(zip(kin_names,kin)),
    "target_expansion_coordinates": dict(zip(exp_names,exp)),
    "reconstructed_functions": ["g(phi)","k1(phi)","k2(phi)","V(phi)","F(phi)"],
    "replay_stage": {
        "description": "After reconstruction the covariant action is held fixed and hi_class free-evolves it.",
        "structural_offset_parameter_smg": 0.0,
        "post_reconstruction_structural_retuning_detected": False,
        "interpretation": "This is a genuine frozen-action dynamical replay."
    },
    "upstream_conditioning": {
        "description": (
            "The frozen action functions are reconstructed from the parameterized late300 target "
            "background and from target kinetic/late-expansion coordinates."
        ),
        "target_background_consumed_by_reconstruction": True,
        "target_coordinates_consumed_by_reconstruction": True,
        "interpretation": (
            "Reproducing late300 after the action is frozen proves replay/dynamical consistency, "
            "but does not by itself make the upstream late300 coordinates parameter-free predictions."
        )
    },
    "complexity_implication": {
        "conditional_replay_delta_k": (
            "Zero new post-reconstruction fit coordinates are introduced by the free replay itself."
        ),
        "model_family_delta_k": (
            "Not established as zero by this replay.  A model-family complexity count must trace "
            "which target/action coordinates were fitted versus independently derived or fixed a priori."
        ),
        "next_required_closure": (
            "Derive enough of the late300 target/action coordinates from SDMC identities or specify "
            "their reproducible priors before claiming complexity-penalized or full Bayesian preference."
        )
    }
}

Path("output").mkdir(exist_ok=True)
Path("output/late300_action_complexity_provenance.json").write_text(
    json.dumps(out,indent=2,sort_keys=True)+"\n"
)
print("LATE300_ACTION_COMPLEXITY_PROVENANCE",json.dumps(out,sort_keys=True))
