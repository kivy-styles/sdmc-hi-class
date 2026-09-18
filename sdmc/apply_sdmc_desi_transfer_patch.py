#!/usr/bin/env python3
from pathlib import Path

path = Path("source/perturbations.c")
text = path.read_text()

old = """        // Thanks to the following lines, (phi,psi) are also stored as sources
        // (Obtained directly in newtonian gauge, infereed from (h,eta) in synchronous gauge).
        // If density transfer functions are requested in the (default) CLASS format,
        // (phi, psi) will be appended to the delta_i's in the final output.
        ppt->has_source_phi = _TRUE_;
        ppt->has_source_psi = _TRUE_;
"""

new = """        // Thanks to the following lines, (phi,psi) are also stored as sources
        // (Obtained directly in newtonian gauge, infereed from (h,eta) in synchronous gauge).
        // If density transfer functions are requested in the (default) CLASS format,
        // (phi, psi) will be appended to the delta_i's in the final output.
        ppt->has_source_phi = _TRUE_;
        ppt->has_source_psi = _TRUE_;

        /*
         * SDMC/DESI diagnostic bridge:
         * expose the synchronous metric derivatives needed to reconstruct
         * the Newtonian-gauge velocity transfer,
         *
         *   alpha = (h' + 6 eta')/(2 k^2),
         *   theta_c,N = k^2 alpha,
         *   theta_b,N = theta_b,S + k^2 alpha.
         *
         * This is output-only: it does not alter the perturbation equations.
         */
        if (ppt->gauge == synchronous) {
          ppt->has_source_h_prime = _TRUE_;
          ppt->has_source_eta_prime = _TRUE_;
        }
"""

if old not in text:
    if "SDMC/DESI diagnostic bridge" in text:
        print("DESI transfer diagnostic patch already applied")
    else:
        raise SystemExit("Expected density-transfer source block not found")
else:
    path.write_text(text.replace(old, new, 1))
    print("Applied SDMC/DESI transfer diagnostic output patch")
