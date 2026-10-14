# SDMC split-sector validation status

This file records the next validation stage for the `sdmc-boltzmann` branch.

The current split-sector design keeps the frozen SDMC background unchanged while separating the early tracker perturbations from the late Horndeski No-Slip scalar. The early tracker is carried by a standard fluid with the frozen homogeneous density history and `cs2_fld = 0.003`. The late scalar is dormant before the release and then follows either the NKp-v2 B0 trajectory or the frozen Kp analytic release.

Required validation sequence:

1. Compile the split-sector patch on a clean checkout.
2. Run the NKp-v2 split-sector background and native hi_class stability checks.
3. Produce TT/TE/EE, lensing and linear P(k) outputs for NKp-v2.
4. Compare NKp-v2 split-sector results with the clustered-tracker-only control and with the rejected one-field B0 completion.
5. Run the frozen Kp late-No-Slip core with the same split early tracker.
6. Only after the Kp core is stable, apply the Kp terminal acoustic operator and compute the core-to-full differential.

No likelihood result should be labeled final until the full TT/TE/EE likelihood is run against the resulting spectra.
