# SDMC candidate late300

**Canonical name:** late300

**Historical provenance:** this candidate was generated and exact-tested under the workflow label `g022-A034`. Historical branch names, run IDs and artifact names are intentionally retained so the exact likelihood chain remains reproducible. From this point forward, the candidate itself is called **late300**.

The canonical numerical definition is stored in `sdmc/late300_candidate.json`.

## Promotion rule

A future SDMC candidate should not replace late300 merely because one likelihood improves. It must first preserve all four exact closure signs relative to optimized local021:

- fair Planck + DESI < 0;
- fair Planck + DESI + Pantheon+ < 0;
- fair Planck + DESI + Union3 < 0;
- fair Planck + DESI + DES-Y5 < 0.

Only after those gates remain closed should improvements in calibrated SH0ES, No-Slip-corrected weak lensing, BAO, BBN, age, stability margin or JWST/FRESCO be counted as a genuine successor.

## Lensing convention

The old compressed (S_8)-only warning is not the physical SDMC shear prediction. For late300, the No-Slip relation
[
\mu(z)=\Sigma(z)=M_{\rm Pl}^2/M_*^2(z),\qquad \eta=1
]
must be propagated into shear. The current repository audit reports an amplitude-level No-Slip-corrected proxy; a full tomographic DES-Y3/HSC/KiDS likelihood remains a future decisive test.
