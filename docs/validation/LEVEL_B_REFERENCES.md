# Level B — independent references

The Level-B test uses independent published data/correlations rather than closed-form CFDX-only solutions.

The turbulent-channel reference is the Moser, Kim & Mansour DNS at Re_tau = 178.12. The source reports y+ and normalized mean velocity U+; the repository keeps a compact deterministic subset for regression. The complete source dataset is available from the UT Austin turbulence database.

Source: https://turbulence.oden.utexas.edu/data/MKM/chan180/profiles/chan180.means

The flat-plate checks use standard ZPG correlations:
Cf_x = 0.664 / sqrt(Re_x) for laminar flow.
Cf_x = 0.0592 / Re_x^0.2 for the fully turbulent engineering correlation.

These are reference oracles. They are not claimed to be a resolved turbulent-channel or flat-plate CFD validation until a multi-cell CFDX driver produces the corresponding field statistics.
