# CFDX user reference

This page is the entry point for task-oriented reference material. The reference is intentionally explicit about capability status.

## Data formats

The canonical semantic separation is:

- `<case>.cfdx.h5`: case/setup definition;
- `<case>.dat.h5`: numerical state/checkpoint;
- `<case>_<time>.vtu`: visualisation/output representation.

The exact schema is documented in the Theory data-model chapter and is versioned.

## Numerical methods

For each numerical method, consult the Theory chapter for mathematics, the Developer Guide for implementation contracts, and the V&V Guide for verification/qualification evidence.

## Validation cases

The benchmark catalogue includes Couette, Poiseuille, Ghia cavity, VMFL036, backward-facing step, NACA0012, thermal and radiation cases. A catalogue entry is not proof that the current implementation is qualified.

## Troubleshooting

Inspect the first failed invariant rather than changing several numerical settings simultaneously. Preserve a reproducible failing case and its logs.

## API and CLI

API/CLI syntax must come from the actual repository interfaces. This reference deliberately does not invent commands or Python symbols that have not been verified against the source tree.

## Documentation/evidence precedence

If documentation and executable evidence disagree, the retained executable evidence determines V&V status; the documentation must then be corrected.
