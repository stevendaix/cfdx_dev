# Structured turbulence setup

The application schema mirrors the C++ turbulence families while keeping selection metadata in one explicit, machine-readable catalogue.

The selector distinguishes three states:
- solver-ready: a complete transport model is exposed by the current solver stack;
- kernel-only: algebraic/length-scale support exists but a complete solver path is not exposed;
- planned: the model is part of the architecture but is not yet implemented as a production transport path.

The GUI should never silently promote a kernel or planned model to a production solver option.

## User-facing hierarchy

The long-term setup flow is:
1. Flow regime — Laminar / RANS / LES / Hybrid RANS-LES.
2. Base model — filtered by regime and implementation status.
3. Wall treatment — selected independently when supported.
4. Corrections — rotation/curvature, production limiter, compressibility, roughness and model-specific options.
5. Transition — separate model coupling, not another turbulence-model enum.
6. Applicability — y+, dimensionality, steady/transient constraints and mesh requirements.
7. Recommendation — transparent reasons and warnings; never an implicit model change.

## Current production dropdown

The production dropdown contains only:
- Laminar
- standard k-epsilon
- RNG k-epsilon
- Realizable k-epsilon
- k-omega
- k-omega SST
- Spalart-Allmaras

Smagorinsky, WALE, DES, DDES and IDDES remain visible to the architecture as kernel-only capabilities but are intentionally excluded from the production choice list until their full solver paths and V&V are complete.

## Future extensibility

New turbulence models should normally add one catalogue entry and their capability metadata. They should not require redesigning the case-selection API. Transition and wall treatment should remain composable dimensions rather than being encoded into names such as SST_TRANSITION_WALLFUNCTION.

See docs/development/TURBULENCE_MODEL_SELECTION_AUDIT.md for the comparative audit and validation promotion rules.