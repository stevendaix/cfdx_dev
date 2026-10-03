# Benchmark Validation Campaigns

Benchmarks are the physical/reference-validation phase after code and solution verification.

Each benchmark must identify an independent reference, operating conditions, uncertainty/expected scatter where available, quantity of interest, comparison metric and pre-declared acceptance criterion.

Current campaign slots:

- Couette
- Poiseuille
- Ghia lid-driven cavity
- VMFL036
- backward-facing step
- NACA0012
- thermal benchmarks
- radiation benchmarks

These names are campaign targets, not claims that the current CFDX implementation passes them. Current status must come from retained execution evidence.


## Planned next pass

The benchmark documentation is intentionally ahead of the execution results. The next implementation/evidence sequence is:

1. **Resynchronise Couette, Poiseuille and Ghia** after the latest solver/AMG changes. Historical values must not be silently reused as current evidence.
2. **Complete VMFL036** only after the native axisymmetric implementation and operator verification are available.
3. **Freeze NACA0012 and BFS reference definitions** before interpreting force/reattachment results.
4. **Close thermal/CHT balance evidence** independently from radiation.
5. **Close radiation model-space evidence** including view-factor checks and angular refinement for DOM.
6. **Expand the turbulence populations** only after the common equation-level V&V matrix is connected to Issue #473.
7. Promote results into docs/validation/15-qualification/ only after the evidence artifact exists and the declared population is complete.

This sequence follows the current qualification authority in Issue #118 and the numerical-method dependency structure in Issue #461.

## Evidence status

Documentation states such as VERIFIED_CONTRACT, EXECUTED, FAIL and NOT_RUN describe execution evidence. They are deliberately separate from PASS and QUALIFIED decisions.

A benchmark page therefore documents **what must be demonstrated**, while the retained campaign artifact determines **what was actually demonstrated**.
