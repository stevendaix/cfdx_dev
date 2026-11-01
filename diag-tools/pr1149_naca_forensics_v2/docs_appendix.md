## Chronology, stencil quality and iteration-62 snapshot

Evidence-only extension of the runtime forensic instrumentation
(2026-11-01), after the N256 replay reproduced the divergence at iteration 62
with a cleanly converged pressure solve (14 Krylov iterations, relative
residual 1.6e-11) around a pressure correction of about 2.97e4. The goal is
to separate a linear-system instability from a RHS/coupling that has
already grown before the pressure solve runs.

### CFDX_NACA_MATRIX_DUMP_AT

Both NACA long tests now also set `CFDX_NACA_MATRIX_DUMP_AT=62`. The
Matrix-Market snapshot in the forensic gate keeps writing the first pressure
corrector of iteration 1 with the unchanged name `<prefix>.cells<N>.mtx` /
`.rhs.csv`, and additionally snapshots the first pressure corrector of every
iteration listed in the comma-separated variable, tagged
`<prefix>.cells<N>.iter62.mtx` / `.rhs.csv`, so the system at the reproduced
divergence can be compared with the first one. Malformed tokens are ignored
and an unset variable restores the iteration-1-only behaviour; no existing
trace is renamed, so no consumer needs adaptation.

### NACA_LS_STENCIL_QUALITY / NACA_LS_STENCIL_WORST

At iteration 1, corrector 0, inside the existing forensic gate, the solver
reconstructs the one-ring internal-neighbour least-squares stencil used by
`compute_gradient_least_squares` under `EXCLUDE_BOUNDARY` and reports the
pure geometric quality of the 3x3 normal equations with zero-valued samples
(values cannot influence rank or conditioning):

    NACA_LS_STENCIL_QUALITY cells=<nc> rank3=<n3> rank2=<n2> rank1=<n1>
        rank0=<n0> samples_min=<s> pivot_ratio_max=<r>
        pivot_ratio_max_cell=<c> pivot_ratio_max_xyz=<x>:<y>:<z>
    NACA_LS_STENCIL_WORST cell=<c> pivot_ratio=<r> xyz=<x>:<y>:<z>

followed by the eight worst-conditioned cells. The NACA meshes are
two-dimensional and coplanar, so rank 2 is the expected stencil and
`StencilQuality.condition_estimate` (finite only for full rank 3) would be
infinite everywhere; the reported pivot ratio `max_pivot/min_pivot` over
the accepted pivots is the in-plane conditioning of the pivoted normal
equations and is meaningful for planar stencils.

### Chronology analyzer

`scripts/validation/analyze_naca_pressure_chronology.py` parses
`NACA_PRESSURE_CORRECTION` chronologies and reports, per metric
(`rhs_linf`, `correction_linf`, `pressure_before_linf`,
`rhs_diagonal_ratio_linf`), the first consecutive pair whose growth exceeds
`--ratio` (default 10) as `FIRST_GROWTH`; the `LARGEST_CLEAN_SOLVE`
(converged status with relative residual <= 1e-8 maximising
`correction_linf`) — a large clean solve means the linear system itself is
well behaved and the RHS/coupling was already extreme before Krylov ran;
and a per-mesh (`rows`) summary. Usage:

    python3 scripts/validation/analyze_naca_pressure_chronology.py <log> --ratio 10

The analyzer is registered in CTest next to the spectrum analyzer and
covered by `tests/python/test_analyze_naca_pressure_chronology.py`.

Everything here is diagnostic evidence: no tolerance, gate, oracle,
convergence criterion or numerical path is modified, and the overrides
remain local to the two NACA long tests; other CFD validation cases do not
inherit them.
