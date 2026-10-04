# N11 — Numerical robustness diagnostics

## Scope of this implementation slice

This PR starts the repository-grounded N11 robustness work requested from the numerical-method maturity roadmap.

It deliberately implements diagnostics and algebraic scaling as **observational infrastructure**. It does not silently alter a matrix, relax a tolerance, change a solver, or fall back to another preconditioner.

The existing repository already provides an explicit constant null-space projector and pressure reference handling. This PR connects that capability to the new diagnostics layer rather than duplicating it.

## Implemented

### Matrix diagnostics

`cfdx/core/linalg/matrix_diagnostics.h` provides:

- structural consistency and finite-value checks;
- empty-row / empty-column detection;
- missing-diagonal and near-zero-diagonal detection;
- isolated-DOF detection;
- duplicate CSR-entry detection;
- diagonal dynamic-range evidence;
- row/column magnitude dynamic ranges;
- minimum Gershgorin diagonal-dominance margin for square matrices;
- sparsity-graph connected-component count;
- explicit pathology classification.

The dynamic-range quantities are **condition indicators**, not claims of an exact condition number.

### Explicit algebraic scaling

The API supports:

- no scaling;
- row scaling;
- column scaling;
- row+column scaling;
- symmetric diagonal scaling for square matrices.

Scaling returns the actual row/column factors and reports zero rows/columns. A zero row is never repaired by inserting an epsilon.

For

\[
A\mathbf{x}=\mathbf{b},
\]

row/column scaling is represented explicitly as

\[
D_r A D_c\hat{\mathbf{x}}=D_r\mathbf{b},
\qquad \mathbf{x}=D_c\hat{\mathbf{x}}.
\]

The implementation currently exposes the matrix factors; production solver wiring remains a separate step so that existing numerical paths are not changed implicitly.

### Null-space compatibility

The existing `NullSpaceProjector` remains the authoritative implementation for:

- explicit constant null spaces;
- projection;
- operator null-space verification;
- RHS compatibility checks.

N11 now exposes the compatibility norm through the diagnostics header.

### Failure classification

The new classifier distinguishes:

- non-finite solver evidence;
- matrix pathology;
- incompatible RHS;
- divergence;
- stagnation;
- maximum-iteration termination.

The classification is diagnostic only. It never changes the requested solver or preconditioner.

## Tests

The existing registered `test_null_space` target now also exercises:

- matrix pathology detection;
- connected-component detection;
- diagonal dynamic-range reporting;
- row/column scaling factors;
- preservation of explicit zero-row behaviour;
- null-space RHS compatibility.

This keeps the first N11 slice inside an already registered CTest target and avoids introducing a second test-registration mechanism.

## Deliberate non-goals

This PR does **not** claim that the complete N11 programme is closed.

The following remain subsequent qualification work:

1. physical-variable/reference-scale normalization;
2. unit-invariance campaigns;
3. full condition-number estimators with a clearly stated mathematical domain;
4. automatic matrix-characteristic-based solver recommendations;
5. production wiring of scaling into every linear-system class;
6. controlled aspect-ratio / skewness / non-orthogonality robustness campaigns;
7. AMG/Schur robustness campaigns using the N8 exact-oracle evidence;
8. initial-condition and coefficient perturbation campaigns;
9. comprehensive BC/null-space compatibility matrix;
10. machine-readable robustness campaign artifacts and final N11 qualification.

Those items should be implemented as separate small PRs with executable evidence rather than being represented as completed by documentation alone.

## Relationship with the existing roadmap

The current repository labels N10 as difficult-mesh robustness and N11 as conservation/boundedness. This PR therefore treats the requested conditioning/scaling work as a **cross-cutting numerical-robustness layer** that complements, rather than duplicates, the existing #423/#480 work.

In particular:

- difficult-mesh campaign ownership remains with the existing mesh-robustness work;
- conservation/boundedness remains independently diagnosed;
- N8 remains responsible for AMG/Schur implementation and solver-specific qualification;
- this layer provides common matrix evidence consumed by those campaigns.

## No silent recovery

The following behaviour remains prohibited:

```
requested AMG
    -> AMG failure
    -> silently use ILU
```

Likewise, this layer never performs:

```
|A_ii| < epsilon -> A_ii = epsilon
```

A pathological matrix is reported as pathological. Any retry or fallback must remain an explicit, configured policy with evidence.
