# N9.6 — Pressure–velocity qualification matrix

N9.6 defines the reproducible qualification matrix for the pressure–velocity
algorithms after N9.5. This document is a qualification contract, not a
claim that every matrix cell is already qualified.

## Status

N9 remains **PARTIAL**.

- N9.0–N9.4: completed on master.
- N9.5: Couette common-system/authoritative-flux equivalence and skew Couette
  correction are completed through merged PR #658.
- N9.6: matrix definition is now explicit; the remaining matrix cells require
  execution and retained evidence before promotion.

N1–N8 are outside this campaign and are not reopened by N9.6.

## Canonical matrix

| Case | SIMPLE | SIMPLEC | PISO | PIMPLE | Fractional Step | COUPLED |
|---|---|---|---|---|---|---|
| Couette | required | required | required | required | required | required |
| Poiseuille | required | required | required | required | required | required |
| Ghia Re=100 | required | required | required | required | required | required |
| skew Couette | required | required | required | required | required | required |

The skew Couette case is the deterministic N9.5 geometry used to verify that
the pressure-gradient operator remains consistent on a non-orthogonal mesh.
Its physical boundaries are preserved exactly.

## Per-cell acceptance contract

A cell is **QUALIFIED** only when all applicable checks pass:

1. the exact case configuration is reproducible;
2. the requested algorithm is the algorithm actually dispatched;
3. the nonlinear convergence contract is satisfied;
4. the independent true momentum residual is satisfied;
5. authoritative conservative face mass fluxes are finite and complete;
6. continuity/conservation gates are satisfied independently;
7. solution equivalence is evaluated against the declared discrete reference
   where the cell is an invariance cell;
8. analytical/reference QoI gates are satisfied where a reference exists;
9. no silent fallback or method substitution occurred;
10. raw machine-readable evidence is retained for the exact HEAD.

Residuals, conservation, equivalence and QoI are separate evidence channels.
No one may substitute for another.

## Evidence record

Each matrix cell must retain at least:

- case and mesh identifier;
- algorithm and numerical scheme;
- solver request and resolved linear plan;
- convergence status and iteration count;
- continuity and normalized continuity;
- independent true momentum residual;
- authoritative face-flux diagnostics;
- velocity/pressure equivalence metrics when applicable;
- QoI and reference error when applicable;
- mesh/refinement identifier;
- runtime and linear-iteration information when available;
- final gate verdict;
- exact commit/HEAD used for execution.

A missing or malformed record is **INCOMPLETE**, never PASS.

## Qualification states

Use only these states in reports:

- **IMPLEMENTED** — production path exists and is executable.
- **VERIFIED** — the relevant numerical/operator contract has independent tests.
- **VALIDATED** — quantitative physical/reference evidence is complete.
- **QUALIFIED** — all required gates for the matrix cell are complete and
  reproducible.

Passing a unit test or a single Couette case does not promote another matrix
cell.

## Promotion policy

The campaign must report the complete matrix, including cells that are not
yet executed or that fail. A failed cell remains evidence of the current
state and is not removed to obtain a green aggregate result.

No tolerance relaxation, iteration-limit inflation, disabled test, fixed-
iteration pass condition, silent fallback, or geometry-specific empirical
correction is permitted.

## Next execution order

1. Couette — retain the N9.5 all-algorithm evidence as the baseline.
2. Poiseuille — execute all six algorithms with the analytical pressure/flow
   reference and independent conservation checks.
3. Ghia Re=100 — execute all six algorithms only through the common N9
   contract; preserve the existing Ghia validation gates.
4. skew Couette — retain the merged N9.5 operator-consistency result and
   extend it to every algorithm only after the exact case is available in the
   campaign.
5. Promote individual cells only after exact-head CI and evidence audit.

Cavity at other Reynolds numbers and external-flow cases remain later
qualification work; they are not silently promoted by this matrix.

## Relationship to N9.5

N9.5 established the critical common-system invariants:

- SIMPLE is the declared segregated reference;
- final authoritative face mass flux is compared face-by-face;
- the normalized flux-equivalence gate is strict;
- pressure is compared modulo its additive gauge;
- the coupled pressure-gradient block uses the same distance-weighted
  interpolation as the production Gauss gradient;
- skew Couette therefore exercises the same discrete pressure operator rather
  than an arithmetic-average surrogate.

PR #658 is the root-cause correction for the skew Couette coupled path. N9.6
does not weaken or replace that gate; it uses it as the foundation for the
broader matrix.

## Definition of done for N9.6

N9.6 is complete only when every required cell in the canonical matrix has
machine-readable evidence and an auditable final verdict on the exact HEAD.
Until then N9 remains **PARTIAL**.
