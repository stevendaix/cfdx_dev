# PCD Schur approximation — N8 implementation contract

## Scope

This document defines the CFDX algebraic implementation of the pressure-convection-diffusion (PCD) Schur inverse approximation.

For the coupled pressure-velocity operator

    A = [ Auu  G ]
        [ D    C ],

the exact pressure Schur complement is

    S = C - D Auu^{-1} G.

For incompressible flow, the PCD approximation used here applies a single
pressure inverse,

    S^{-1} ~= -Fp^{-1}.

This supersedes an earlier contract, `S ~= -Kp Mp^{-1} Fp` and therefore
`S^{-1} ~= -Fp^{-1} Mp Kp^{-1}`, which applied two inverse pressure operators
where one is required. That contract, the measurements that refuted it, and the
gauge defect found alongside it are retained below as evidence.

The implementation is an algebraic Schur approximation. It does not claim that
a particular pressure discretisation is universally correct.

## Operator contract

The caller supplies three pressure-space operators:

- Mp: pressure mass operator;
- Kp: pressure diffusion/Laplacian operator;
- Fp: pressure convection-diffusion operator.

The caller also supplies independent solves for Kp and Fp.

The apply sequence is exactly:

1. pressure = -Fp^{-1} rhs.

**Exactly one pressure inverse is applied.** No dense inverse is formed and no
pressure operator is inferred from Auu, G or D.

This separation is intentional: PCD depends on the pressure convection and
diffusion discretisations, so constructing those operators implicitly from the
saddle-point blocks would hide a numerical assumption.

### Why a single inverse

The single inverse is required by the production algebra, not chosen for
convenience. The previously documented composition

    S^{-1} ~= -Fp^{-1} Mp Kp^{-1}

applies two inverse pressure operators where the Schur complement requires one.
On a uniform mesh Mp is the uniform scalar diag(cell_volumes), so that
composition collapses to a uniform factor times Fp^{-1} Kp^{-1}: a squared
pressure inverse. Measured on the production 512-unknown coupled matrix, it
produced an action of norm 51510.5 against an exact-Schur action of norm 0.70831,
with a cosine of -0.029 between them, and the coupled solve stalled at a relative
residual of 0.442.

Composition comparison, each substituted into the real production path, relative
residual after 2000 iterations:

| Schur inverse action | with gauge pin | without gauge pin |
| --- | --- | --- |
| `-Fp^{-1} Mp Kp^{-1}` (previous) | 0.442160 | 0.011159 |
| `-Mp^{-1} Fp Kp^{-1}` | 0.441942 | not converged |
| `-Fp^{-1} Kp^{-1} Mp` | 0.442160 | 0.011159 |
| `-Kp^{-1}` (pure diffusion) | 0.441942 | converged, 17 iterations |
| `-Fp^{-1}` | 0.441942 | converged, 17 iterations |
| `-Mp^{-1}` (diagonal) | 0.441983 | 0.003620 |

The governing rule is the number of inverses, not which operators sit between
them. Every double-inverse ordering fails; every single-inverse ordering
converges.

### The ordering is provisional

`-Fp^{-1}` is used because it retains the convection-diffusion content that
distinguishes PCD from a pure pressure-Laplacian preconditioner. **The choice
between `Fp` and `Kp` as the single inverse is not yet qualified.**

Both converge the production acceptance case, because the mass flux is zero at
its first coupled iteration, which leaves Fp and Kp nearly proportional. That
near-equality is a property of that case and is not evidence that the ordering is
immaterial: once the mass flux is non-zero, Fp and Kp differ, and
`pcd_single_inverse_ordering_is_not_vacuous_under_convection` in
`tests/unit/test_pcd_schur.cpp` measures that separation rather than assuming it.

Accordingly the claim is:

> A single-inverse action is required by the production algebra. The exact PCD
> operator ordering remains to be qualified on a non-zero-convection case.

## Sign convention

CFDX defines the saddle-point Schur complement as

    S = C - D Auu^{-1} G.

For the usual incompressible convention with C = 0, the Schur complement is
negative pressure-Laplacian-like. The PCD inverse consequently carries the
explicit leading minus sign.

The unit test does not check the action against the contract above, because a
test written from the contract cannot detect a wrong contract. It checks it
against an oracle built from the definition of S, in
`pcd_approximates_exact_schur_better_than_a_double_inverse`: the Schur
complement is formed as S = C - D Auu^{-1} G, inverted densely, and the PCD
action is compared against S^{-1} rhs and against the double-inverse
composition. Measured on that construction: single-inverse relative error
0.621938, double-inverse relative error 4.69451.

## Lifecycle

setup() validates:

- saddle-point block dimensions;
- pressure operator dimensions;
- finite pressure-operator values;
- presence of both pressure solve callbacks;
- pressure null-space dimension when a policy is configured.

update_values(blocks) preserves the common SchurApproximation contract.

Because BlockOperator intentionally contains only the saddle-point blocks,
PCD exposes an explicit update_pressure_values(...) overload for its three
additional pressure operators. This path:

- accepts value-only updates when the CSR graphs are unchanged;
- rejects any pressure-operator graph change;
- updates the retained operator references only after all validation succeeds.

A symbolic graph change therefore requires a fresh setup()/new PCD object
rather than silently reusing stale symbolic assumptions.

## Null-space policy

When a pressure null-space projector is supplied:

- an incompatible input right-hand side is rejected;
- compatible input is projected before the pressure solves;
- intermediate and final pressure vectors are projected;
- no incompatible component is silently manufactured or accepted.

The policy is optional because not every pressure formulation has the same
gauge/null-space structure.

## Reference-cell gauge

The reference cell is imposed on the pressure operators themselves: the caller
replaces that operator row with a unit row, so `Kp` and `Fp` are nonsingular and
their inverse action is defined for any right-hand side. The action must
therefore return a gauge component when its input carries one.

The action does not pin the reference cell on input or output. Pinning it is a
defect, not a policy: because no preconditioned Krylov vector can then carry a
gauge component, the coupled gauge row `p_ref = reference_value` can never be
corrected, and its residual stays frozen at its initial magnitude forever. On
the production case the iterate starts from `p.fill(0.5)` while the reference
value is 0, and the resulting residual breakdown was:

    total_norm=0.500247  gauge_row=384  gauge_rank=0  gauge_resid=0.5
    rank=0 row=384 block=GAUGE abs=0.5
    rank=1 row=392 block=pres  abs=0.0134799

`0.5` is exactly the initial pressure value, and it dominated the total residual.
This was independent of the composition: with the pin in place every candidate
ordering, including a diagonal scaling, stalled at approximately 0.442.

`pcd_gauge_component_is_not_frozen` in `tests/unit/test_pcd_schur.cpp` is the
regression guard: it requires a non-zero input gauge component to produce a
non-zero output gauge component.

## Current verification

`tests/unit/test_pcd_schur.cpp` provides deterministic unit coverage for:

- the single-inverse action, and its separation from the double-inverse
  composition it previously implemented;
- linearity of the action;
- the reference-cell gauge not being frozen;
- agreement with an oracle built from the definition of S, independent of the
  PCD contract;
- measurability of the Fp-versus-Kp ordering once convection is present;
- method name/registration surface;
- dimension rejection;
- pressure graph-change rejection;
- value-only pressure updates;
- null-space compatibility.

This is implementation and algebraic verification, not production qualification.

## Production path evidence

The production coupled path is exercised by `test_n8_pressure_velocity_matrix`
through the `COUPLED/PCD/upwind/bounded` case, which carries the same physical
gate set as the `BlockSchur` and `MGR` cases.

It failed when first exercised. Two independent defects were found and corrected;
the evidence for both is retained below rather than discarded, because the
failure mode is not recoverable from the passing run alone.

1. The production branch was unreachable. `preconditioner_model_catalog()`
   marked `pcd` as `ModelAvailability::Planned`, and
   `select_linear_solver` rejects a request whose model is not `Available`
   before it reaches the coupled guard that admits `pcd`. The wiring compiled but
   no request could ever select it. The catalog entry is now `Available` so the
   branch is reachable and measurable. That flag is the dispatch gate; it is not
   a maturity claim. `numerical_method_registry.h` still records
   `preconditioner.pcd` as `Planned` and `NUMERICAL_METHOD_CAPABILITY_MATRIX.json`
   still records `status: partial`. Neither is advanced by the dispatch change.

2. With the branch reachable, the coupled solve does not converge:

       solve_steady_incompressible: coupled momentum-continuity solve did not
       converge (status=max_iter_reached, iterations=2000,
       residual=0.500000, relative=0.441942)

   The preconditioner applies successfully on every Krylov vector; the status is
   `MAX_ITER_REACHED`, not `NOT_APPLICABLE`, so this is stagnation rather than a
   rejected operator.

   The PCD branch does not call `setup(A)` on the preconditioner itself, unlike
   the BlockSchur branch. The setup is performed by the Krylov driver instead:
   the production call resolves to the `SparseMatrix` overload of `solve_fgmres`
   (`gmres_solver.h`), whose `setup_preconditioner` parameter defaults to `true`
   and which invokes `preconditioner->setup(A)` before the first apply. This is
   checkable rather than assumed. An uninitialised `PcdSchurApproximation` has a
   null `blocks_`, so `apply()` returns false, `CoupledBlockSchurAMGPreconditioner::apply()`
   propagates that false, and the driver returns `SolverStatus::NOT_APPLICABLE`
   with the iteration count at the point of failure. That would appear as
   `status=3` with a near-zero iteration count, not `status=1` at 2000.

   Measured, on the same binary, before and after adding an explicit
   `setup(A)` plus failure check to the PCD branch:

       status=1 iterations=2000 residual=0.500247 relative=0.442160
       status=1 iterations=2000 residual=0.500247 relative=0.442160

   Bit-for-bit identical. Had the driver not been setting the preconditioner up,
   the explicit call would have changed the outcome from `status=3` to
   `status=1`. The explicit call has been kept because it is what makes a PCD
   setup failure report `PcdSchurApproximation::last_error()` instead of an
   opaque `NOT_APPLICABLE`; it is a diagnosability fix, not the enabler.

### Localisation

Two independent defects, measured on the assembled 512-unknown production
coupled matrix (8x16 cells).

**Defect 1 — the gauge pin freezes the reference pressure.** `apply()` pins both
its input right-hand side and its output pressure at
`pressure_reference_cell_`. The coupled system fixes that component by the
gauge row `p_ref = reference_value`, and it is reachable: the Krylov iterate
starts from `p.fill(0.5)` while `reference_value` is 0. A preconditioner that
never returns a gauge component therefore leaves the gauge-row residual frozen
at its initial magnitude forever. Measured, with the pin removed and the
arithmetic otherwise unchanged:

| | relative residual after 2000 iterations |
| --- | --- |
| pinned | 0.442160 |
| unpinned | 0.011159 |

The residual is entirely this row. Row-by-row breakdown of the final residual:

    total_norm=0.500247  gauge_row=384  gauge_rank=0  gauge_resid=0.5
    rank=0 row=384 block=GAUGE abs=0.5
    rank=1 row=392 block=pres  abs=0.0134799
    rank=2 row=393 block=pres  abs=0.00257616

`0.5` is exactly the initial pressure value. No momentum row and no other
pressure row contributes comparably.

**Defect 2 — the action is a squared inverse, not a Schur inverse.** The
documented action evaluates `-Fp^{-1} Mp Kp^{-1}`. Since `Mp` is
`diag(cell_volumes)`, on this mesh `Mp = (1/128) I`, so the composition
reduces to a uniform factor times `Fp^{-1} Kp^{-1}`: two inverse pressure
operators where the Schur complement requires one. Measured against
`ExactSchurApproximation` on the same matrix with an independently converged
`Auu` solve:

| quantity | value |
| --- | --- |
| `\|\|exact Schur inverse action\|\|` | 0.70831 |
| `\|\|PCD action\|\|` | 51510.5 |
| relative error | 7.27e4 |
| cosine between the two actions | -0.0292 |

The magnitude gap is explained by the double inverse, not primarily by the
viscosity sensitivity of `Fp`: applying the inverse of a mildly
ill-conditioned pressure operator twice amplifies by roughly its condition
number. The two inner solves are not implicated — `Kp` and `Fp` both converge,
at 75-82 iterations to a relative true residual of 6e-11 to 1e-10, and the
assembled operators are well formed (`Kp` diagonal 1 to 5, `Fp` diagonal 0.25
to 1, both 128x128 with 590 nonzeros).

**Composition comparison.** Each candidate was substituted into the real
production path, with and without the gauge pin, and the full coupled solve was
run. Relative residual after 2000 iterations:

| Schur inverse action | pinned | unpinned |
| --- | --- | --- |
| `-Fp^{-1} Mp Kp^{-1}` (implemented) | 0.442160 | 0.011159 |
| `-Mp^{-1} Fp Kp^{-1}` | 0.441942 | not converged |
| `-Fp^{-1} Kp^{-1} Mp` | 0.442160 | 0.011159 |
| `-Kp^{-1}` (ideal incompressible Schur inverse) | 0.441942 | **converged, 17 iterations, all gates pass** |
| `-Fp^{-1}` | 0.441942 | **converged, 17 iterations, all gates pass** |
| `-Mp^{-1}` (diagonal) | 0.441983 | 0.003620 |

Two conclusions follow, and they are different from the ones the failed run
suggested:

- With the pin in place, **every** composition stalls at approximately 0.442,
  including a diagonal scaling. The composition is therefore not the dominant
  cause; the gauge pin is.
- With the pin removed, the **single-inverse** forms `-Kp^{-1}` and `-Fp^{-1}`
  converge in 17 iterations and pass the complete physical gate set, matching
  `BlockSchur` and `MGR`, with algorithm invariance `max_abs_dU = 1.66e-07`.
  Every **double-inverse** form still fails.

The contract for this operator is therefore a single inverse. `-Fp^{-1}` is
implemented, to retain the convection-diffusion content that distinguishes PCD
from a pure pressure-Laplacian preconditioner; the ordering relative to `Kp` is
provisional, as set out above.

`-Mp^{-1} Fp Kp^{-1}` also contains two inverses and does not converge, so it is
not a valid replacement either.

No tolerance, iteration limit or restart value was altered to obtain any of
these numbers, or in the correction.

### Post-correction result on the production matrix

With the gauge pin removed and the action reduced to a single `-Fp^{-1}` inverse,
the `COUPLED/PCD` case passes every gate in the shared physical gate set:

| quantity | BlockSchur | MGR | PCD |
| --- | --- | --- | --- |
| nonlinear iterations | 17 | 17 | 17 |
| coupled linear true residual (relative) | 9.59e-11 | 9.59e-11 | 9.59e-11 |
| profile L2 | 6.23e-10 | 6.23e-10 | 2.86e-09 |
| algorithm invariance `max_abs_dU` | 1.65e-07 | 1.65e-07 | 1.69e-07 |
| gates failed | 0 | 0 | 0 |

`PHASE9_ACCEPTANCE: PASS`, `MODEL_SUMMARY successful=8 failed=0`.

This is a passing case for one configuration, on one mesh, with a zero initial
mass flux. It is not qualification.

## Qualification still required

PCD must not be marked N8-qualified. It passes the coupled acceptance case on
one mesh, so the following remain open on top of the algebraic verification
above:

- qualification of the `Fp` versus `Kp` ordering on a case with non-zero mass
  flux, where the two operators genuinely differ;
- a production case with non-zero convection, which the acceptance matrix does
  not currently contain because its coupled cases start from a zero velocity
  field;
- representative Couette, Poiseuille and Ghia matrices;
- controlled skew/non-orthogonal matrices;
- independent true-residual measurement;
- FGMRES iteration and convergence evidence on the production path;
- pressure Kp/Fp solve setup/solve cost;
- memory/NNZ accounting;
- null-space/gauge coverage on production pressure systems;
- numerical-update lifecycle in the production dispatcher;
- measured comparison against the exact-Schur oracle, promoted from the
  diagnostic measurement above into an enforced gate;
- an evidence-based acceptance envelope before automatic selection.

No arbitrary approximation-quality threshold is introduced by this implementation.
