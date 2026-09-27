# Audit C15–C19 and repository cleanup

Date: 2026-09-27

This document records the current audit state after the consolidation of TVD/Rosseland work in PR #459 and the repository cleanup merged by PR #458.

## C15/C19 — CUDA stopping and breakdown semantics

### C15 — stopping criterion

The CUDA Poisson path now computes a tolerance from the relative tolerance and the right-hand-side norm:

- `tol_abs = tolerance * max(||b||, 1e-15)`
- the initial and iterative residuals are evaluated on the device;
- convergence is based on the computed residual against `tol_abs`;
- non-finite residuals are reported as `DIVERGED`.

This removes the audited absolute-residual-only criterion.

**Evidence level:** source audit only on the current repository state. The repository CI does not currently provide an nvcc build or GPU execution.

### C19 — Krylov breakdown

The CUDA CG path explicitly checks:

- `pAp > 0`;
- `pAp` finite;
- `alpha` finite;
- `rz` non-zero and finite;
- `rz_new` finite;
- residual finite after every update.

A breakdown is returned as `NOT_APPLICABLE`, with the iteration count preserved. It is therefore not silently converted to `MAX_ITER_REACHED`.

**Evidence level:** source audit only; no hardware CUDA execution is claimed.

### Required CUDA evidence still missing

The following remain separate hardware/toolchain evidence items:

1. compile with nvcc using `CFDX_ENABLE_GPU=ON`;
2. execute the CUDA Poisson solver on real NVIDIA hardware;
3. compare CUDA vs CPU solution and independent `b-Ax` residual;
4. exercise a deliberate breakdown/non-SPD case and verify the reported status;
5. verify no CPU fallback when GPU execution is explicitly requested;
6. record iterations, residual, wall time and device memory use.

No CUDA capability is marked fully validated until these are available.

## C16 — VTU

The dedicated `test_vtu_writer_validation` remains registered and checks topology, fields and metadata on a polyhedral VTU case.

**Status:** validated by repository test coverage; exact-head CI remains the acceptance evidence.

## C17 — absolute thresholds

The audit correction distinguishes absolute residuals from relative criteria and keeps physical/solver thresholds explicit rather than inflating tolerances.

**Status:** validated by the existing numerical-test contract; no tolerance relaxation is introduced by #459.

## C18 — MPI exception/control flow

The existing MPI validation exercises the relevant try/catch/control-flow path with multiple ranks.

**Status:** validated by the existing MPI test evidence; no CUDA dependency.

## C19 relation to C15

C15 is the stopping-criterion contract; C19 is the breakdown/status contract. They must remain independently testable. A finite residual below tolerance is convergence; a non-finite or invalid Krylov denominator is a solver breakdown, not a successful iteration-budget exhaustion.

---

# Cleanup status

## D11 — tracked build artefacts

A repository-tree audit on current `master` finds no tracked `build*/`, object, shared-library, log, temporary or CMake build artefacts matching the historical D11 patterns.

**Status:** cleanup completed/verified at repository-tree level. `.gitignore` remains the preventive control.

## D12 — dead/duplicate code

PR #458 removed the previously identified orphaned legacy files, including obsolete solver/preconditioner sources, obsolete visualization files and unregistered/broken tests.

The cleanup deliberately retained active files still referenced by CMake or current architecture, including:

- `equation_of_state.h`;
- `block_preconditioner.h`;
- `execution_policy.h`;
- `mesh_partitioner.h`;
- `local_time_stepping.h`.

No numerical implementation is changed by the cleanup.

**Status:** merged in #458.

## D13 — documentation

The historical audit contains older snapshots and must not be interpreted as the current source of truth when its commit references pre-date #458/#459.

This document is the current focused evidence for C15–C19 and D11/D12. Issue #34 remains the roadmap source of truth and should be synchronized after the corresponding PR merges.

## Acceptance rule

- No CUDA hardware evidence is claimed from CPU-only CI.
- No tolerance is weakened to obtain a green result.
- No test is disabled.
- No silent GPU-to-CPU fallback is accepted.
- Numerical changes and repository cleanup remain separately auditable even when consolidated in one PR.
