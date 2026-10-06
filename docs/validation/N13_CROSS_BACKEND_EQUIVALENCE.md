# N13 — Cross-backend numerical equivalence

N13 proves that execution paths representing the same discrete problem remain
numerically equivalent within explicit, method-specific tolerances. It does not
equate bitwise identity, convergence speed, performance, or memory use.

## Contract

For fields a and b:

- absolute L1 = sum(abs(a-b));
- absolute L2 = sqrt(sum((a-b)^2));
- absolute Linf = max(abs(a-b));
- relative L2 = L2(a-b) / max(L2(a), scale_floor);
- relative Linf = Linf(a-b) / max(Linf(a), scale_floor).

The implementation is in src/cfdx/core/numerics/backend_equivalence.h.

Default N13 numerical gates are relative L2 <= 1e-12 and relative Linf <= 1e-12.
These are equivalence gates, not global CFD convergence tolerances.

## Serial versus MPI

test_n13_mpi_equivalence runs the canonical two-cell Poisson case on two ranks and
compares the reconstructed global solution with the analytical discrete oracle
[0.25, 0.75]. It also compares normal and deterministic MPI reductions.

## Assembled versus matrix-free

test_matrix_free_fv_operator compares algebraic operator application before any
solver-level conclusion. Multiple vectors are used so the comparison is not a
single-state coincidence.

## CPU versus CUDA

test_n13_cuda_equivalence executes a real CUDA kernel and compares it with a CPU
reference. A non-CUDA build or missing CUDA device is BLOCKED, never PASS.
Functional execution and numerical equivalence are independent gates. Performance
and memory behaviour belong to N15.

## Deterministic versus normal reductions

The deterministic reduction fixes contribution order. Normal MPI reduction keeps
the MPI implementation's normal execution order. N13 verifies numerical
equivalence for the canonical reduction; it does not claim bitwise reproducibility
for normal mode.

## Restart equivalence

test_n13_restart_equivalence verifies steady field preservation, global-ID remapping,
and preservation of current/previous transient states. Multistep restart must
preserve the complete temporal history, not only the latest field.

## N to M MPI restart

The existing parallel-HDF5 path maps checkpoint data by global cell ID. When
parallel HDF5 is available, the production writer/reader tests exercise a 2 to 3
rank restart. Redistribution correctness is kept separate from GPU and performance
claims.

## Evidence

Every campaign must retain commit SHA, compiler/build, MPI rank count, CUDA device
when applicable, mesh/numerics, tolerance, field/operator metrics, conservation,
true residual, and restart metadata.

The machine-readable campaign runner is scripts/run_n13_equivalence.py.

## Definition of Done

N13 requires executable evidence for all six requested axes, no silent CPU fallback,
no tolerance relaxation used to hide a defect, and explicit BLOCKED status whenever
required hardware is unavailable.


## Evidence runner hardening

The campaign runner is an evidence gate, not only a collection of process exit codes.

It now records:

- exact Git commit SHA;
- CMake build type and C++ compiler identity/version;
- MPI implementation/version;
- CUDA device, driver and toolkit when available;
- GitHub runner provenance when executed in Actions;
- captured test diagnostics and extracted equivalence metrics;
- explicit PASS, FAIL or BLOCKED state for every required test;
- explicit production 2-rank -> 3-rank MPI/HDF5 restart campaign when parallel HDF5 is present.

Unavailable infrastructure remains BLOCKED and is never converted into PASS. CPU-only CI may therefore produce
`PASS_WITH_BLOCKED_AXES`; hardware qualification runs must use `--require-cuda` and, when applicable,
`--require-parallel-hdf5`.

The CUDA campaign is exposed through the manually dispatched workflow
`.github/workflows/n13-cuda-equivalence.yml`. It requires a GitHub GPU runner and never falls back to CPU.
