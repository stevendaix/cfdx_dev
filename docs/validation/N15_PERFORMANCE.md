# N15 — Performance and scalability

N15 provides a reproducible performance evidence layer for CFDX. It measures
cost, scaling and memory without converting performance measurements into
numerical-qualification claims.

The executable runner is \`scripts/n15_benchmark.py\`.

## Provenance

Every retained run records Git SHA, UTC timestamps, command, working directory,
CPU/memory information, GPU inventory when available, compiler, build type,
MPI/CUDA descriptors, mesh size and numerical settings.

Application metrics can additionally report setup, assembly, AMG setup, linear
solve, nonlinear solve, I/O, communication, iteration counts and true residuals.
The runner never invents missing internal timings.

## Single run

Use:
    python scripts/n15_benchmark.py run --case benchmarks/n15/poisson_baseline.json --output artifacts/n15/run.json -- ./build/cfdx_poisson

## Strong scaling

Use:
    python scripts/n15_benchmark.py scale --kind strong --values 1,2,4,8,16 --case benchmarks/n15/poisson_baseline.json --output artifacts/n15/strong.json -- mpirun -n \${CFDX_N15_MPI_RANKS} ./build/cfdx_poisson

The runner exports CFDX_N15_MPI_RANKS. Strong scaling keeps the problem fixed.

For reference rank p0:
S(p) = T(p0) / T(p)
E(p) = S(p) / (p / p0)

A two-point result is diagnostic only; a scaling claim needs a documented
campaign population and stable hardware.

## Weak scaling

Set mesh.cells_per_rank. The runner records cells = cells_per_rank * ranks
while increasing ranks.

## Memory

Peak child RSS is retained where the operating system exposes it. Production
campaigns should also emit component-level memory for mesh/geometry, fields,
assembled matrix, AMG hierarchy, Krylov vectors and temporary storage. MPI
campaigns must retain per-rank memory.

## AMG and solver decomposition

The application is responsible for real phase timers. N15 retains setup,
AMG setup/solve, linear iterations, nonlinear iterations and true residuals.
Iteration counts are never converted into synthetic time.

## Matrix-free versus assembled

Both modes require identical mesh, physics, numerical settings and convergence
criteria. N13 numerical-equivalence evidence is a prerequisite for comparison.

## GPU

GPU speedup requires real CUDA execution and N13 numerical equivalence. Measure
end-to-end runtime, transfers, communication, setup, I/O and memory where
available. CPU fallback is never a PASS. No CUDA hardware means BLOCKED.

## Mixed precision

Compare against an FP64 reference and retain solution/QoI differences,
conservation and true-residual evidence before reporting speedup.

## Regression

The regression command compares retained result JSON files and records both Git
SHAs. The default threshold is diagnostic rather than a universal CI gate
because shared runners, CPU frequency and MPI placement introduce variance.

## CI and limitations

Ordinary CI validates the runner and its contracts. Dedicated MPI/GPU/large
memory campaigns require appropriate hardware and retained artifacts.

This PR cannot execute a genuine production-scale MPI campaign, GPU speedup
campaign, AMG phase breakdown, matrix-free/assembled campaign or mixed-precision
qualification unless the required executable paths and hardware are available.
It therefore does not fabricate PASS results and leaves N15 PARTIAL.
