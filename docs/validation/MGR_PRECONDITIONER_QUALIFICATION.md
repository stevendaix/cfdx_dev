# Native MGR-style reduction preconditioner qualification

## Scope

CFDX contains a dependency-free MGR-style reduction preconditioner inspired by the
reduction hierarchy exposed by HYPRE MGR. It is not a HYPRE runtime wrapper and
does not claim feature or performance parity.

For an explicit fine/coarse partition,

A = [A_FF A_FC; A_CF A_CC],  S~ = A_CC - A_CF M_F^-1 A_FC.

The implementation supports two F-relaxation choices:
- diagonal/Jacobi: M_F = diag(A_FF);
- ILU(0): M_F is the existing CFDX ILU(0) factorization of the F block.

HYPRE's current MGR documentation explicitly exposes per-level C/F definitions,
multiple reduction levels and configurable F-relaxation, including ILU/direct
variants. CFDX now mirrors these concepts at the native API level while keeping
the implementation smaller and dependency-free. citeturn2search4turn2search5

## Six-point qualification campaign

### 1. Couette: MGR versus existing Block-Schur

test_phase9_acceptance now runs both:
- COUPLED/BlockSchur/upwind/bounded;
- COUPLED/MGR/upwind/bounded.

Both use the same 4N coupled matrix and the same physical acceptance gates:
analytical velocity profile, Umax, transverse velocity, pressure uniformity,
conservative continuity, corrected-flux continuity and momentum residuals.

The test also reports iteration counts and algorithm-invariance max_abs_dU.
MGR is not made the automatic production default by this PR.

### 2. Quantitative qualification

Acceptance evidence includes independently:
- true b-Ax residual;
- conservative continuity;
- momentum-equation residual, including component and location diagnostics;
- analytical solution error where an oracle exists;
- Krylov/linear iteration counts;
- absence of hidden fallback or hierarchy rebuild.

A converged Krylov status alone is not an acceptance criterion.

### 3. Stronger F-relaxation

The native MGR API accepts FineRelaxation::Diagonal or FineRelaxation::ILU0.
ILU(0) reuses the existing CFDX factorization implementation; MGR does not
duplicate an ILU kernel.

The unit qualification compares both modes on a genuinely coupled F block and
records their true residuals.

### 4. Multilevel MGR

The constructor accepts additional local C/F partitions. Each additional
partition recursively reduces the previous coarse operator. The final retained
operator is handled by native AMG.

This is a configurable native reduction hierarchy, not an automatic coarsening
heuristic. Each supplied hierarchy level is validated.

### 5. Nested/generic block MGR

The implementation uses arbitrary index sets rather than hard-coded
velocity/pressure block types. The coupled solver currently supplies the natural
4N partition F=(u_x,u_y,u_z), C=p.

The unit campaign also exercises a second C/F level on a generic coupled system.
This is the foundation for later velocity/pressure/temperature nesting without
duplicating the reduction algorithm.

PETSc's FieldSplit API similarly treats velocity, pressure and temperature as
arbitrary fields/splits and permits explicit index sets, supporting this generic
design. citeturn1search2

### 6. MPI/GPU boundary

The native MGR implementation remains CPU/local-CSR only in this PR.

The existing CFDX MPI stack has partition/halo infrastructure, and the GPU stack
has explicit CUDA execution paths. MGR does not silently fall back to CPU or to a
different preconditioner when those paths are requested. Distributed coarse
assembly, coarse-process reduction, device kernels and CPU/GPU numerical
equivalence remain explicit follow-up work.

This boundary is intentional: CPU emulation is not evidence of CUDA equivalence,
and a local preconditioner test is not evidence of MPI partition invariance.

## Lifecycle contract

- setup() builds the reduction hierarchy and F-relaxation state.
- update_values() refreshes numerical values only when relevant CSR patterns are
  unchanged.
- A changed reduced or ILU F-operator pattern is rejected explicitly and
  requires setup().
- No hierarchy is silently rebuilt from update_values().

## Current coupled policy

Automatic coupled solver selection remains CoupledBlockSchur. MGR is an
explicitly requested alternative until the complete Couette and broader coupled
campaign provides quantitative evidence.

PETSc documents that Schur-preconditioner quality depends strongly on the
approximation used for the eliminated block; an explicit diag(A_FF) Schur
approximation is only effective when that diagonal is a good approximation to
A_FF. citeturn1search0turn1search3

## Acceptance status

Keep this PR unmerged until exact-head CI is green and the validation-total
campaign completes with the MGR Couette case passing all physical gates.