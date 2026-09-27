# Native MGR-style reduction preconditioner qualification

## Scope

CFDX now contains a dependency-free, one-level MGR-style reduction preconditioner for coupled systems. The implementation is inspired by the reduction hierarchy described by HYPRE MGR, but is not a HYPRE wrapper and does not claim feature or performance parity.

For an explicit fine/coarse partition,

A = [A_FF A_FC]\n    [A_CF A_CC]

The current coarse operator is

S~ = A_CC - A_CF diag(A_FF)^-1 A_FC.

The application is:

1. F-relaxation: z_F = diag(A_FF)^-1 r_F;
2. restrict the coarse residual: r_C - A_CF z_F;
3. solve S~ z_C with native AMG;
4. prolongate z_C and apply the F back-substitution z_F <- z_F - diag(A_FF)^-1 A_FC z_C.

## CFDX integration

The coupled 4N solver can select MGR explicitly through PreconditionerModel::MGR. The default automatic coupled policy remains the existing coupled Block-Schur preconditioner until MGR has passed the full quantitative campaign.

For the current 4N ordering, the natural first reduction is:

- F variables: u_x, u_y, u_z;
- C variables: pressure.

The reference-pressure row remains part of the assembled system, so the reduction does not change the gauge contract.

## Lifecycle contract

- setup() constructs the reduced CSR graph and AMG hierarchy.
- update_values() refreshes numeric coefficients when the reduced graph is unchanged.
- A changed reduced CSR pattern is rejected explicitly and requires a new setup().
- The implementation does not silently rebuild a hierarchy during update_values().

## Qualification

The unit test checks valid C/F partitioning, successful setup, finite non-zero application, numeric refresh with unchanged graph, and changed-topology rejection.

The next acceptance layer is a coupled CFD campaign with independent b - A x, continuity, momentum residual, solution-oracle and iteration-count diagnostics on Couette/Poiseuille/cavity or an appropriate manufactured coupled matrix.

## Current limitations

This first implementation intentionally does not claim multiple MGR reduction levels, adaptive C/F selection, F-relaxation beyond diagonal Jacobi, distributed/MPI MGR, coarse-level process reduction, GPU execution, or HYPRE/PETSc runtime compatibility. These are separate implementation and validation slices.
