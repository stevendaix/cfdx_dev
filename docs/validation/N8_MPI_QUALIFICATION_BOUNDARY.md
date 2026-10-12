# N8 MPI Qualification Boundary

## Scope

N8 does not claim distributed MPI qualification for the AMG or coupled-Schur production preconditioners in the current tree.

The solver model catalogue is the authoritative capability declaration:

- CG, GMRES and FGMRES are declared MPI-capable.
- Native AMG and Smoothed Aggregation AMG are serial-only.
- Coupled Block Schur, PCD and MGR are serial-only.
- The MPI boundary therefore covers distributed Krylov/reduction infrastructure, not distributed AMG/Schur production solves.

This is a qualification boundary, not an implementation-quality claim. A serial-only declaration must not be interpreted as evidence that the method is numerically wrong; it means that CFDX has not qualified that production path under MPI.

## Evidence

The merge-gated test_n8_mpi_boundary is launched with exactly two MPI ranks when MPI is enabled. It checks:

1. the declared MPI capability of CG/GMRES/FGMRES;
2. the explicit serial-only boundary for Native AMG, Smoothed Aggregation AMG, Coupled Block Schur, PCD and MGR;
3. deterministic MPI reduction on the same communication primitive used by the existing MPI numerical-equivalence tests.

The existing test_phase8_mpi_poisson and test_n13_mpi_equivalence remain the numerical MPI evidence for the supported distributed Poisson/Krylov and reduction paths.

The N8 Schur and AMG qualification tests remain serial evidence. They are not relabelled as MPI-qualified by this test.

## Non-claims

This boundary does not claim:

- MPI convergence equivalence for AMG or Schur;
- MPI scaling for AMG or Schur;
- distributed PCD/MGR/Block-Schur qualification;
- production MPI qualification merely because a Krylov reduction is distributed.

Any future distributed AMG/Schur work must add an actual MPI execution and independently recomputed numerical-equivalence evidence before changing the capability declaration.
