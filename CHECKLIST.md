# CFDX — Checklist complète synchronisée avec l'issue #34

L'issue #34 est la source de vérité pour l'exécution. Cette checklist indique le niveau
d'achèvement réel et ne transforme pas une simple présence de code en validation.

Legend: [ ] missing / [~] implemented but validation incomplete / [x] validated foundation

## M0 Data / Mesh / Operators / Algebra

- [x] HDF5 round-trip and corruption tests
- [x] Versioned schema metadata
- [x] Deterministic topology and mesh hashes
- [x] OpenFOAM/Gmsh/meshio import foundations
- [x] Geometry/quality/validator foundations
- [x] Gauss and least-squares gradient foundations
- [~] Full gradient/interpolation/divergence/laplacian MMS matrix
- [~] Full mesh-family and skew/non-orthogonal convergence matrix
- [x] CSR, CG, BiCGStab, GMRES foundations
- [~] Complete preconditioner robustness/equivalence matrix

## MPI / GPU / OOC

- [x] Real multi-rank partition CI test
- [~] Ghost/halo operator equivalence
- [~] Parallel HDF5
- [~] N→M restart
- [x] ExecutionPolicy and memory-planning foundations
- [~] Complete CUDA FVM kernels
- [~] CPU/GPU equivalence
- [~] GPU-OOC end-to-end CFD beyond VRAM

## Physics

- [~] M1 incompressible + pressure-velocity coupling
- [~] M1 canonical benchmark campaign
- [~] M2 turbulence solver-level validation
- [~] M3 thermal/CHT solver-level validation
- [~] M4 radiation solver-level validation
- [x] M5 bounded VOF reference kernel
- [~] M5 full PLIC/interface/surface-tension/contact-angle solver
- [x] M6 mesh-motion/volume-gate foundation
- [~] M6 remeshing/topology-change/field-transfer solver
- [x] M7 partitioned coupling/Aitken/work foundation
- [~] M7 coupled fluid/structure benchmark

## Application / Reproducibility

- [x] Fluent-like case lifecycle
- [x] CLI check/info/inspect/run/convert/benchmark
- [x] Python orchestration API
- [x] Reproducibility manifest
- [x] HDF5 schema versioning
- [~] Full compiler/CUDA/MPI/GPU/CPU metadata persistence
- [~] Deterministic/performance evidence
- [ ] GUI

## Validation

- [x] Analytical/reference verification framework
- [x] Numerical-model verification matrix
- [x] VMFL matrix/reporting foundation
- [~] Executable VMFL solver-level coverage
- [~] Manufactured solutions and observed-order campaign
- [~] Independent OpenFOAM comparisons
- [~] Performance/memory/profiling benchmark matrix
- [ ] Final spec/architecture audit
