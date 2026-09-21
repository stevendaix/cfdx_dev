# CFDX — Suivi des tâches

Dernière mise à jour : 2026-09-21

La feuille de route opérationnelle est l'issue #34. Les statuts ci-dessous distinguent
l'existence du code de sa vérification numérique et de sa validation solver-level.

## 🟢 Terminé / intégré

- [x] HDF5 mesh/field round-trip and integrity hardening
- [x] OpenFOAM/Gmsh/meshio import foundations
- [x] Geometry and mesh validation foundations
- [x] Core FVM operators and least-squares gradient foundation
- [x] CSR / Vector / CG / BiCGStab / GMRES foundations
- [x] MPI partition test on real multiple ranks
- [x] Runtime memory planner / partition semantics
- [x] GPU execution-policy and OOC foundations
- [x] M1 incompressible foundation and physical residual diagnostics
- [x] M2 turbulence model foundations
- [x] M3 energy / CHT foundations and energy-state update
- [x] M4 radiation foundations
- [x] Fluent-like application workflow
- [x] VMFL/numerical verification reporting foundation
- [x] Versioned HDF5 schema + topology/mesh hashes
- [x] CLI: check/info/inspect/run/convert/benchmark
- [x] Python orchestration API
- [x] M5/M6/M7 kernel foundations with explicit scope

## 🟡 En cours — validation / completion

- [ ] M0 complete exit criteria: all mesh families + MMS + conservation + refinement
- [ ] Full operator verification on orthogonal/non-orthogonal/skewed meshes
- [ ] Full linear-algebra robustness matrix and preconditioner equivalence
- [ ] MPI serial/parallel equivalence, parallel HDF5, N→M restart
- [ ] Full CUDA kernels and CPU/GPU numerical equivalence
- [ ] End-to-end GPU-OOC beyond-VRAM CFD
- [ ] M1 solver-level canonical validation and additional benchmark families
- [ ] M2 solver-level RANS/wall-treatment validation
- [ ] M3 full multi-region CHT validation
- [ ] M4 coupled radiation/thermal solver validation
- [ ] VMFL001–VMFL078 executable solver-level coverage
- [ ] Performance/memory/profiling baseline
- [ ] Deterministic vs performance mode evidence
- [ ] VOF full PLIC/interface reconstruction/contact-angle solver integration
- [ ] Dynamic mesh remeshing/topology-change conservative field transfer
- [ ] Full fluid/structure solver coupling and FSI benchmark suite
- [ ] Final spec 1–100 audit

## 🔴 Explicitly not yet complete

- [ ] GUI/3-D visualization layer
- [ ] Full production-grade CUDA implementation for every FVM kernel
- [ ] Full production-grade VOF/dynamic-mesh/FSI application solvers

## Governance

Every remaining item is executed through issue #34 as:
audit → implementation → unit test → numerical/physical validation → quantitative oracle →
CI check → PR audit → merge → roadmap update.
