# CFDX — Suivi des tâches

Dernière mise à jour : 2026-09-19

## À faire

- [ ] Phase 0 — Data model (Case, Mesh, Field, BoundaryField, HDF5, Hash)
- [ ] Phase 1 — Import (OpenFOAM, Gmsh, meshio)
- [ ] Phase 2 — Geometry (face/cell geometry, quality, validation)
- [ ] Phase 3 — Operators (gradient, interpolation, divergence, laplacien, flux)
- [ ] Phase 4 — Scalar solver (Poisson, Laplace)
- [ ] Phase 5 — Linear algebra (CSR, CG, BiCGStab, GMRES, preconditioners)
- [ ] Phase 6 — MPI (domain decomposition, ghost cells, halo exchange, parallel HDF5)
- [ ] Phase 7 — GPU (CUDA backend, device storage, memory planner, GPU kernels)
- [ ] Phase 8 — GPU OOC (tile manager, halo manager, working sets, async transfers)
- [ ] Phase 9 — Incompressible / Module 1 (Navier-Stokes, SIMPLE, PISO, PIMPLE, Rhie-Chow)
- [ ] Phase 10 — Physics M2-M7 (turbulence, thermal/CHT, radiation, VOF, dynamic mesh, FSI)
- [ ] Infrastructure (CMake, Python API/CLI, profiling, logging, CLI commands, HDF5 inspection, schema versioning, reproductibilité, déterminisme, précision)
- [ ] Validation (tests analytiques, conservation, comparaison OpenFOAM, benchmarks, profiling)

## En cours

- Aucune tâche en cours

## Fait

- [x] Spécification technique CFDX v0.7 (sections 1 à 99)
- [x] Repository `cfdx_dev` créé sur GitHub (https://github.com/stevendaix/cfdx_dev)
- [x] Push de `spec.md` sur le remote `origin`
- [x] `SPEC_GUIDE.md` — guide de la spécification
- [x] `CHECKLIST.md` — checklist complète des tâches
- [x] `TASKS.md` — suivi des tâches (ce fichier)