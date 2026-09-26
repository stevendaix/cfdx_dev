# CFDX — Checklist complete

Organisee par phase d'apres la roadmap (§95) et les criteres de sortie du Module 0 (§96).

Legende : [ ] a faire / [~] en cours / [x] fait

---

## Phase 0 — Data model

- [x] Creer un `case.cfdx.h5` minimal (via pipeline.py)
- [x] Recharger le fichier (via pipeline.py)
- [x] Reconstruire le mesh depuis le fichier
- [x] Reconstruire les fields
- [x] Calculer les hashes (topology, mesh, case)
- [x] Verifier l'integrite (case_setup.json validates)

## Phase 1 — Import

- [x] Import OpenFOAM (constant/polyMesh : points, faces, owner, neighbour, boundary)
- [x] Import Gmsh (.msh)
- [x] Import meshio (.vtu, etc.)
- [x] Import Fluent legacy ASCII (.cas/.dat)
- [x] Import Fluent HDF5 (.cas.h5) via pyfluent CaseFile (no license)
- [x] Import STAR-CCM+ (.sim) — metadata extraction only, topology blocked
- [x] Import SU2 (.su2/.cfg)
- [x] Import Code_Saturne (.xml/.py) — config import, mesh requires neutral export
- [x] Pipeline commun de normalisation
- [x] Gap Analysis reporting (Markdown + JSON)

## Phase 2 — Geometry

- [ ] Calcul des centres de face
- [ ] Calcul des aires et vecteurs surface Sf
- [ ] Calcul des centres de cellule
- [ ] Calcul des volumes de cellule
- [ ] Calcul du skewness
- [ ] Calcul de la non-orthogonalite
- [ ] Validator de mesh (topologie + geometrie + qualite + conservation)

## Phase 3 — Operators

- [ ] Gradient (Gauss linear)
- [ ] Interpolation cell->face (linear, upwind, limited)
- [ ] Divergence
- [ ] Laplacien (orthogonal, non-orthogonal corrected, uncorrected)
- [ ] Flux
- [ ] Surface integrate
- [ ] Volume integrate
- [ ] Tests de conservation (flux = 0 sur volume fermé)

## Phase 4 — Scalar solver

- [ ] Solveur Poisson
- [ ] Solveur Laplace
- [ ] Tests analytiques (cube, tetrahedron, pyramid, prism, polyhedron arbitraire)

## Phase 5 — Linear algebra

- [ ] SparseMatrix (CSR)
- [ ] Vector
- [ ] LinearSystem
- [ ] Solveur CG
- [ ] Solveur BiCGStab
- [ ] Solveur GMRES
- [ ] Preconditionneur Jacobi
- [ ] Preconditionneur Gauss-Seidel
- [ ] Preconditionneur ILU
- [ ] Preconditionneur AMG (plus tard)

## Phase 6 — MPI

- [ ] Domain decomposition
- [ ] Ghost cells
- [ ] Halo exchange
- [ ] Parallel HDF5
- [ ] Restart independant du nombre de ranks

## Phase 7 — GPU

- [ ] CUDA backend
- [ ] Device storage
- [ ] Memory planner
- [ ] GPU kernels (gradient, divergence, laplacien, etc.)
- [ ] Full GPU mode (H2D initial, iterations GPU, D2H checkpoint)

## Phase 8 — GPU OOC

- [ ] Tile manager
- [ ] Halo manager
- [ ] Working sets
- [ ] Pinned buffer pool
- [ ] Async transfers (H2D, D2H)
- [ ] Double buffering
- [ ] Mode GPU out-of-core (probleme > VRAM)

## Phase 9 — Incompressible (Module 1)

- [ ] Navier-Stokes
- [ ] Continuity
- [ ] Couplage pression-vitesse
- [ ] Rhie-Chow
- [ ] SIMPLE
- [ ] SIMPLEC
- [ ] PISO
- [ ] PIMPLE

## Phase 10 — Physics (M2-M7)

- [ ] M2 Turbulence (RANS k-epsilon, k-omega, SST, puis LES, DES)
- [ ] M3 Thermal / CHT
- [ ] M4 Radiation
- [ ] M5 VOF
- [ ] M6 Dynamic mesh
- [ ] M7 FSI

## Infrastructure

- [ ] Build system (CMake)
- [ ] Python API / CLI
- [ ] Profiling hooks (mesh, geometry, assembly, gradient, divergence, laplacien, matrix-vector, preconditioner, solver, MPI, H2D, D2H, halo, I/O)
- [ ] Logging (ERROR, WARNING, INFO, DEBUG, TRACE)
- [ ] CLI commands (check, info, run, convert, inspect, benchmark)
- [ ] Inspection HDF5 (h5ls, h5dump)
- [ ] Versionnage du schema HDF5
- [ ] Reproductibilite (version, compiler, CUDA, MPI, GPU, CPU, precision, schemes, hashes)
- [ ] Determinisme (deterministic mode vs performance mode)
- [ ] Precision (float32, float64, mixed precision)

## Validation

- [ ] Tests analytiques gradient (constant, lineaire, quadratique)
- [ ] Tests analytiques laplacien (constant, lineaire, quadratique)
- [ ] Tests de conservation
- [ ] Comparaison OpenFOAM (mesh, BC, discretisation, parametres physiques)
- [ ] Benchmark (wall time, CPU time, memory, VRAM, bandwidth, FLOP/s, MPI, GPU utilization, PCIe)
- [ ] Memory benchmark (mesh, fields, matrix, preconditioner, temporaires, peak)
- [ ] Profiling
- [ ] Tests sur cube, tetrahedron, pyramid, prism, polyhedron arbitraire