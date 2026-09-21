# CFDX — Roadmap

Basé sur la spécification CFDX v0.7 (§95) et l'agent prompt (§39).

## Phase 0 — Repository Bootstrap

- [x] Repository GitHub créé (https://github.com/stevendaix/cfdx_dev)
- [x] Spécification v0.7 commitée
- [x] docs/development créés (ROADMAP, BACKLOG, DECISIONS, AGENT_STATE)
- [x] Arborescence source créée (src/cfdx/...)
- [x] Système de build (CMake)
- [x] HDF5 statique construit (third_party/hdf5)
- [x] pybind11 installé
- [x] Tests unitaires de base

## M0.1 — Mesh topology

- [x] M0.1-T01 Point storage
- [x] M0.1-T02 Face CSR connectivity
- [x] M0.1-T03 Owner/neighbour storage
- [x] M0.1-T04 Cell CSR connectivity
- [x] M0.1-T05 Boundary patches
- [x] M0.1-T06 Topology validation

## M0.2 — Geometry

- [x] M0.2-T01 Face geometry (centre, area, Sf, normal)
- [x] M0.2-T02 Cell geometry (centre, volume)
- [x] M0.2-T03 Skewness
- [x] M0.2-T04 Non-orthogonality

## M0.3 — Mesh quality

- [x] M0.3-T01 Mesh validator (topology + geometry + quality + conservation)

## M0.4 — Fields

- [x] M0.4-T01 Field<T, Location> template
- [x] M0.4-T02 Field metadata
- [x] M0.4-T03 StorageHandle (Host/Device/WorkingSet)

## M0.5 — Boundary fields

- [x] M0.5-T01 BoundaryField / PatchField

## M0.6 — Interpolation

- [x] M0.6-T01 Cell→Face interpolation (linear, upwind, limited)

## M0.7 — FVM operators

- [x] M0.7-T01 Gauss gradient
- [x] M0.7-T02 Divergence
- [x] M0.7-T03 Laplacian (orthogonal + non-orthogonal corrected)
- [x] M0.7-T04 Flux
- [x] M0.7-T05 Surface/Volume integrate
- [x] M0.7-T06 Source term linearization (Su + Sp*φ)

## M0.8 — Linear algebra

- [x] M0.8-T01 SparseMatrix (CSR)
- [x] M0.8-T02 Vector
- [x] M0.8-T03 LinearSystem
- [x] M0.8-T04 CG solver
- [x] M0.8-T05 BiCGStab solver
- [x] M0.8-T06 GMRES(m) solver (restarted)
- [x] M0.8-T07 Preconditioners: Jacobi, Red-Black Gauss-Seidel, ILU(0)
- [x] M0.8-T08 AMG wrapper (Hypre/ML) ← CRITIQUE pour Poisson
- [x] M0.8-T09 Block preconditioners (Schur complement for coupled systems)
- [x] M0.8-T10 Matrix-free operator evaluation (prep for M0.12 GPU/large-scale)

## M0.9 — HDF5 / Export (v4 Memory-Traffic-First architecture applied)

> **Note architecture v4** : Memory Planner (lifetime + reuse + budget + ledger), Mesh Reordering (RCM + SFC), 6 KPIs (bytes/cell stored/iteration, peak RAM/VRAM/cell, CPU-GPU/MPI traffic/iteration). Voir `docs/development/BACKLOG.md` et `memory_planner.md` (§1-§4).


- [x] M0.9-T01 HDF5 writer (mesh, fields, meta)
- [x] M0.9-T02 HDF5 reader
- [x] M0.9-T03 Round-trip validation
- [x] M0.9-T04 Hash computation
- [x] M0.9-T05 Lightweight VTU writer (XML, no VTK dependency)

## M0.10 — Import/export

- [x] M0.10-T01 OpenFOAM importer
- [x] M0.10-T02 Gmsh importer
- [x] M0.10-T03 meshio importer (directory exists, implementation complete)
- [x] M0.10-T04 CFDX case writer (complete)
- [x] M0.10-A1 Memory Planner v4 (lifetime + reuse + budget + ledger)
- [x] M0.10-A3 Mesh Reordering v4 (RCM + SFC + GPUOptimized + PartitionAware)

## M0.11 — MPI

- [x] M0.11-T01 Domain decomposition (real implementation)
- [x] M0.11-T02 Ghost cells (real implementation)
- [x] M0.11-T03 Halo exchange (real implementation)

## M0.12 — Execution abstraction

- [x] M0.12-T01 CPU backend (stub)
- [x] M0.12-T02 ExecutionPolicy enum
- [x] M0.12-T03 Memory planner (v4)

## M0.13 — Temporal Discretization

- [x] M0.13-T01 Time integration schemes (Euler implicit/explicit, Crank-Nicolson, BDF2)
- [x] M0.13-T02 Local time stepping (stub) (steady-state acceleration)
- [x] M0.13-T03 Adaptive time step control (stub) (CFL-based)

## M0.14 — Thermodynamics & Transport Properties

- [x] M0.14-T01 Equation of State (stub) (Incompressible constant, Ideal Gas)
- [x] M0.14-T02 Transport models (stub) (Constant, Sutherland viscosity, Fourier conduction)
- [x] M0.14-T03 Multi-component mixture (stub) (for future scalars)

## M0.15 — Diagnostics & Execution Monitoring

- [x] M0.15-T01 Logging system (stub) (spdlog integration)
- [x] M0.15-T02 Residual monitor (stub) & convergence criteria (L2/Linf norms)
- [x] M0.15-T03 Performance profiling hooks (stub) (timers for bottlenecks)

## M1 — Incompressible laminar

- [ ] Navier-Stokes
- [ ] Continuity
- [ ] Pressure-velocity coupling
- [ ] SIMPLE / SIMPLEC / PISO / PIMPLE
- [ ] Rhie-Chow

## M2 — Turbulence

- [ ] RANS (k-epsilon, k-omega, SST)
- [ ] LES / DES

## M3 — Thermal / CHT

- [ ] Energy equation
- [ ] Conduction / convection
- [ ] Conjugate heat transfer

## M4 — Radiation

- [ ] Surface radiation
- [ ] Participating media
- [ ] View factors / DOM / P1

## M5 — VOF

- [ ] Volume fraction
- [ ] Interface reconstruction
- [ ] Surface tension / contact angle
- [ ] Compressive schemes

## M6 — Dynamic mesh

- [ ] Mesh motion / deformation
- [ ] Remeshing
- [ ] Topology changes

## M0.16 — Couplage Pression-Vitesse (SIMPLE / PISO)

- [ ] M0.16-T01 : Algorithme SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) avec sous-relaxation.
- [ ] M0.16-T02 : Algorithme PISO (pour les écoulements transitoires).
- [ ] M0.16-T03 : Correction de non-orthogonalité explicite dans l'équation de pression.

## M0.17 — Turbulence RANS

- [ ] M0.17-T01 : Modèle de viscosité turbulente (Boussinesq).
- [ ] M0.17-T02 : Modèle k-epsilon standard.
- [ ] M0.17-T03 : Modèle k-omega SST (standard industriel pour couches limites et décollements).
- [ ] M0.17-T04 : Fonctions d'amortissement de paroi (Wall functions : nutUSpaldingWallFunction, etc.).

## M0.18 — Accélération Matérielle (GPU / CUDA)

- [ ] M0.18-T01 : Abstraction des kernels (macros ou templates CPU / CUDA / HIP).
- [ ] M0.18-T02 : Portage du produit matrice-vecteur (SpMV) et des opérateurs gradient/flux sur GPU.
- [ ] M0.18-T03 : Gestion mémoire unifiée / transfert explicite Host ↔ Device pour les champs (Field SoA optimisé).

## M0.19 — AMR / Maillage Mobile

- [ ] M0.19-T01 : Raffinement / déraffinement basé sur estimateur d'erreur (gradient pression / vorticité).
- [ ] M0.19-T02 : Interpolation conservative des champs lors du remaillage.
- [ ] M0.19-T03 : Support des mailles pendantes (hanging nodes) ou raffinement octree.

## M0.20 — Robustesse Production (Restart / MMS)

- [ ] M0.20-T01 : Système Checkpoint / Restart binaire (via HDF5, état complet du solveur).
- [ ] M0.20-T02 : Suite de tests MMS automatisée (convergence L2 / L_inf, vérification ordre 2 du schéma).
- [ ] M0.20-T03 : Benchmarks de régression (cavité entraînée Re=1000, comparaison données Ghia et al.).

## M7 — FSI

- [ ] Fluid-structure interface
- [ ] Partitioned / monolithic coupling

## M8 — Reacting Flows / Combustion

- [ ] Species transport
- [ ] Finite-rate chemistry
- [ ] Flame models

## M9 — Discrete Phase Model (DPM)

- [ ] Lagrangian particle tracking
- [ ] One/two-way coupling

## M10 — Optimization / Adjoint

- [ ] Discrete adjoint differentiation
- [ ] Shape optimization