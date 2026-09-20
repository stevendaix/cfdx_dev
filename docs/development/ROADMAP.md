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
- [ ] M0.8-T06 GMRES(m) solver (restarted)
- [ ] M0.8-T07 Preconditioners: Jacobi, Red-Black Gauss-Seidel, ILU(0)
- [ ] M0.8-T08 AMG wrapper (Hypre/ML) ← CRITIQUE pour Poisson
- [ ] M0.8-T09 Block preconditioners (Schur complement for coupled systems)
- [ ] M0.8-T10 Matrix-free operator evaluation (prep for M0.12 GPU/large-scale)

## M0.9 — HDF5 / Export (v4 Memory-Traffic-First architecture applied)

> **Note architecture v4** : Memory Planner (lifetime + reuse + budget + ledger), Mesh Reordering (RCM + SFC), 6 KPIs (bytes/cell stored/iteration, peak RAM/VRAM/cell, CPU-GPU/MPI traffic/iteration). Voir `docs/development/BACKLOG.md` et `memory_planner.md` (§1-§4).


- [x] M0.9-T01 HDF5 writer (mesh, fields, meta)
- [x] M0.9-T02 HDF5 reader
- [x] M0.9-T03 Round-trip validation
- [ ] M0.9-T04 Hash computation
- [x] M0.9-T05 Lightweight VTU writer (XML, no VTK dependency)

## M0.10 — Import/export

- [x] M0.10-T01 OpenFOAM importer
- [x] M0.10-T02 Gmsh importer
- [ ] M0.10-T03 meshio importer
- [ ] M0.10-T04 CFDX case writer

## M0.11 — MPI

- [ ] M0.11-T01 Domain decomposition
- [ ] M0.11-T02 Ghost cells
- [ ] M0.11-T03 Halo exchange

## M0.12 — Execution abstraction

- [ ] M0.12-T01 CPU backend
- [x] M0.12-T02 ExecutionPolicy enum (stub)
- [x] M0.12-T03 Memory planner (v4 skeleton)

## M0.13 — Temporal Discretization

- [x] M0.13-T01 Time integration schemes (Euler implicit/explicit, Crank-Nicolson, BDF2)
- [ ] M0.13-T02 Local time stepping (steady-state acceleration)
- [ ] M0.13-T03 Adaptive time step control (CFL-based)

## M0.14 — Thermodynamics & Transport Properties

- [ ] M0.14-T01 Equation of State (Incompressible constant, Ideal Gas)
- [ ] M0.14-T02 Transport models (Constant, Sutherland viscosity, Fourier conduction)
- [ ] M0.14-T03 Basic multi-component mixture (for future scalars)

## M0.15 — Diagnostics & Execution Monitoring

- [ ] M0.15-T01 Logging system (spdlog integration)
- [ ] M0.15-T02 Residual monitor & convergence criteria (L2/Linf norms)
- [ ] M0.15-T03 Performance profiling hooks (timers for bottlenecks)

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