# CFDX — Roadmap

## Phase 0 — Repository Bootstrap

- [x] Repository, specification, documentation and CMake
- [x] HDF5 / pybind11 infrastructure
- [x] Unit-test infrastructure

## M0 — Numerical and infrastructure foundation

- [x] M0.1 Mesh topology
- [x] M0.2 Geometry
- [x] M0.3 Mesh quality
- [x] M0.4 Fields
- [x] M0.5 Boundary fields
- [x] M0.6 Interpolation
- [x] M0.7 FVM operators
- [x] M0.8 Linear algebra
- [x] M0.9 HDF5 / export
- [ ] **M0.10 Import/export — CI validation in progress**
- [x] M0.11 MPI infrastructure
- [x] M0.12 Execution abstraction
- [x] M0.13 Temporal discretization
- [x] M0.14 Thermodynamics / transport
- [x] M0.15 Diagnostics / monitoring

### M0.10 Import/export acceptance criteria

- [x] OpenFOAM native `constant/polyMesh` reader
- [x] meshio universal bridge
- [x] Gmsh routed through meshio
- [x] Common higher-order cells reduced to corner topology
- [x] Polyhedral topology where supported by meshio
- [x] Boundary / physical-group mapping
- [x] Topology validation before acceptance
- [x] OpenFOAM regression geometry
- [x] meshio regression geometry
- [x] Cube / tetra / pyramid / wedge geometry matrix
- [ ] Release CI green
- [ ] DebugSanitizers CI green

## M1 — Incompressible laminar

- [x] Navier-Stokes momentum equation
- [x] Continuity
- [x] Diffusion/convection operator foundation
- [x] Pressure-correction matrix foundation
- [x] SIMPLE/SIMPLEC control and relaxation primitives
- [x] PISO/PIMPLE correction primitives
- [x] Rhie-Chow face-flux correction primitive
- [x] Cavity validation harness
- [x] Poiseuille analytical validation
- [x] Full nonlinear SIMPLE/PISO/PIMPLE solver loop

## M2 — Turbulence

- [x] RANS k-epsilon eddy-viscosity model foundation
- [x] RANS k-omega SST eddy-viscosity model foundation
- [x] LES Smagorinsky eddy-viscosity model
- [x] DES length-scale model foundation
- [x] Transport-equation assembly and wall treatment
- [x] Channel/flat-plate validation harness

## M3 — Thermal / CHT

- [x] Energy convection/conduction operator foundation
- [x] Conductive face heat flux
- [x] Fluid/solid interface conductance and heat flux
- [x] Full transient energy solver
- [x] Multi-region CHT coupling
- [x] Thermal validation cases

## M4 — Radiation

- [x] Blackbody and gray-surface radiation
- [x] Two-surface net radiation exchange
- [~] View-factor validation — bounds/closure plus area-weighted reciprocity; geometric view-factor computation remains open
- [x] P1 source-term primitive
- [x] P1 constant-property scalar solve
- [~] DOM quadrature validation — isotropic first/second moments are checked; angular refinement remains open
- [~] Participating-media DOM transport — gray constant-property transport is implemented; production diffuse-gray wall BCs and spatially varying optical properties remain open
- [~] Radiation/energy coupling — infrastructure and equilibrium gate exist; non-trivial heat-transfer benchmark remains open
- [~] Rosseland — conductivity helper exists; dedicated nonlinear energy solve and boundary treatment remain open
- [ ] Spectral/non-gray radiation
- [ ] General geometric S2S/view-factor computation

### Phase 12 audit
See `docs/development/PHASE12_RADIATION_AUDIT.md`. The phase is **IMPLEMENTED / VALIDATION IN PROGRESS**, not fully green. Component-level analytical checks must not be interpreted as complete solver-level radiation validation.

### M1-M4 implementation gate

The implementation stack is now present end-to-end:
- finite-volume momentum/continuity and pressure correction;
- nonlinear SIMPLE/PISO/PIMPLE iteration;
- k-epsilon and SST transport plus wall closures;
- transient energy and spatially varying thermal boundary values;
- two-region CHT interface matching and heat-flux balance;
- participating-media DOM transport;
- radiation/energy outer coupling;
- unit, analytical and solver-level regression tests.

The remaining acceptance criterion is execution of the complete CI/DebugSanitizers regression on the final PR head. Feature completion is not treated as evidence of numerical validation; benchmark tolerances remain explicit in the validation tests.


## Current audit update — 2026-09-21

The M1-M4 implementation is undergoing a full numerical and code audit. Analytical Level-A benchmarks, conservation checks, radiation normalization, pressure-correction scaling, turbulence production units, mesh I/O bounds checks, MPI halo exchange and legacy convection assembly are being hardened before the stack is considered validated.


## Phase 8 — Poisson backend completion — 2026-09-23

- [x] 8.1–8.8 CPU Poisson/Laplace analytical, boundary-condition, conservation and refinement validation
- [x] 8.9 MPI distributed Poisson backend
  - matrix-free distributed CG on owned cell unknowns;
  - explicit cell halo exchange at MPI interfaces;
  - deterministic or fast global Krylov reductions;
  - serial/MPI solution-equivalence regression on a two-rank mesh;
  - global residual verification after convergence.
- [x] 8.10 CUDA Poisson backend
  - CSR matrix and Krylov vectors resident on the GPU during iteration;
  - CUDA SpMV, Jacobi preconditioning, vector updates and device reductions;
  - explicit host/device transfer only at setup and final solution retrieval;
  - CUDA-device availability handling and numerical regression.
