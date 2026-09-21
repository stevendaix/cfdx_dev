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

- [ ] Navier-Stokes momentum equation
- [ ] Continuity
- [ ] Diffusion discretization
- [ ] Convection discretization
- [ ] Pressure-velocity coupling
- [ ] SIMPLE / SIMPLEC
- [ ] PISO / PIMPLE
- [ ] Rhie-Chow
- [ ] Cavity validation
- [ ] Poiseuille validation

## M2 — Turbulence

- [ ] RANS k-epsilon
- [ ] RANS k-omega SST
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

- [ ] Fluid-structure interaction
- [ ] Partitioned / monolithic coupling

## M8 — Reacting flows / combustion

- [ ] Species transport
- [ ] Finite-rate chemistry
- [ ] Flame models

## M9 — Discrete phase model

- [ ] Lagrangian particle tracking
- [ ] One/two-way coupling

## M10 — Optimization / adjoint

- [ ] Discrete adjoint
- [ ] Shape optimization
