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
- [x] View-factor matrix validation
- [x] P1 source-term primitive
- [x] DOM quadrature validation primitive
- [x] Full participating-media transport solve
- [x] Radiation/energy coupling validation

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
