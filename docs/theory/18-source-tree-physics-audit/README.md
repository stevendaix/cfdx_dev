> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 18 — CFDX Source-Tree Physics and Numerical Audit

**Status: AUDIT REGISTER — generated from the current PR #530 repository tree.**

This chapter is the index for the file-by-file Theory audit. A row records the source family and its scientific responsibility; it is not a claim that the implementation is qualified.

## Mesh

| File | Scientific role |
|---|---|
| core/mesh/point.* | point coordinates |
| core/mesh/face.* | polygon connectivity |
| core/mesh/cell.* | cell-to-face connectivity |
| core/mesh/ownership.* | owner/neighbour orientation |
| core/mesh/boundary.* | boundary patch topology |
| core/mesh/mesh.* | aggregate topology and validation |
| core/mesh/index_types.h | index/offset types |
| core/mesh/solver_mesh.h | solver-facing mesh view |
| core/mesh/sfc_ordering.h | space-filling-curve ordering |

## Geometry

| File | Scientific role |
|---|---|
| core/geometry/face_geometry.* | face centre, area and oriented surface vector |
| core/geometry/cell_geometry.* | cell centre and signed/absolute volume |
| core/geometry/geometry_cache.h | reusable derived geometry |
| core/geometry/mesh_quality.* | skewness/non-orthogonality metrics |
| core/geometry/mesh_validator.* | topology/geometry/quality/conservation checks |

## Fields

| File | Scientific role |
|---|---|
| core/field/field.* | typed field values, location and metadata |
| core/field/storage.* | host/device storage state and ownership |

The field contract must always state location, dimension, units and precision. A cell field and a face field are different numerical objects.

## Boundary

| File | Scientific role |
|---|---|
| core/boundary/boundary.h | boundary aggregate |
| boundary_condition.h | condition interface |
| boundary_constraint.h | mathematical constraint representation |
| boundary_field.h | boundary field values |
| boundary_role.h | physical/geometric role |
| boundary_validator.h | consistency checks |
| flow_boundary_conditions.h | flow-specific conditions |
| legacy_patch_field_adapter.h | migration compatibility |
| mathematical_condition.h | mathematical condition |
| value_provider.h | imposed/evaluated values |

The boundary audit must trace each condition from physical statement to matrix/RHS contribution.

## FVM and numerics

| File | Scientific role |
|---|---|
| core/fvm/least_squares_gradient.h | LS/WLS gradient assembly |
| core/numerics/gradient.h | gradient dispatch and implementations |
| gradient_stencil.h | stencil/boundary contract |
| interpolation.h | face interpolation |
| convection.h | convective discretisation |
| convection_assembly.* | convection matrix/flux assembly |
| divergence.h | divergence operator |
| laplacian.h | diffusion/Laplacian operator |
| flux.h | flux abstraction |
| source_term.h | source discretisation |
| temporal.h | temporal discretisation |
| conservation.h | conservation diagnostics/contracts |
| numerical_method_contract.h | common numerical contract |
| numerical_method_registry.h | method capability registry |
| numerical_method_selection.h | case-level method selection |
| matrix_free.h | matrix-free operator path |
| newton.h | nonlinear linearisation support |
| integrate.h | numerical integration |

The Theory audit must keep gradient, interpolation, flux, diffusion, source and time integration separate.

## Linear algebra

The current tree contains sparse_matrix/vector/linear_system primitives, CG/BiCGStab/GMRES and stationary solvers, Krylov controls/reductions, linear operator and dispatch layers, AMG, Chebyshev smoothing, MGR, block/Schur/field-split infrastructure, null-space handling, matrix-free operators, mixed precision and MPI overlap.

Representative files are:

- core/linalg/linear_system.h
- linear_solver_dispatch.h
- linear_solver_models.h
- cg_solver.h
- bicgstab_solver.h
- gmres_solver.h
- amg_preconditioner.h
- amg_wrapper.h
- mgr_preconditioner.h
- exact_schur.h
- block_schur.h
- coupled_amg_schur.h
- native_fieldsplit.h
- null_space.h
- matrix_free_fv_operator.h
- mixed_precision.h

The audit must distinguish solver convergence, true residual, energy contraction and physical convergence.

## Physics

The current physics family includes:

- incompressible.h
- finite_volume_transport.h
- steady_incompressible_solver.h
- pressure_velocity.h and pressure_velocity_algorithms.h
- thermal.h and energy_solver.h
- cht_solver.h
- radiation*.h
- turbulence*.h
- spalart_allmaras.h
- sst_solver.h
- wall_distance*.h
- adaptive_cfl.h
- local_time_stepping.h
- low_storage_time_integration.h
- conservation_boundedness.h
- axisymmetric_fvm.h and axisymmetric_vmfl036_solver.h
- low_mach.h and compressible_flux.h
- boussinesq.h
- fsi.h
- vof.h
- multi_component.h
- poisson/poisson.h

Each model requires its own governing equation, assumptions, source terms, numerical coupling, verification and validation population.

## I/O and data model

The current HDF5 family is:

- io/hdf5/schema.h
- case_hdf5_io.h/.cpp
- hdf5_reader.h/.cpp
- hdf5_writer.h/.cpp
- mini_json.h

Mesh import/output families include:

- io/mesh/mesh_importer.*
- io/openfoam/openfoam_importer.*
- io/gmsh/gmsh_importer.*
- io/vtu/vtu_writer.*

These are data-model transformations and must document mapping and information loss.

## Audit completion rule

A file is not considered fully audited until its row can be expanded with exact tests, benchmark evidence, mathematical assumptions, limitations and references. Chapter 18 is therefore a living register, not a fabricated completion claim.

The immediate next audit order is:

    mesh → geometry → fields → boundary → FVM/numerics
    → physics → linear algebra → solver → I/O

This order follows the computational dependency graph.


## Scientific explanation standard

For every important equation, algorithm or data-model rule, document the motivation; definitions/units/sign conventions; derivation or formal rationale; discrete/FVM representation where applicable; actual CFDX execution path; example; error/limitation/sensitivity analysis; exact implementation files; executable verification; benchmark/V&V evidence; and stable bibliography/reference identifiers. Tables summarize explanations and do not replace them. Keep **Implemented / Verified / Validated / Qualified** distinct; documentation maturity never promotes numerical maturity.