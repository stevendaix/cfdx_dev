# CFDX Boundary Conditions Architecture

## Status

Target architecture for issue #446.

This document defines the contract before solver migration. The existing `PatchField` remains a numerical compatibility layer during migration.

## Design rule

CFDX separates five conceptual layers:

1. **Boundary** — patch identity, faces and geometry/topology.
2. **BoundaryCondition** — physical intent.
3. **ValueProvider** — data/profile generation.
4. **MathematicalCondition** — mathematical constraint.
5. **Discretization** — FVM realization.

Execution is a separate concern: the same mathematical constraint must be lowerable to CPU, GPU and MPI execution without changing its physical meaning.

## Boundary

A boundary identifies a set of mesh faces. It may carry a semantic `BoundaryRole`:

- inlet
- outlet
- wall
- symmetry
- periodic
- interface

The role is metadata and must never be used to infer a field-specific physical condition.

A patch named `inlet` is therefore not implicitly a `VelocityInlet`.

## Physical boundary condition

A `BoundaryCondition` expresses physical intent and can produce constraints for multiple fields.

Examples:

- `VelocityInlet`
- `PressureOutlet`
- `NoSlip`
- thermal wall
- radiation wall
- periodic coupling

The physical layer must not depend on SIMPLE, SIMPLEC, PISO, PIMPLE or a particular linear-system implementation.

The preferred flow is:

```text
BoundaryCondition
      |
      v
BoundaryConstraintSet
      |
      +---- U -> MathematicalCondition
      +---- p -> MathematicalCondition
      +---- T -> MathematicalCondition
      |
      v
FVM discretization
```

## Mathematical conditions

The initial mathematical vocabulary is:

- Dirichlet
- Neumann
- Robin
- Flux
- Mixed / flux-dependent
- Coupled
- Periodic

These are solver-independent contracts. Numerical assembly decides how a condition contributes to a particular discretized equation.

Pressure-outlet backflow is explicitly a mixed/flux-dependent condition:

```text
outflow  -> Neumann/extrapolation
backflow -> prescribed backflow state
```

The switching and linearisation policy must be defined by the physical condition and discretization contract, not hidden inside `PatchField`.

## Value providers

A value provider represents data, not physical semantics.

The target configuration-side vocabulary includes:

- constant
- time-dependent
- spatial profile
- table/CSV
- expression
- mapped/interpolated field
- coupled value

A provider may be rich and dynamic on the configuration side, but kernels must not require virtual dispatch.

The target execution flow is:

```text
configuration provider
        |
        v
validation / lowering
        |
        v
compact execution representation
        |
        +---- CPU
        +---- GPU
        +---- MPI
```

The exact GPU representation is deliberately left to implementation and benchmarking. No `virtual evaluate()` API is allowed in a performance-critical GPU kernel.

## Multiphysics composition

A boundary is compositional. Flow, thermal, turbulence, radiation and future species/interface conditions are independent components.

Do not create a class hierarchy such as:

```text
Wall
ThermalWall
RadiationWall
ThermalRadiationWall
...
```

A heated radiating wall is instead a boundary with flow + thermal + radiation constraints.

## Well-posedness

The case validator operates on the generated mathematical constraints.

Each solver declares requirements such as:

- required fields;
- acceptable constraint classes;
- global constraints;
- reference/gauge requirements;
- interface compatibility.

Examples of diagnostics include:

- missing pressure reference for a pure-Neumann incompressible problem;
- incompatible periodic pairing;
- missing thermal constraint;
- incompatible coupled-interface constraints.

The validator reports configuration errors before matrix assembly whenever the required information is available.

## Dimensions

Boundary values are dimensioned and checked against the target field before numerical execution.

Persisted dimensions use a structured exponent vector rather than arbitrary unit strings in hot numerical paths.

## Persistence

Boundary-condition configuration belongs to the authoritative case:

```text
.cfdx.h5 = case/setup source of truth
.dat.h5  = numerical state/checkpoint
VTU      = visualization
```

The BC schema is versioned under:

```text
/boundary_conditions/
    schema_version
```

BC configuration is not restart-only runtime state.

## Legacy migration

The old `PatchField` API remains available while the new contract is introduced.

A `LegacyPatchFieldAdapter` will preserve mathematical semantics. It must not infer physical BCs from patch names.

For example:

```text
U @ inlet + fixedValue
        |
        v
Dirichlet(U)
```

not:

```text
patch name = inlet
        |
        v
VelocityInlet
```

This rule prevents accidental semantic changes during migration.

## Solver independence

The same physical BC definition must be usable by SIMPLE, SIMPLEC, PISO, PIMPLE, COUPLED and future formulations.

A solver consumes `BoundaryConstraintSet` and chooses the appropriate discretization/assembly path.

## Initial C++ contract

The implementation is expected to evolve around:

```text
src/cfdx/core/boundary/
    boundary.h
    boundary_role.h
    boundary_condition.h
    boundary_constraint.h
    mathematical_condition.h
    value_provider.h
    boundary_validator.h
    legacy_patch_field_adapter.h
```

The first implementation should remain lightweight and header-oriented where practical. It must not pull solver, CUDA or GUI dependencies into the core boundary contract.

## Validation gates

Architecture acceptance is not numerical qualification.

The implementation must retain existing quantitative tests, especially:

- Poisson mixed/Neumann regression;
- Couette;
- Poiseuille.

No tolerance inflation, disabled test, or silent backend fallback is acceptable.

See #118 and #442 for the engineering qualification contract.
