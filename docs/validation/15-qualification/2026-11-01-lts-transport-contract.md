# Lot 3A - dimensional local CFL and transactional transport fields (2026-11-01)

Status: unit tests PASS on a **local reconstruction** of PR #1155 lot 1 + lot 2.
This reconstruction is **not** a verified copy of the remote GitHub HEAD.
No CFD physical benchmark was promoted by this patch.

## LTS-01: local CFL is a cell-mass-scaled quantity

Previous formula: `dt = CFL * rho / sum(abs(m_dot_f))`.  This has units
`seconds / cubic metre` and is only correct for unit-volume cells.  Existing
unit tests used a one-cubic-metre mesh, hiding the missing factor.

Corrected formula, with `m_dot` in kg/s, density in kg/m^3 and V in m^3:

`CFL_c = dt_c * sum_f(abs(m_dot_f)) / (rho_c * V_c)`

`dt_c = clamp(CFL_target * rho_c * V_c / sum_f(abs(m_dot_f)), dt_min, dt_max)`.

A zero total convective flux selects explicit `dt_max`, as before.  Nonpositive,
nonfinite or nonrepresentable cell volumes, density, face flux, controls or
computed time steps now reject the update without mutating the previous `dt`.
The global reduction rejects nonfinite and nonpositive local dt values.

The original Mesh-based signature remains available; it computes oriented cell
volumes before dispatching.  For performance-sensitive repeated calls,
**prefer the new `std::vector<double> cell_volumes` overload with cached volumes**.
The old helper was observed only in tests within the provided source snapshot:
this fixes a public numerical utility but does not prove a speedup in the
production COUPLED solver. `growth_limit` remains a validated but otherwise
unused control; no silent change to its historical behavior was introduced.

Test: one cubic metre, eight cubic metres and one eighth cubic metre at the
same speed and density, plus explicit volumes and transactional failures.
Before the fix: `L=2 m` gave `0.0625 s` instead of `0.5 s`.

## THERMO-02: no partial material-field mutation on bad input

The previous `compute_transport_fields` loop wrote viscosity/conductivity/
diffusivity for preceding cells before a later cell raised an exception.
Pressure-field or output-field size mismatches could also fail only after
partial writes.  Array outputs sharing the same field could be overwritten.

The two overloads now validate scalar shape/size, reject aliases, compute all
properties into a single `O(n_cells)` scratch array, and only then write the
three output fields.  This preserves previous outputs on input/physics errors.
The Prandtl and Schmidt property functions also reject unrepresentable values
instead of returning infinite conductivities/diffusivities.

**Performance tradeoff:** one contiguous scratch allocation of three doubles
per cell avoids recomputing the Sutherland law twice; this has not yet been
profiled at large production scale. The functions had no detected production
call site in the inspected source tree and should not be advertised as a fix
for NACA or SST convergence.

## Executable tests

- `test_local_cfl_physical_contract`: 5 cases (volume scaling, cached volume,
  invalid volume, transactional invalid mass flux, global dt rejection).
- `test_transport_field_contract`: 6 cases (late invalid cell, EOS failure,
  shape, aliases, coefficient overflow, healthy oracle).
- Also run `test_local_time_stepping`, `test_transport_models`,
  `test_numerical_model_verification`, plus the previous lot energy/SST tests.

No existing qualification tolerances or physical gates were modified.
