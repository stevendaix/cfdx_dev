# N9-S10 — pressure-reference policies

N9-S10 makes the pressure level an explicit solver policy for pure-Neumann incompressible cases.

## Policies

- `REFERENCE_CELL`: apply a uniform pressure shift so the configured reference cell equals `pressure_reference_value`.
- `ZERO_MEAN`: apply a uniform shift so the volume-weighted domain mean pressure is zero.
- `NONE`: do not apply a pressure-level normalization. For a pure-Neumann segregated solve, a compatible constant null-space projection must be explicitly requested.

A fixed-pressure boundary condition always takes precedence because it defines a physical pressure level; the gauge policy must not overwrite that physical boundary condition.

The pressure gradient, velocity, conservative flux, continuity, and momentum solution are unchanged by the uniform gauge shift. The policy therefore controls only the pressure level, not the pressure-gradient physics.

## Qualification

`test_n9_s10_pressure_reference_policies` exercises all three policies on the production incompressible solver and checks:

- configured reference value;
- zero volume-weighted mean while preserving pressure differences;
- preservation of an existing pressure level when no gauge normalization is requested.

The CTest entry is labelled `n9;qualification;validation;long;validation-total`.

## Acceptance

N9-S10 is not a claim that N9 is closed. N9 remains partial until S7, S8, S9, S10 and the final matrix/evidence audit are complete.

No tolerance relaxation, disabled test, silent fallback, or gauge-specific empirical correction is permitted.
