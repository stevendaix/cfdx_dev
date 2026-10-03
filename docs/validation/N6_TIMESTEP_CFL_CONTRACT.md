# N6 — Time-step, CFL and pseudo-time control contract

## Scope

N6 closes the time-step control contract defined by issue #461. It covers deterministic CFL/time-step decisions, geometry-aware characteristic length, auditable change reasons, and transactional rollback of a failed nonlinear step.

N6 does not claim compressible wave-speed CFL, dual-time stepping, or a production nonlinear retry loop where the solver itself has not yet exposed that integration boundary. Those remain outside the current qualification claim.

## CFL definition

For a cell with volume V_i, volumetric face fluxes Phi_f, and physical time step Delta t:

CFL_i = Delta t / V_i * sum_f |Phi_f|.

The reported global convective CFL is max_i(CFL_i), and the limiting cell is reported explicitly.

The implementation requires finite fluxes and positive finite cell volumes. Invalid geometry or data is rejected rather than converted to a fallback value.

## Characteristic length

For geometry-driven time-step estimates, CFDX uses:

h_i = 2 V_i / A_i,

where A_i is the sum of the areas of the cell faces.

This hydraulic-diameter-style geometric scale is deliberately sensitive to thin/high-aspect-ratio cells and avoids the isotropic V^(1/3) assumption. It is a control/diagnostic scale and does not replace the face-flux CFL calculation.

## Deterministic controller

The adaptive controller proposes:

Delta t_new = clip(
    Delta t * sqrt(CFL_target / CFL_measured),
    Delta t * shrink_limit,
    Delta t * growth_limit
),

followed by the absolute [Delta t_min, Delta t_max] bounds.

Every recorded decision contains the old timestep, new timestep, measured CFL and explicit reason. The controller has no random or hidden state-dependent policy.

## Rollback contract

A rejected nonlinear step must not mutate the physical state or multistep history before retry.

TimeStepRollback snapshots and restores:

- physical cell state;
- phi_prev;
- phi_curr;
- has_prev;
- dt_prev.

The deterministic retry controller reduces the timestep by shrink_limit and enforces dt_min and max_retries.

## Verification

tests/unit/test_timestep_control.cpp verifies:

1. characteristic length reacts to a thin/high-aspect-ratio cell;
2. local and global CFL are computed from the documented flux/volume definition;
3. identical inputs produce identical timestep decisions;
4. timestep bounds and deterministic rollback decisions are enforced;
5. rollback restores the physical state and complete temporal history;
6. rollback cannot occur without an active transaction.

## Qualification boundary

The N6 acceptance claim is limited to the implemented deterministic control and diagnostic layer. Compressible acoustic CFL, dual-time stepping and full production solver retry integration remain explicitly unclaimed until their production contracts exist.
