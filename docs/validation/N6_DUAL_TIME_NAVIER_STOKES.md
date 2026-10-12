# N6 production dual-time Navier–Stokes

## Scope

This increment connects the dual-time physical-time contract to the production incompressible finite-volume momentum equations.

The existing nonlinear pressure-velocity loop is used as the pseudo-time solve. The physical derivative is assembled into the momentum matrix, so SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step and coupled pressure-velocity paths keep their existing production flux and conservation machinery.

## Variable-step BDF2

For accepted history U^n,U^(n-1) and omega = dt_n/dt_(n-1), CFDX uses the verified variable-step BDF2 coefficients from cfdx/core/numerics/temporal.h:

  a0=(1+omega)^2, a1=omega^2, b=(1+omega) dt_n, denominator=1+2 omega.

The production momentum matrix receives M_P = rho V (1+2 omega)/b and the temporal history source receives b_P = rho V (a0 U^n - a1 U^(n-1))/b.

The first accepted step uses backward Euler and creates the BDF2 history.

## Physical-step acceptance

A physical step is accepted only if:

1. the production nonlinear pressure-velocity solve converges;
2. the independent velocity temporal-error estimate is <= 1;
3. U and p are finite/admissible.

The embedded temporal estimator compares the production BDF2 solution with a backward-Euler solution from the same accepted state.

## Retry transaction

A rejected step restores U, p, previous accepted velocity, history-valid flag, previous accepted dt, physical time and accepted step count.

The retry reduces dt by the configured bounded shrink factor. No rejected attempt advances physical time or accepted history.

## Restart equivalence

DualTimeNavierStokesCheckpoint contains the complete temporal state required by the lifecycle. The executable N6 test restores an accepted checkpoint into two independent production solver lifecycles and requires machine-level agreement of U, p and physical time after the same continuation step.

## V&V boundary

The temporal-order oracle remains the production-FVM MMS campaign from N5 (#552/#573). This N6 campaign verifies the missing production pressure-velocity-coupled integration boundary, variable-step history use, physical retry transaction and restart equivalence.

This does not yet claim compressible acoustic CFL, MPI/GPU numerical equivalence, or full external-flow transient qualification.

No tolerance relaxation or validation disabling is used.