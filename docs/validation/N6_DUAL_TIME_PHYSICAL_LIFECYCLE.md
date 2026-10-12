# N6 physical-time retry and restart lifecycle

This increment builds on the field-level dual-time contract and adds the missing physical-time transaction boundary.

## Contract

For each requested physical step:

1. clamp the requested physical time step to configured bounds;
2. invoke the transient step callback with the accepted state and trial time step;
3. run the existing dual-time field integrator;
4. accept only when its temporal, nonlinear and admissibility gates pass;
5. advance physical time and accepted-step index only after acceptance;
6. on a rejected dual-time step, shrink the time step deterministically and retry;
7. never mutate accepted temporal history as part of a rejected attempt.

The pseudo-time nonlinear solve remains owned by DualTimeFieldIntegrator; this layer does not introduce a second pseudo-time algorithm.

## Restart contract

DualTimePhysicalCheckpoint captures:

- accepted current field;
- accepted previous field;
- whether BDF2 history is valid;
- previous accepted physical time step;
- physical time;
- accepted step index.

Restoring this object reseeds the dual-time history without changing the numerical scheme.

## Scope boundary

This increment closes the physical-time transaction/retry and restart lifecycle around the existing field-level contract. It does not claim full Navier–Stokes production dual-time qualification.

Still required for full N6 production qualification:

- coupling the dual-time residual to the production Navier–Stokes/FVM residual;
- pressure-velocity coupling inside the dual-time nonlinear solve;
- production transient boundary-condition adapters;
- end-to-end transient MMS temporal-order evidence;
- restart/retry equivalence on the production solver;
- compressible wave-speed CFL when compressible production support exists.

No tolerance is relaxed and no validation is disabled.
