# Skill: Temporal Discretisation

## Purpose
Guide implementation and verification of physical-time and pseudo-time discretisation.

## Scope
- Explicit and implicit time integration.
- First/second-order schemes and multistep/startup treatment.
- Stability restrictions and CFL-like conditions.
- Physical time versus pseudo-time.
- Checkpoint/restart consistency.
- Temporal refinement and Richardson/MMS verification.

## Method
1. Define the semi-discrete equation and time operator mathematically.
2. State the formal temporal order and startup order.
3. Identify stability restrictions and nonlinear/linear solve tolerances separately.
4. Keep pseudo-time acceleration distinct from physical-time accuracy.
5. Verify restart by comparing uninterrupted and restarted trajectories.
6. Use temporal refinement with spatial error controlled or analytically removed.

## Required evidence
- Manufactured solution or modal/periodic analytical test.
- Temporal refinement with measured slope.
- Stability test at documented time-step/CFL limits.
- Restart equivalence test.
- Explicit reporting of solver tolerances so iteration error is not mistaken for temporal error.

## Review rules
A converged nonlinear solve does not prove temporal order. Do not change tolerances, startup treatment, or time-step definitions merely to obtain a desired slope.
