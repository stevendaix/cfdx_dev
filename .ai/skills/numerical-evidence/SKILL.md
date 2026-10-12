# Numerical evidence

## Purpose

Provide a common evidence discipline for CFDX numerical development.

## Evidence ladder

Keep these claims separate:

- **Implemented:** code path exists and is reachable.
- **Tested:** automated checks exercise it.
- **Verified:** numerical properties are demonstrated against an exact/analytical/discrete reference.
- **Validated:** results agree with an appropriate physical/reference benchmark within stated uncertainty.
- **Qualified:** the acceptance criteria defined by the CFDX roadmap are met.

## Required report fields

For a numerical result record:

- method and configuration;
- mesh/refinement family;
- boundary and initial conditions;
- solver/numerical settings;
- metric and norm;
- reference solution;
- observed error/order or conservation/boundedness measure;
- tolerance and stopping criteria;
- commit SHA/PR;
- limitations and known failure modes.

## Gate discipline

A green CI run proves only the configured checks passed. It does not automatically prove numerical qualification.

Never obtain a green result by disabling a required test, inflating a tolerance, weakening an acceptance criterion, or silently changing the problem definition.
