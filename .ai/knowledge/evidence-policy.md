# CFDX AI Evidence Policy

AI-generated conclusions must identify the evidence supporting them.

## Evidence categories

- implementation: code exists;
- test: automated test exercises behavior;
- verification: quantitative comparison against a known/independent result;
- validation: comparison demonstrating suitability for the intended physical/numerical scope;
- qualification: acceptance against the project's explicit maturity criteria.

These categories are not interchangeable.

## Numerical changes

The agent should seek mathematical formulation, discrete formulation, implementation, regression tests and independent V&V evidence as appropriate.

## CI

A green CI result is evidence that the configured CI checks passed. It is not by itself proof of numerical qualification.

## Missing evidence

Missing or unavailable evidence must be reported explicitly rather than inferred.