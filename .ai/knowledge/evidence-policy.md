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
## Claim discipline

Agents must distinguish observed facts from interpretation and unresolved risk. They must never fill an evidence gap with model memory or assumption. When evidence is unavailable, contradictory or produced by a failed tool, the report must say so explicitly.

## Status vocabulary

Use precise status terms such as `implemented`, `tested`, `verified`, `validated`, `qualified`, `partial`, `experimental`, `optional`, `blocked` and `unknown` according to the evidence actually available. Do not use `pass` as a synonym for qualification.

## Negative evidence

Failed tests, non-convergent refinements, skipped checks, unsupported paths and known limitations are first-class evidence and must be preserved until explicitly resolved.
