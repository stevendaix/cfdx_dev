---
name: cfdx-benchmark-reference
description: Execute CFDX benchmark and independent-reference comparisons without conflating regression and validation.
version: 1
---

# Principles

- Identify the published/independent reference and its assumptions.
- Record geometry, boundary conditions, physical parameters, mesh and numerical settings.
- Distinguish same-formulation, analogous and independent implementations.
- Compare predefined quantities of interest rather than visually selecting favorable results.
- Report differences, uncertainty and known reference limitations.
- Do not treat agreement with one benchmark as general validation.

# Evidence

A benchmark comparison is validation evidence only for the scope it actually exercises. Internal regression tests remain distinct from independent reference comparisons.
