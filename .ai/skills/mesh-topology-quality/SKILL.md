# Mesh, topology and quality

## Purpose

Audit and reason about mesh topology, geometry, metrics, and quality as numerical inputs to CFDX.

## Required reasoning

1. Establish cell/face/vertex ownership and adjacency from the repository data model.
2. Check face orientation and outward-normal conventions.
3. Check cell volumes, face areas, centroids, and geometric consistency.
4. Identify non-orthogonality, skewness, aspect ratio, degeneracy, and boundary-layer resolution issues when relevant.
5. Distinguish topological validity from numerical quality.
6. Trace imported/generated mesh data to the operators that consume it.

## Numerical contract

Geometric identities used by a discretisation must be checked independently where possible. A mesh that parses successfully is not necessarily suitable for a requested discretisation or convergence study.

## Verification

Use simple analytical geometries, affine mappings, manufactured geometry checks, mesh-quality diagnostics, and refinement families. Record the mesh definition and quality metrics with benchmark evidence.

## Anti-patterns

- Do not treat parser success as mesh qualification.
- Do not invent mesh-quality thresholds without identifying the affected numerical method.
- Do not hide degenerate cells or faces by filtering them without reporting the change.
