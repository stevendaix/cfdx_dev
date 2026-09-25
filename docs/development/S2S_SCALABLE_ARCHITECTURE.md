# S2S scalable radiation architecture

## Current limitations

The dense centroid view-factor matrix requires O(N²) memory. The deterministic ray tracer also previously scanned every target and blocker triangle for every ray.

## Phase 1 implemented by this PR

- CPU AABB/BVH hierarchy for target and blocker surfaces.
- Ray/AABB traversal with nearest-hit queries.
- Existing deterministic ray sampling and public API remain unchanged.
- Regression coverage compares BVH nearest-hit results with brute-force triangle intersection.
- BVH stores triangle indices and does not duplicate triangle geometry. The
  referenced triangle vector must outlive the BVH and remain unchanged while
  it is queried; rvalue construction/building is rejected to prevent an
  immediately dangling geometry reference.

This removes the full triangle scan from the visibility path without changing the numerical estimator.

## Next phases

1. Separate the radiation operator from its storage: dense for small systems, CSR for sparse interactions, and matrix-free for very large systems.
2. Add reciprocity A_i F_ij = A_j F_ji, enclosure closure diagnostics, controlled correction, and memory estimates.
3. Add CPU parallel execution and memory-aware working sets.
4. Add CUDA-resident triangle SoA/BVH and batched ray traversal, with explicit hardware gating and no silent CPU fallback.
5. Validate dense/sparse/matrix-free equivalence, memory scaling, ray-count convergence, and CPU/GPU equivalence when CUDA hardware is available.

The dense matrix remains a small-case option; it should not be the industrial-scale default.
