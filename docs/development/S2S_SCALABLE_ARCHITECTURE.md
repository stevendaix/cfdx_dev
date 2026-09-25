# S2S scalable radiation architecture

## Why

The current S2S prototype contains two complementary estimators:

- a dense centroid view-factor matrix, O(N^2) memory;
- deterministic cosine-weighted ray tracing, whose visibility stage previously scanned every target and blocker triangle for every ray.

The second path is the important scalable direction, but it must not retain an O(N_triangles) traversal cost.

## Phase 1 implemented by this PR

- CPU BVH/AABB hierarchy for target and blocker surfaces.
- Ray/AABB traversal with nearest-hit queries.
- The existing deterministic ray sampler and public API remain unchanged.
- A regression test compares BVH nearest-hit results with the original brute-force triangle intersection.
- The BVH stores triangle indices rather than duplicating triangle geometry.

This changes the ray visibility search from a full triangle scan to spatial traversal while preserving the existing numerical estimator.

## Phase 2: view-factor storage

The next step is to make the radiation operator independent from its storage:

1. Dense storage for small surface counts.
2. CSR/sparse storage with a configurable cutoff for large but localized interactions.
3. Reciprocity diagnostics using A_i F_ij = A_j F_ji.
4. Enclosure closure diagnostics and controlled correction.
5. Memory estimates exposed before allocation.

The dense matrix must remain an explicit small-case option, not the default for industrial-scale surfaces.

## Phase 3: matrix-free S2S

For very large surfaces, avoid storing F entirely:

    G = F J

becomes an operator application implemented by ray tracing. The surface/BVH data stays resident and only the required radiation fields are exchanged.

The API should therefore expose an operator abstraction rather than forcing callers to depend on a dense matrix.

## Phase 4: CPU/GPU execution

The CUDA implementation should use the same surface representation and BVH concepts:

- host/device resident triangle SoA;
- GPU-resident BVH;
- batched ray generation/traversal;
- device-side accumulation;
- explicit Working Set accounting;
- no silent CPU fallback.

GPU availability/equivalence remains a hardware-gated validation item.

## Validation gates

Every implementation phase must retain:

- analytical enclosure checks;
- reciprocity;
- closure;
- ray-count convergence;
- CPU brute-force/BVH equivalence;
- dense/sparse/matrix-free operator equivalence;
- memory scaling;
- CPU/GPU equivalence when CUDA hardware is available.

The existing S2S V&V remains unchanged by Phase 1.
