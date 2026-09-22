# CFDX — Issue #34 evidence matrix

Updated: 2026-09-22

This document is the working evidence map for issue #34. A checked implementation item is not considered validated until the required numerical/physical evidence and CI gate are present.

## Evidence levels

| Level | Meaning | Acceptance evidence |
|---|---|---|
| U | Unit | deterministic API/invariant test |
| A | Analytical | exact manufactured/analytical oracle |
| M | Manufactured | measured L2/Linf error against MMS source |
| S | Solver | end-to-end solver residual + physical invariant |
| R | Reference | independent OpenFOAM/published/reference comparison |
| C | Convergence | observed order from at least 3 useful resolutions where applicable |
| G | Regression | prior accepted result remains unchanged |
| CI | Integration | clean-checkout required jobs pass |

## Baseline snapshot

- Target: master
- Baseline SHA: c5bd51afa1cde629800e4526e131ec42dfe70593
- Issue: #34
- Rule: no roadmap item is closed from code presence alone.

## M0 evidence map

| Scope | Code | Existing tests/evidence | Gap tracked by #34 |
|---|---|---|---|
| HDF5 case/fields | src/cfdx/io/hdf5/* | round-trip/integrity tests | full case reconstruction + corruption matrix |
| OpenFOAM import | src/cfdx/io/openfoam/* | regression geometry | cross-import equivalence + parser robustness |
| Gmsh/meshio | src/cfdx/io/gmsh/*, Python bridge | geometry matrix | independent equivalence matrix |
| Geometry cache | src/cfdx/core/geometry/geometry_cache.h | orientation/quality regressions | all mesh families + malformed topology |
| Gradient | src/cfdx/core/numerics/gradient.h | unit/constant/linear foundations | quadratic + MMS + convergence |
| Divergence | src/cfdx/core/numerics/divergence.h | conservation foundation | analytical/MMS + skew/non-orthogonal |
| Laplacian | src/cfdx/core/numerics/laplacian.h | constant-field foundation | analytical quadratic + non-orthogonal + convergence |
| Interpolation | src/cfdx/core/numerics/interpolation.h | limiter bounds | face-value analytical/MMS verification |
| Matrix-free | src/cfdx/core/numerics/matrix_free.h, FVM operator | equivalence tests exist in part | complete assembled↔matrix-free verification |
| Linear algebra | src/cfdx/core/linalg/* | solver foundations | robustness/reference matrix |
| MPI | src/cfdx/core/parallel/* | multi-rank foundation | serial/MPI equivalence + restart |
| Runtime/GPU | execution/memory policies | foundations | full GPU/OOC equivalence remains later phases |

## Current blockers

1. Full operator verification is not yet complete; implementation presence is insufficient.
2. Non-orthogonal/skewed discretisation needs quantitative evidence.
3. Matrix-free equivalence must be demonstrated against the assembled operator, not only API presence.
4. Phase completion must remain synchronized with issue #34 rather than legacy roadmap checkboxes.

## Merge gate

For each PR: implementation → unit test → analytical/MMS/physical validation where applicable → quantitative result → CI → code/spec audit → merge → roadmap update.
