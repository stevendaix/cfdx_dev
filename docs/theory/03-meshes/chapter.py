# %% [markdown]
"""# Meshes

## Scientific purpose

This chapter is written as executable literate documentation. It must teach the mathematical chain rather than merely list software features.

## Required structure

1. Physical motivation and problem definition.
2. Variables, dimensions, units and sign conventions.
3. Governing equations and assumptions.
4. Control-volume formulation.
5. Discrete formulation and algebraic consequences.
6. Consistency, conservation, stability, boundedness and accuracy.
7. CFDX implementation mapping.
8. Executable verification and benchmark evidence.
9. Limitations and improvement paths.
10. References.

The status vocabulary is strict: Implemented, Verified, Validated, Qualified. A code path is not evidence of numerical correctness.

## Core equations

- S_f=n_f A_f
- V_P=(1/3) sum_f A_f r_f.n_f for a closed polyhedron
- sum_f S_f=0
- d_f=|C_N-C_P|
- non-orthogonality angle = acos((d_f.S_f)/(|d_f||S_f|))
- skewness measures offset between face centroid and centre-to-centre intersection

## Scientific checks

The executable cells below are intentionally deterministic. They check mathematical identities or toy/reference systems only. They must not be presented as CFDX solver qualification unless the corresponding CFDX evidence is explicitly linked.

## CFDX traceability

The final chapter must identify the exact source files, tests and V&V artifacts corresponding to each equation. If an implementation is partial, the documentation must say so explicitly.

## References

Use the central BibTeX bibliography. Add stable identifiers (DOI, publisher, standard or project URL) to the authoritative record rather than duplicating metadata here.
"""

# %%
from __future__ import annotations
import numpy as np

# Deterministic mathematical sanity checks
x = np.linspace(0.0, 1.0, 5)
assert np.all(np.isfinite(x))
Sf=np.array([[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.],[0.,0.,1.],[0.,0.,-1.]])
assert np.allclose(Sf.sum(axis=0),0.0)
