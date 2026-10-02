# %% [markdown]
"""# Multiphysics

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

- F(y)=0
- J(y) delta y=-F(y)
- J=[F_u,u F_u,p;F_p,u F_p,p]
- segregated iteration y^(k+1)=S(y^k)
- monolithic solve J delta y=-F
- interface residuals enforce conservation between subdomains

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
