# %% [markdown]
"""# Gradients and Reconstruction

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

- int_V grad(phi)dV = int_boundary phi n dA
- grad(phi)_P approx (1/V_P) sum_f phi_f S_f
- phi_N-phi_P = grad(phi)_P.d_N + O(h^2)
- A=sum_N w_N d_N d_N^T; b=sum_N w_N d_N Delta_phi_N; A g=b
- rank(A)<d implies an underdetermined gradient stencil
- kappa(A) measures sensitivity of the least-squares solve
- phi_f=(1-w)phi_P+w phi_N

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
d = np.array([[1.,0.,0.],[0.,1.,0.],[-1.,0.,0.],[0.,-1.,0.]])
g = np.array([2.,-3.,1.])
A = d.T @ d
b = d.T @ (d @ g)
assert np.allclose(np.linalg.lstsq(A, b, rcond=None)[0], g)
