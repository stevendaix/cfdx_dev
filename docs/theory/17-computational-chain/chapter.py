# %% [markdown]
"""# Computational Chain and Code Map

Scientific executable chapter. Follow physical motivation, definitions, assumptions, governing equations, discretisation, numerical properties, CFDX implementation, verification evidence, limitations and references. Implementation is not automatically verification, validation or qualification.

## Core equations and contracts

- case -> schema -> topology -> geometry -> fields -> BC -> numerics -> physics -> matrix -> solver -> diagnostics -> output
- each arrow is a contract boundary
- A x=b is produced by discretisation and consumed by the linear solver
- diagnostics are independently recomputed where required
- V&V evidence closes the chain back to requirements

## Evidence rule

Every implementation claim must map to actual source files and tests. Every numerical result must identify configuration, mesh/resolution, solver settings and software revision. Missing evidence is a documented gap.
"""

# %%
from __future__ import annotations
import numpy as np
x=np.linspace(0.0,1.0,5)
assert np.all(np.isfinite(x))
stages=["case","mesh","fields","numerics","physics","matrix","solver","diagnostics","output"]
assert len(stages)==9
