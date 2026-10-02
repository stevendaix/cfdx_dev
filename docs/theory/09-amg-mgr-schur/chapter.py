# %% [markdown]
"""# AMG, MGR and Schur

Scientific executable chapter. Follow physical motivation, definitions, assumptions, governing equations, discretisation, numerical properties, CFDX implementation, verification evidence, limitations and references. Implementation is not automatically verification, validation or qualification.

## Core equations and contracts

- A=[A_uu G;D C]
- S=C-D A_uu^-1 G
- S_tilde=C-D M_b^-1 G
- A_H=R A_h P
- V-cycle = pre-smooth + restrict + coarse solve + prolongate + post-smooth
- energy contraction q=||e_new||_A/||e_old||_A

## Evidence rule

Every implementation claim must map to actual source files and tests. Every numerical result must identify configuration, mesh/resolution, solver settings and software revision. Missing evidence is a documented gap.
"""

# %%
from __future__ import annotations
import numpy as np
x=np.linspace(0.0,1.0,5)
assert np.all(np.isfinite(x))
A=np.array([[2.,1.],[1.,3.]])
b=np.array([1.,2.])
x=np.linalg.solve(A,b)
assert np.allclose(A@x,b)
