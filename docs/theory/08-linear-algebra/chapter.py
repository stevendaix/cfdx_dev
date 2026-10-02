# %% [markdown]
"""# Linear Algebra

## Scientific objective

The discretised problem is A*x=b. CG applies to symmetric positive-definite systems; GMRES treats general nonsymmetric systems; BiCGStab provides a short-recurrence nonsymmetric alternative. Preconditioning changes the algebraic problem to improve convergence. Linear residual reduction is not itself a physical accuracy metric.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

import numpy as np; A=np.array([[4.,1.],[1.,3.]]); b=np.array([1.,2.]); x=np.linalg.solve(A,b); assert np.linalg.norm(b-A@x)<1e-14
