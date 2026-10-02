# %% [markdown]
"""# Gradients and Reconstruction

## Scientific objective

Green-Gauss follows from the divergence theorem: integral_V grad(phi)=integral_boundary phi*n, hence grad(phi)_P is approximately (1/V_P) sum_f phi_f*S_f. Least squares starts from phi_N-phi_P=g.d_N+O(h^2) and minimises sum w_N*(delta_phi_N-g.d_N)^2, giving A*g=b with A=sum w*d*d^T and b=sum w*d*delta_phi. Rank deficiency must be diagnosed.

## Executable contract

This deterministic cell is a minimal mathematical sanity check. Quantitative studies must be linked to V&V evidence.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

import numpy as np; d=np.array([[1.,0.],[0.,1.],[-1.,0.]]); g=np.array([2.,-3.]); A=d.T@d; b=d.T@(d@g); assert np.allclose(np.linalg.solve(A,b),g)
