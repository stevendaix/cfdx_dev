# %% [markdown]
"""# Pressure-Velocity Coupling

## Scientific objective

In incompressible flow, momentum determines velocity while continuity requires div(u)=0. Segregated algorithms construct pressure corrections from linearised momentum relations; SIMPLE uses under-relaxed corrections, while PISO/PIMPLE add correction stages. Coupled methods expose the saddle-point block structure. Verification requires continuity, momentum, mass imbalance and relaxation sensitivity.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

import numpy as np; A=np.array([[4.,1.],[1.,3.]]); b=np.array([1.,2.]); x=np.linalg.solve(A,b); assert np.linalg.norm(A@x-b)<1e-14
