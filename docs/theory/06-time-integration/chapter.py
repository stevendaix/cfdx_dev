# %% [markdown]
"""# Time Integration

## Scientific objective

After spatial discretisation, M du/dt=R(u,t). Backward Euler is first order; Crank-Nicolson is second order for smooth problems under its assumptions; BDF2 uses (3u[n+1]-4u[n]+u[n-1])/(2 dt). Temporal order must be measured with spatial error controlled.

## Executable contract

This deterministic cell is a minimal mathematical sanity check. Quantitative studies must be linked to V&V evidence.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

import numpy as np; dt=.1; u=1.; rate=-1.; u=u/(1-dt*rate); assert np.isfinite(u)
