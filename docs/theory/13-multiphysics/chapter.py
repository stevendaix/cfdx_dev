# %% [markdown]
"""# Multiphysics

## Scientific objective

A coupled model is F(y)=0 with y containing all physical fields. Its Jacobian has blocks J_ij describing direct coupling. Segregated coupling is a fixed-point iteration y_(m+1)=G(y_m), optionally relaxed; monolithic methods solve the coupled block system. Interface equations must explicitly enforce compatibility and conservation.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

x=1.0; x=0.5*x+0.5; assert x==1.0
