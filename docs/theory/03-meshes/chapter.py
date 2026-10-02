# %% [markdown]
"""# Meshes

## Scientific objective

Topology defines points, faces, owner/neighbour and cell adjacency. Geometry derives centres, oriented areas and volumes. Internal orientation requires S_N=-S_P. A closed polyhedron must have consistent surface closure and positive physical volume. Import must preserve orientation and boundary identity or report loss explicitly.

## Executable contract

The code cell below is deliberately small: it is a deterministic algebraic sanity check. Chapter-specific quantitative experiments must be added beside this source and linked to the relevant V&V evidence.

## Status

Implementation, verification, validation and qualification are separate claims. This chapter does not promote CFDX maturity.
"""

# %%
from __future__ import annotations

import numpy as np; S=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float); assert np.allclose(S.sum(0),0)