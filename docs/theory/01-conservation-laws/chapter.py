# %% [markdown]
"""# Conservation Laws

## Scientific objective

For density q, flux F and source s, d/dt integral_V(q)dV + integral_boundary(F.n)dS = integral_V(s)dV. The divergence theorem gives dt(q)+div(F)=s. Finite volume preserves the integral statement, so internal face contributions must cancel pairwise.

## Executable contract

The code cell below is deliberately small: it is a deterministic algebraic sanity check. Chapter-specific quantitative experiments must be added beside this source and linked to the relevant V&V evidence.

## Status

Implementation, verification, validation and qualification are separate claims. This chapter does not promote CFDX maturity.
"""

# %%
from __future__ import annotations

import numpy as np; S=np.array([[1.,0.,0.],[-1.,0.,0.]]); assert np.allclose(S.sum(axis=0),0.)