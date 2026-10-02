# %% [markdown]
"""# Foundations

## Scientific objective

A continuum field phi(x,t) is represented on a domain Omega. Conservation laws, constitutive laws, initial data and boundary data define the mathematical problem. Nondimensionalisation introduces characteristic scales and dimensionless groups such as Re=rho*U*L/mu. The numerical chain is PDE -> integral balance -> discrete operators -> algebraic system -> iterative solution -> diagnostics.

## Executable contract

The code cell below is deliberately small: it is a deterministic algebraic sanity check. Chapter-specific quantitative experiments must be added beside this source and linked to the relevant V&V evidence.

## Status

Implementation, verification, validation and qualification are separate claims. This chapter does not promote CFDX maturity.
"""

# %%
from __future__ import annotations

re=1.0*1.0*1.0/0.01; assert re==100.0