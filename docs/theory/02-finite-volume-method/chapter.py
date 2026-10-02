# %% [markdown]
"""# Finite-Volume Method

## Scientific objective

Integrating a conservation law over cell P gives d/dt integral_VP(phi)dV + sum_f(F_f.S_f)=integral_VP(S)dV. Geometry enters through volume and oriented face vectors; reconstruction supplies face states; assembly produces A*x=b. Conservation is structural and independent of nonlinear convergence.

## Executable contract

The code cell below is deliberately small: it is a deterministic algebraic sanity check. Chapter-specific quantitative experiments must be added beside this source and linked to the relevant V&V evidence.

## Status

Implementation, verification, validation and qualification are separate claims. This chapter does not promote CFDX maturity.
"""

# %%
from __future__ import annotations

gamma=2.; area=3.; distance=4.; assert abs(gamma*area/distance-1.5)<1e-14