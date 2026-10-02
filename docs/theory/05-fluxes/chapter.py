# %% [markdown]
"""# Fluxes

## Scientific objective

Convective scalar flux is mdot_f*phi_f with mdot_f=rho_f*u_f.S_f. Diffusive flux is Gamma_f*grad(phi)_f.S_f. A scheme includes mass-flux definition, face reconstruction, limiting and boundary treatment. Internal-face antisymmetry, accuracy and boundedness are separate properties.

## Executable contract

This deterministic cell is a minimal mathematical sanity check. Quantitative studies must be linked to V&V evidence.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

mdot=2.; phi=.75; owner=mdot*phi; neighbour=-owner; assert abs(owner+neighbour)<1e-14
