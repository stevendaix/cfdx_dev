# %% [markdown]
"""# Data Model and File Formats

Scientific executable chapter. Follow physical motivation, definitions, assumptions, governing equations, discretisation, numerical properties, CFDX implementation, verification evidence, limitations and references. Implementation is not automatically verification, validation or qualification.

## Core equations and contracts

- case setup is declarative state
- HDF5 datasets map topology, fields and metadata to typed arrays
- restart state must preserve sufficient history for the selected temporal scheme
- VTU stores visualization-oriented field data rather than authoritative case setup
- provenance links result to case, mesh, numerical configuration and software revision

## Evidence rule

Every implementation claim must map to actual source files and tests. Every numerical result must identify configuration, mesh/resolution, solver settings and software revision. Missing evidence is a documented gap.
"""

# %%
from __future__ import annotations
import numpy as np
x=np.linspace(0.0,1.0,5)
assert np.all(np.isfinite(x))
payload=np.array([1.,2.,3.])
assert payload.dtype.kind == "f"
assert payload.shape == (3,)
